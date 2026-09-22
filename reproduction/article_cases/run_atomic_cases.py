#!/usr/bin/env python3
"""Reproduce direct-EHT numerical datasets for the atomic article cases."""
from pathlib import Path
import json, os, runpy
import numpy as np
import pandas as pd
from ase.io import read

from hopping3d_eht.eht import graphene_pi_eht
from hopping3d_eht.article_materials import bn_pi_eht, w2o6_atom_site_eht
from hopping3d.core import build_finite_graph, simulate_periodic_diffusion, simulate_first_passage
from hopping3d.tensor import principal_decomposition

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent/'results'; OUT.mkdir(exist_ok=True)
FULL=os.environ.get('FULL','0')=='1'
REPS=10 if FULL else 2
BN_WALKERS=1200 if FULL else 250
WO_TRAJ=1200 if FULL else 250

def carbon_structure():
    cif=ROOT/'reproduction/inputs/stacked_graphene_4layers.cif'
    if not cif.exists():
        runpy.run_path(str(ROOT/'reproduction/inputs/build_layered_carbon.py'),run_name='__main__')
    return read(cif)

def carbon_stack_point():
    a=carbon_structure(); e=graphene_pi_eht(a,cutoff_A=3.75); H=e['Heff']; intra=[]; inter=[]
    for i in range(len(a)):
        for j in range(i+1,len(a)):
            v=np.asarray(a.get_distance(i,j,mic=False,vector=True)); r=float(np.linalg.norm(v)); J=abs(float(H[i,j]))
            if r<=1.75: intra.append(J)
            elif r<=3.75 and abs(float(v[2]))>1.0: inter.append(J)
    ji=float(np.median(intra)); jp=float(np.median(inter))
    return {'n_atoms':len(a),'J_intralayer_median_eV':ji,'J_interlayer_median_eV':jp,
            'J_ratio':jp/ji,'EHT_interlayer_rate_scale':(jp/ji)**2,
            'interpretation':'topological rescue is unchanged for any positive interlayer edge; EHT fixes a low kinetic scale'}

def bn_point():
    a=read(ROOT/'examples/article_systems/bn/BN_bulk.cif').repeat((2,2,2)); p=bn_pi_eht(a); rows=[]; Ds=[]
    cfg={'temperature_K':300.,'rate_model':'miller_abrahams','nu0_Hz':1e13,'xi0_A':1e9,
         'edge_transfer_integrals_eV':p['edge_transfer_integrals_eV'],'J_ref_eV':p['J_ref_eV'],
         'strict_edge_transfer_integrals':True,'initial_distribution':'equilibrium'}
    chem={'site_offsets_array_eV':[0.0]*len(a)}
    for r in range(REPS):
        o=simulate_periodic_diffusion(a,3.8,cfg,chem,n_walkers=BN_WALKERS,observation_time_s=5e-11,seed=2100+r)
        D=np.asarray(o['D_m2_s']); vals,_=principal_decomposition(D); Ds.append(D)
        rows.append({'realization':r,'mean_hops':o['mean_hops'],'D1':vals[0],'D2':vals[1],'D3':vals[2],'anisotropy_D1_D3':vals[0]/vals[2]})
    pd.DataFrame(rows).to_csv(OUT/'BN_EHT_tensor_raw.csv',index=False)
    vals,_=principal_decomposition(np.mean(Ds,axis=0))
    return {'n_atoms':len(a),'J_ref_intralayer_eV':p['J_ref_eV'],'median_interlayer_J_eV':p['median_interlayer_J_eV'],
            'EHT_interlayer_rate_scale':p['median_interlayer_rate_scale'],'raw_B_minus_N_level_eV':p['raw_B_minus_N_level_eV'],
            'site_energy_use':'common site energy retained for the geometric tensor benchmark',
            'ensemble_D_principal_m2_s':[float(x) for x in vals],'ensemble_anisotropy_D1_D3':float(vals[0]/vals[2]),
            'mean_hops':float(np.mean([x['mean_hops'] for x in rows]))}

def w2o6_point():
    a=read(ROOT/'examples/article_systems/w2o6/W2O6.cif'); p=w2o6_atom_site_eht(a); g=build_finite_graph(a,3.4); sy=np.asarray(a.get_chemical_symbols())
    base={'direction':'X','temperature_K':300.,'rate_model':'miller_abrahams','nu0_Hz':1e13,'xi0_A':1e9,
          'edge_transfer_integrals_eV':p['edge_transfer_integrals_eV'],'J_ref_eV':p['J_ref_WO_eV'],
          'strict_edge_transfer_integrals':True,'n_trajectories':WO_TRAJ,'max_time_s':1e-5,'max_steps':10000,
          'electrode_width_A':2.0,'gap_activation':'none'}
    raw=[]
    for mode in ('zero_site_contrast','EHT_raw_site_contrast'):
        offsets=np.zeros(len(a)) if mode=='zero_site_contrast' else np.asarray([p['screening_site_offsets_eV'][s] for s in sy],float)
        for r in range(REPS):
            o=simulate_first_passage(a,g,.10,base,{'site_offsets_array_eV':offsets.tolist()},rng_seed=6100+100*r+(0 if mode.startswith('zero') else 1),record_paths=100)
            c={'O':0,'W':0}
            for q in o['paths']:
                for i in q['indices']: c[str(sy[i])]+=1
            tot=max(c['O']+c['W'],1)
            raw.append({'mode':mode,'realization':r,'T_eff':o['T_eff'],'reduced_tau':1e13*o['mean_first_passage_s'],
                        'mean_hops':o['mean_hops'],'mean_tortuosity':o['mean_tortuosity'],'W_visit_fraction':c['W']/tot})
    df=pd.DataFrame(raw); df.to_csv(OUT/'W2O6_EHT_transport_raw.csv',index=False); summary={}
    for mode,d in df.groupby('mode'):
        summary[mode]={c+'_mean':float(d[c].mean()) for c in ['T_eff','reduced_tau','mean_hops','mean_tortuosity','W_visit_fraction']}
    return {'edge_counts':p['edge_counts'],'J_ref_WO_eV':p['J_ref_WO_eV'],'median_OO_rate_weight':p['median_OO_rate_weight'],
            'OO_rate_weight_range':p['OO_rate_weight_range'],'raw_W5d_minus_O2p_eV':p['raw_W5d_minus_O2p_eV'],
            'transport':summary,'site_energy_note':p['transport_note'],'pair_note':p['pair_note']}

report={'method':'direct EHT parametrization of the three atomic article benchmarks','electronic_backend':'direct Extended-Huckel only',
        'carbon_stack':carbon_stack_point(),'BN_tensor':bn_point(),'W2O6_chemistry':w2o6_point(),
        'scope':'semiempirical parameter screening for transport mechanisms; not first-principles material prediction'}
(OUT/'atomic_summary.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
assert report['carbon_stack']['EHT_interlayer_rate_scale']>0
assert report['BN_tensor']['ensemble_anisotropy_D1_D3']>1
assert report['W2O6_chemistry']['edge_counts']['W-W']==0
print('Atomic article cases PASSED')
