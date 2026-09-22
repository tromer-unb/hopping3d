"""Automatic multi-adsorbate parametrization for graphene/molecular interfaces."""
from __future__ import annotations
from collections import Counter
import numpy as np
from ase import Atoms
from ase.neighborlist import neighbor_list,natural_cutoffs
from .alignment import K_BENZENE_PI,molecule_graphene_alignment_eht
from .molecular import graphene_benzene_pi_interface,conjugated_pi_dimer_couplings
from .interface import (bonded_components,_component_record,_plane_normal,_unwrapped_component_positions,
    _interface_geometry,_local_interface_cluster,_parameter,BENZENE_REORGANIZATION_EV)

def detect_graphene_adsorbates(atoms:Atoms):
    comps=bonded_components(atoms); records=[_component_record(atoms,c) for c in comps]; host_pos=None
    for k,c in enumerate(comps):
        sy={atoms[int(i)].symbol for i in c}
        if sy=={'C'} and len(c)>=12: host_pos=k; break
    if host_pos is None: raise ValueError('No graphene-like all-carbon host detected')
    host_idx=comps[host_pos]; ads=[c for k,c in enumerate(comps) if k!=host_pos and len(c)>=2]
    if not ads: raise ValueError('No molecular adsorbates detected')
    host=atoms[host_idx]; host.set_pbc(atoms.pbc); cut=natural_cutoffs(host,mult=1.15); ii,_=neighbor_list('ij',host,cut); deg=np.bincount(ii,minlength=len(host)); _,rms=_plane_normal(host.positions)
    if np.median(deg)<2.5 or rms>=0.25: raise ValueError('Carbon host fails graphene-like checks')
    for c in ads:
        if not {atoms[int(i)].symbol for i in c}.issubset({'C','H'}): raise ValueError('Current coverage mode supports separated conjugated C/H adsorbates')
    return host_idx,ads,{'host_type':'graphene-like carbon pi sheet','n_adsorbates':len(ads),'components':records,'host_component':int(host_pos),'host_median_coordination':float(np.median(deg)),'host_plane_rms_A':float(rms)}

def _component_atoms(atoms,idx):
    a=atoms[np.asarray(idx,int)].copy(); a.positions=_unwrapped_component_positions(atoms,idx); a.set_pbc(False); return a

def _nearest_image_pair(atoms,idx_a,idx_b):
    A=_component_atoms(atoms,idx_a); B=_component_atoms(atoms,idx_b); ca=A.positions.mean(axis=0); cb=B.positions.mean(axis=0); dr=cb-ca; cell=np.asarray(atoms.cell.array,float); pbc=np.asarray(atoms.pbc,bool)
    if abs(np.linalg.det(cell))>1e-10:
        f=np.linalg.solve(cell.T,dr)
        for a in range(3):
            if pbc[a]: f[a]-=np.round(f[a])
        mic=f@cell; B.positions+=mic-dr
    return A,B

def _surface_area_nm2(atoms,host_idx):
    hn,_=_plane_normal(atoms.positions[np.asarray(host_idx,int)]); cand=[]
    for a in range(3):
        v=np.asarray(atoms.cell[a],float); L=np.linalg.norm(v)
        if L>1e-8 and atoms.pbc[a]: cand.append((abs(float(np.dot(v/L,hn))),a))
    if len(cand)<2: return None
    axes=[a for _,a in sorted(cand)[:2]]; return float(np.linalg.norm(np.cross(atoms.cell[axes[0]],atoms.cell[axes[1]]))/100.0)

def parameterize_coverage_interface(atoms:Atoms,patch_radius_A=8.0,electron_offset_eV=None):
    host_idx,ads,det=detect_graphene_adsorbates(atoms); host=atoms[host_idx].copy(); area=_surface_area_nm2(atoms,host_idx)
    result={'schema':'hopping3d-coverage-0.1','detection':det,'adsorbates':[],'molecule_pairs':[],'external_executables_required':[],
            'warnings':['Energy alignment remains screening-only unless independently calibrated.']}
    if area is not None: result['coverage']={'surface_area_nm2':area,'molecules_per_nm2':float(len(ads)/area)}
    align_ref=None; dyn_kin=[]; dyn_kout=[]
    for m,idx in enumerate(ads):
        mol=_component_atoms(atoms,idx); geom=_interface_geometry(atoms,host_idx,idx); symbols=Counter(mol.get_chemical_symbols()); formula='C'+(str(symbols['C']) if symbols['C']!=1 else '')
        if symbols.get('H',0): formula+='H'+(str(symbols['H']) if symbols['H']!=1 else '')
        rec={'id':m,'formula':formula,'geometry':geom,'parameters':{}}
        if formula=='C6H6':
            rec['type']='benzene'; combo,gi,mi,patch=_local_interface_cluster(atoms,host_idx,idx,patch_radius_A); eint=graphene_benzene_pi_interface(combo,gi,mi,K_BENZENE_PI)
            Vh=np.asarray(eint['hole_coupling_matrix_eV'],float); Ve=np.asarray(eint['electron_coupling_matrix_eV'],float); Sh=float(np.sum(Vh*Vh)); Se=float(np.sum(Ve*Ve))
            rec['parameters'].update({'lambda_hole_eV':BENZENE_REORGANIZATION_EV['hole'],'lambda_electron_eV':BENZENE_REORGANIZATION_EV['electron'],
                'hole_interface_strength_eV2':Sh,'electron_interface_strength_eV2':Se,'hole_J_eff_eV':float(np.sqrt(Sh/2.0)),'electron_J_eff_eV':float(np.sqrt(Se/2.0))})
            if align_ref is None: align_ref=molecule_graphene_alignment_eht(host,mol)
            if electron_offset_eV is not None:
                from .reservoir import graphene_spectral_strengths,electron_reservoir_rates
                spec=graphene_spectral_strengths(combo[gi],Ve,mode='coherent'); rr=electron_reservoir_rates(spec,float(electron_offset_eV),BENZENE_REORGANIZATION_EV['electron'],degeneracy=2,T_K=300.0)
                dyn_kin.append(float(rr['k_reservoir_to_molecule_s-1'])); dyn_kout.append(float(rr['k_molecule_to_reservoir_s-1']))
        else:
            rec['type']='conjugated C/H molecule (unclassified)'; result['warnings'].append(f'Adsorbate {m}: no shipped lambda/interface calibration for {formula}.')
        result['adsorbates'].append(rec)
    for i in range(len(ads)):
        for j in range(i+1,len(ads)):
            A,B=_nearest_image_pair(atoms,ads[i],ads[j]); pair=A+B; ca=[k for k,s in enumerate(A.get_chemical_symbols()) if s=='C']; cb=[k for k,s in enumerate(B.get_chemical_symbols()) if s=='C']; dcc=float(np.linalg.norm(A.positions[ca][:,None,:]-B.positions[cb][None,:,:],axis=2).min()); centers=float(np.linalg.norm(B.positions.mean(axis=0)-A.positions.mean(axis=0)))
            rec={'i':i,'j':j,'center_distance_A':centers,'min_CC_A':dcc,'hole_J_eff_eV':None,'electron_J_eff_eV':None}
            if result['adsorbates'][i].get('type')=='benzene' and result['adsorbates'][j].get('type')=='benzene':
                c=conjugated_pi_dimer_couplings(pair,range(12),range(12,24),K=K_BENZENE_PI,inter_cutoff_A=6.0); rec['hole_J_eff_eV']=float(c['J_hole_eV']); rec['electron_J_eff_eV']=float(c['J_electron_eV'])
            result['molecule_pairs'].append(rec)
    if align_ref is not None:
        result['screening_alignment']={'delta_mu_molecule_minus_graphene_eV':float(align_ref['delta_mu_molecule_minus_graphene_eV']),'hole_localization_offset_eV':float(align_ref['hole_localization_offset_eV']),'electron_localization_offset_eV':float(align_ref['electron_localization_offset_eV']),'confidence':'low-to-medium; isolated-fragment screening'}
    result['transport_template']={'host_model':'electronic_reservoir','molecular_degeneracy':2,'energy_offset_eV':None if electron_offset_eV is None else float(electron_offset_eV),'pair_lambda_rule':'lambda_MM = lambda_i + lambda_j for identical molecules','recommended_solver':'reservoir + molecular CTMC'}
    if electron_offset_eV is not None and len(dyn_kin)==len(ads) and all(a.get('type')=='benzene' for a in result['adsorbates']):
        from .coverage import molecular_marcus_rate,molecular_excursion_metrics
        n=len(ads); K=np.zeros((n,n),float); lam=BENZENE_REORGANIZATION_EV['electron']; lam_mm=2.0*lam
        for q in result['molecule_pairs']:
            i,j=int(q['i']),int(q['j']); J=q['electron_J_eff_eV']
            if J is not None and J>0: K[i,j]=K[j,i]=molecular_marcus_rate(J,0.0,lam_mm,300.0)
        met=molecular_excursion_metrics(np.asarray(dyn_kin),np.asarray(dyn_kout),K,_coverage_centers(atoms,ads)); result['electron_coverage_dynamics']={'energy_offset_eV':float(electron_offset_eV),'temperature_K':300.0,'lambda_GM_eV':lam,'lambda_MM_eV':lam_mm,'metrics':met,'approximation':'local coherent EHT reservoir patch per adsorbate + pairwise molecular Marcus network'}
    result['parameter_provenance']={'GM_interface_couplings':{'source':'direct local coherent pi-EHT interface block','confidence':'medium','recommended_use':'single-geometry interface kinetics'},'MM_pair_couplings':{'source':'direct calibrated pi-EHT frontier-subspace projection','confidence':'medium','recommended_use':'molecular-network Marcus kinetics'},'benzene_lambda':{'source':BENZENE_REORGANIZATION_EV['source'],'confidence':'medium','recommended_use':'low-cost Marcus kinetics'},'level_alignment':{'source':'portable EHT + offline PAH calibration','confidence':'low-to-medium','recommended_use':'screening unless environment-corrected'}}
    result['provenance']={'interface_couplings':'direct local coherent pi-EHT','molecule_pair_couplings':'direct pi-EHT frontier subspaces','lambda':'shipped offline benzene calibration when recognized','xTB_required':False,'scope':'semiempirical screening; multi-adsorbate interface is pairwise-factorized'}
    return result

def render_coverage_markdown(result):
    d=result['detection']; lines=['# Hopping3D EHT coverage report','',f"Detected host: **{d['host_type']}**",f"Detected adsorbates: **{d['n_adsorbates']}**",'']
    if 'coverage' in result:
        c=result['coverage']; lines += [f"Surface area: {c['surface_area_nm2']:.3f} nm²",f"Coverage: {c['molecules_per_nm2']:.3f} molecules/nm²",'']
    lines += ['## Adsorbates','','| id | type | height (Å) | tilt (deg) | J_e (meV) | J_h (meV) |','|---:|---|---:|---:|---:|---:|']
    for a in result['adsorbates']:
        p=a['parameters']; je=p.get('electron_J_eff_eV'); jh=p.get('hole_J_eff_eV'); lines.append(f"| {a['id']} | {a.get('type',a['formula'])} | {a['geometry']['adsorption_height_A']:.3f} | {a['geometry']['tilt_deg']:.2f} | {1000*je if je is not None else float('nan'):.3f} | {1000*jh if jh is not None else float('nan'):.3f} |")
    lines += ['', '## Molecular pairs','','| i-j | center distance (Å) | min C-C (Å) | J_e (meV) | J_h (meV) |','|---|---:|---:|---:|---:|']
    for p in result['molecule_pairs']:
        je=p['electron_J_eff_eV']; jh=p['hole_J_eff_eV']; lines.append(f"| {p['i']}-{p['j']} | {p['center_distance_A']:.3f} | {p['min_CC_A']:.3f} | {1000*je if je is not None else float('nan'):.3f} | {1000*jh if jh is not None else float('nan'):.3f} |")
    lines += ['', '## Transport template','',f"Host model: **{result['transport_template']['host_model']}**",f"Recommended solver: **{result['transport_template']['recommended_solver']}**"]
    lines += ['Energy offset remains unset until independently calibrated.',''] if result['transport_template']['energy_offset_eV'] is None else [f"User-supplied electron energy offset: **{result['transport_template']['energy_offset_eV']:.6g} eV**",'']
    if 'electron_coverage_dynamics' in result:
        m=result['electron_coverage_dynamics']['metrics']; lines += ['## Electron coverage dynamics','',f"- probability of at least one intermolecular hop: {m['probability_any_inter_molecular_hop']:.6f}",f"- expected intermolecular hops per excursion: {m['expected_inter_molecular_hops_per_excursion']:.6g}",f"- probability of returning through another molecule: {m['probability_return_via_different_molecule']:.6f}",f"- RMS lateral molecular migration: {m['rms_lateral_migration_A']:.3f} Å",'']
    return '\n'.join(lines)+'\n'

def _coverage_centers(atoms,ads):
    raw=[_component_atoms(atoms,idx).positions.mean(axis=0) for idx in ads]
    if not raw: return np.empty((0,3))
    out=[np.asarray(raw[0],float)]; cell=np.asarray(atoms.cell.array,float); pbc=np.asarray(atoms.pbc,bool); invertible=abs(np.linalg.det(cell))>1e-10
    for c in raw[1:]:
        dr=np.asarray(c,float)-out[0]
        if invertible:
            f=np.linalg.solve(cell.T,dr)
            for a in range(3):
                if pbc[a]: f[a]-=np.round(f[a])
            dr=f@cell
        out.append(out[0]+dr)
    return np.asarray(out,float)
