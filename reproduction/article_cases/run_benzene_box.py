#!/usr/bin/env python3
"""Reproduce the 100-benzene periodic-box EHT + Marcus benchmark."""
from pathlib import Path
import json, os
import numpy as np
from ase import Atoms
from ase.io import write
from hopping3d_eht.molecular import benzene_monomer,conjugated_pi_dimer_couplings
from hopping3d.core import build_periodic_graph,simulate_periodic_diffusion
from hopping3d.tensor import principal_decomposition

OUT=Path(__file__).resolve().parent/'results'; OUT.mkdir(exist_ok=True)
K=2.417879719037213; A=7.0; SHAPE=(5,5,4); T=300.; LAMBDA=.20
FULL=os.environ.get('FULL','0')=='1'; NW=1400 if FULL else 300; SEEDS=[11,22,33,44] if FULL else [11]; NSTRUCT=10 if FULL else 2
base=benzene_monomer()

def random_rotation(rng):
    u1,u2,u3=rng.random(3); q=np.array([np.sqrt(1-u1)*np.sin(2*np.pi*u2),np.sqrt(1-u1)*np.cos(2*np.pi*u2),np.sqrt(u1)*np.sin(2*np.pi*u3),np.sqrt(u1)*np.cos(2*np.pi*u3)])
    x,y,z,w=q
    return np.array([[1-2*(y*y+z*z),2*(x*y-z*w),2*(x*z+y*w)],[2*(x*y+z*w),1-2*(x*x+z*z),2*(y*z-x*w)],[2*(x*z-y*w),2*(y*z+x*w),1-2*(x*x+y*y)]])

def build_molecules(mode,seed=1234):
    rng=np.random.default_rng(seed); mols=[]; centers=[]; atoms=Atoms(cell=np.array(SHAPE)*A,pbc=True)
    for ix in range(SHAPE[0]):
      for iy in range(SHAPE[1]):
       for iz in range(SHAPE[2]):
        m=base.copy(); m.positions=m.positions@(np.eye(3) if mode=='aligned' else random_rotation(rng)).T
        c=np.array([ix*A,iy*A,iz*A]); m.positions+=c; mols.append(m); centers.append(c); atoms+=m
    atoms.set_cell(np.array(SHAPE)*A); atoms.pbc=True
    return mols,np.asarray(centers),atoms

def coarse_graph(centers):
    s=Atoms('C'*len(centers),positions=centers,cell=np.array(SHAPE)*A,pbc=True); return s,build_periodic_graph(s,A+0.05)

def pair_J(mols,graph):
    J={}; vals=[]; seen=set(); minsep=1e9
    for i,edges in enumerate(graph):
      for e in edges:
        key=(min(i,e.j),max(i,e.j))
        if key in seen: continue
        seen.add(key); j=e.j; A0=mols[i].copy(); B0=mols[j].copy(); A0.positions-=A0.positions.mean(0); B0.positions-=B0.positions.mean(0); B0.positions+=e.dr
        d=A0+B0; sep=float(np.linalg.norm(A0.positions[:,None,:]-B0.positions[None,:,:],axis=2).min()); minsep=min(minsep,sep)
        q=conjugated_pi_dimer_couplings(d,range(12),range(12,24),K=K,inter_cutoff_A=8.0); k=f'{key[0]}-{key[1]}'; J[k]=float(q['J_hole_eV']); vals.append([i,j,*e.dr,J[k],q['J_electron_eV']])
    return J,np.asarray(vals,float),minsep

def transport(sites,Jmap):
    cfg={'temperature_K':T,'rate_model':'marcus','reorganization_energy_eV':LAMBDA,'edge_transfer_integrals_eV':Jmap,'strict_edge_transfer_integrals':True,'initial_distribution':'equilibrium'}; chem={'site_offsets_array_eV':[0.0]*len(sites)}
    tobs=1e-7; pilot=simulate_periodic_diffusion(sites,A+0.05,cfg,chem,n_walkers=150,observation_time_s=tobs,seed=777)
    while pilot['mean_hops']<100 and tobs<1.0: tobs*=10; pilot=simulate_periodic_diffusion(sites,A+0.05,cfg,chem,n_walkers=150,observation_time_s=tobs,seed=777)
    Ds=[]; hops=[]
    for seed in SEEDS:
        o=simulate_periodic_diffusion(sites,A+0.05,cfg,chem,n_walkers=NW,observation_time_s=tobs,seed=seed); Ds.append(o['D_m2_s']); hops.append(o['mean_hops'])
    D=np.mean(np.asarray(Ds),axis=0); v,_=principal_decomposition(D); return D,np.asarray(v),float(np.mean(hops)),tobs

report={'method':'100-molecule periodic benzene box; calibrated pi-EHT pair couplings + Marcus KMC','shape':SHAPE,'n_molecules':100,'spacing_A':A,'lambda_eV':LAMBDA,'lambda_role':'controlled external Marcus parameter, not supplied by EHT'}
mols,c,atom=build_molecules('aligned'); sites,g=coarse_graph(c); J,tab,minsep=pair_J(mols,g); Dt,D,h,t=transport(sites,J)
write(OUT/'benzene_box_100_aligned.cif',atom); np.savetxt(OUT/'pair_couplings_aligned.csv',tab,delimiter=',',header='i,j,dx_A,dy_A,dz_A,J_hole_eV,J_electron_eV',comments='')
report['aligned']={'n_edges':len(J),'min_inter_molecular_atom_distance_A':minsep,'J_hole_median_meV':1000*float(np.median(list(J.values()))),'D_principal_m2_s':D.tolist(),'anisotropy_D1_D3':float(D[0]/D[-1]),'mean_hops':h}
DT=[]; medJ=[]
for rr in range(NSTRUCT):
    mols,c,atom=build_molecules('isotropic',1234+rr); sites,g=coarse_graph(c); J,tab,_=pair_J(mols,g); dt,dv,h,t=transport(sites,J); DT.append(dt); medJ.append(1000*np.median(list(J.values())))
    if rr==0: write(OUT/'benzene_box_100_isotropic.cif',atom); np.savetxt(OUT/'pair_couplings_isotropic.csv',tab,delimiter=',',header='i,j,dx_A,dy_A,dz_A,J_hole_eV,J_electron_eV',comments='')
vv,_=principal_decomposition(np.mean(DT,axis=0)); report['isotropic']={'n_structural_realizations':NSTRUCT,'n_edges_per_realization':300,'J_hole_median_meV_mean':float(np.mean(medJ)),'ensemble_D_principal_m2_s':vv.tolist(),'ensemble_anisotropy_D1_D3':float(vv[0]/vv[-1])}
(OUT/('benzene_box_summary_full.json' if FULL else 'benzene_box_summary.json')).write_text(json.dumps(report,indent=2)+'\n'); print(json.dumps(report,indent=2)); print('BENZENE BOX 100 PASSED')
