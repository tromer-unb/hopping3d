#!/usr/bin/env python3
"""Higher-statistics BN tensor anisotropy sweep for the real 3D crystal geometry."""
import os, math
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from hopping3d.core import load_structure,simulate_periodic_diffusion,unit_vector
from hopping3d.tensor import principal_decomposition,dual_unit_vectors,directional_value
HERE=Path(__file__).resolve().parent; ROOT=HERE.parent; OUT=ROOT/'generated'; OUT.mkdir(exist_ok=True)
full=os.environ.get('FULL','0')=='1'; atoms=load_structure(ROOT/'inputs/BN_bulk.cif',repeat=(2,2,2)); normal=unit_vector(atoms.cell[2]).tolist(); G=dual_unit_vectors(atoms.cell)
etas=[0.01,0.02,0.05,0.10,0.20,0.50,1.0]; reps=12 if full else 7; walkers=1800 if full else 900; obs=1.2e-10 if full else 8e-11
base=dict(temperature_K=300.0,gap_eV=0.0,xi0_A=2.5,alpha_gap=0.0,beta_T=0.0,nu0_Hz=1e13,gap_activation='none',pair_scales={'B-N':1.0,'B-B':0.05,'N-N':0.05},layer_normal=normal,interlayer_threshold_A=1.0,intralayer_cutoff_A=1.9,interlayer_cutoff_A=3.8,initial_distribution='equilibrium')
rows=[]
for eta in etas:
  for r in range(reps):
    out=simulate_periodic_diffusion(atoms,3.8,dict(base,interlayer_scale=eta),{},n_walkers=walkers,observation_time_s=obs,seed=1010000+10000*int(round(eta*100))+r)
    D=np.asarray(out['D_m2_s']); vals,vecs=principal_decomposition(D); dA=directional_value(D,G[:,0]); dB=directional_value(D,G[:,1]); dC=directional_value(D,G[:,2])
    rows.append(dict(interlayer_scale=eta,realization=r,D1=vals[0],D2=vals[1],D3=vals[2],D_planeA=dA,D_planeB=dB,D_planeC=dC,anisotropy_principal=.5*(vals[0]+vals[1])/vals[2],anisotropy_plane=.5*(dA+dB)/dC,mean_hops=out['mean_hops']))
raw=pd.DataFrame(rows); raw.to_csv(OUT/'BN_tensor_production_raw.csv',index=False); summary=[]
for eta,d in raw.groupby('interlayer_scale'):
  row=dict(interlayer_scale=eta,n_realizations=len(d))
  for c in ['D1','D2','D3','anisotropy_principal','anisotropy_plane']:
    x=d[c].to_numpy(); row[c+'_mean']=x.mean(); row[c+'_std']=x.std(ddof=1); row[c+'_sem']=x.std(ddof=1)/math.sqrt(len(x))
  summary.append(row)
s=pd.DataFrame(summary); s.to_csv(OUT/'BN_tensor_production_summary.csv',index=False)
fig,axes=plt.subplots(1,2,figsize=(8.2,3.5)); axes[0].errorbar(s.interlayer_scale,s.anisotropy_principal_mean,yerr=1.96*s.anisotropy_principal_sem,marker='o',capsize=2,label='principal axes'); axes[0].errorbar(s.interlayer_scale,s.anisotropy_plane_mean,yerr=1.96*s.anisotropy_plane_sem,marker='s',capsize=2,label='crystal plane normals')
for c in ['D1','D2','D3']: axes[1].errorbar(s.interlayer_scale,s[c+'_mean'],yerr=1.96*s[c+'_sem'],marker='o',capsize=2,label=c)
for ax in axes: ax.set_xscale('log'); ax.set_yscale('log'); ax.grid(alpha=.25); ax.set_xlabel(r'$\eta_\perp$'); ax.legend(frameon=False,fontsize=8)
axes[0].set_ylabel(r'$D_{in}/D_{out}$'); axes[1].set_ylabel(r'$D$ (m$^2$/s)'); fig.tight_layout(); fig.savefig(OUT/'figure_BN_tensor_production.png',dpi=300,bbox_inches='tight'); plt.close(fig); print(s.to_string(index=False))
