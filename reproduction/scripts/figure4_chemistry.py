#!/usr/bin/env python3
"""Higher-statistics empirical chemistry map on W2O6 geometry."""
import os, math
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from hopping3d.core import load_structure,build_finite_graph,simulate_first_passage
HERE=Path(__file__).resolve().parent; ROOT=HERE.parent; OUT=ROOT/'generated'; OUT.mkdir(exist_ok=True)
full=os.environ.get('FULL','0')=='1'; atoms=load_structure(ROOT/'inputs/W2O6.cif'); graph=build_finite_graph(atoms,3.4); syms=atoms.get_chemical_symbols()
deltas=[0.0,0.25,0.50,1.0]; homos=[1.0,0.3,0.1,0.03]; reps=20 if full else 10; ntraj=1200 if full else 500; npaths=120 if full else 60
base=dict(direction='X',temperature_K=300.0,gap_eV=0.0,xi0_A=3.0,alpha_gap=0.0,beta_T=0.0,nu0_Hz=1e13,max_time_s=1e-9,max_steps=5000,n_trajectories=ntraj,electrode_width_A=2.0,gap_activation='none')
rows=[]
for de in deltas:
  chem={'site_offsets_eV':{'O':-de/2,'W':de/2}}
  for homo in homos:
    for r in range(reps):
      cfg=dict(base,pair_scales={'O-W':1.0,'O-O':homo,'W-W':homo}); out=simulate_first_passage(atoms,graph,.10,cfg,chem,rng_seed=910000+10000*int(round(de*100))+100*int(round(homo*100))+r,record_paths=npaths)
      counts={'O':0,'W':0}
      for p in out['paths']:
        for i in p['indices']: counts[syms[i]]+=1
      tot=max(sum(counts.values()),1); rows.append(dict(delta_site_eV=de,homobond_scale=homo,realization=r,T_eff=out['T_eff'],reduced_tau=1e13*out['mean_first_passage_s'],mean_tortuosity=out['mean_tortuosity'],O_visit_fraction=counts['O']/tot,W_visit_fraction=counts['W']/tot))
raw=pd.DataFrame(rows); raw.to_csv(OUT/'W2O6_chemistry_production_raw.csv',index=False); summary=[]
for (de,h),d in raw.groupby(['delta_site_eV','homobond_scale']):
  row=dict(delta_site_eV=de,homobond_scale=h,n_realizations=len(d))
  for c in ['T_eff','reduced_tau','mean_tortuosity','O_visit_fraction','W_visit_fraction']:
    x=d[c].to_numpy(float); row[c+'_mean']=x.mean(); row[c+'_std']=x.std(ddof=1); row[c+'_sem']=x.std(ddof=1)/math.sqrt(len(x))
  summary.append(row)
s=pd.DataFrame(summary); s.to_csv(OUT/'W2O6_chemistry_production_summary.csv',index=False)
fig,axes=plt.subplots(1,2,figsize=(7.8,3.5))
for h,d in s.groupby('homobond_scale'):
  axes[0].errorbar(d.delta_site_eV,d.reduced_tau_mean,yerr=1.96*d.reduced_tau_sem,marker='o',capsize=2,label=f'w_homo={h:g}'); axes[1].errorbar(d.delta_site_eV,d.W_visit_fraction_mean,yerr=1.96*d.W_visit_fraction_sem,marker='o',capsize=2,label=f'w_homo={h:g}')
axes[0].set_yscale('log'); axes[0].set(xlabel=r'$\Delta\epsilon_{W-O}$ (eV)',ylabel=r'$\nu_0\langle\tau_{FP}\rangle$'); axes[1].set(xlabel=r'$\Delta\epsilon_{W-O}$ (eV)',ylabel='W fraction on successful paths',ylim=(-.03,.55))
for ax in axes: ax.grid(alpha=.25)
axes[0].legend(frameon=False,fontsize=7); fig.tight_layout(); fig.savefig(OUT/'figure_W2O6_chemistry_production.png',dpi=300,bbox_inches='tight'); plt.close(fig); print(s.to_string(index=False))
