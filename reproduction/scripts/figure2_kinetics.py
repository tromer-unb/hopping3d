#!/usr/bin/env python3
"""Conditional first-passage kinetics at p_v=0.25 after topology is separated from kinetics."""
import os
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from hopping3d.core import load_structure,apply_random_vacancies,build_finite_graph,direction_vector,electrode_indices
from hopping3d.graph import active_adjacency,reachable_from_targets
from hopping3d.survival import simulate_detailed
HERE=Path(__file__).resolve().parent; ROOT=HERE.parent; OUT=ROOT/'generated'; OUT.mkdir(exist_ok=True)
full=os.environ.get('FULL','0')=='1'; pv=.25; etas=[.01,.02,.05,.10,.20,1.0]; reps=20 if full else 10; ntraj=800 if full else 300
base=load_structure(ROOT/'inputs/stacked_graphene_4layers.cif'); rows=[]
for r in range(reps):
  atoms=apply_random_vacancies(base,pv,seed=760000+r); graph=build_finite_graph(atoms,3.75)
  u=np.asarray(atoms.positions)@direction_vector('X',atoms); src,umin,umax,w=electrode_indices(u,2.0,None); drain=np.where(u>=umax-w)[0]
  for eta in etas:
    adj=active_adjacency(graph,[0,0,1],eta,1.75,3.75,1.0); reach=reachable_from_targets(adj,drain)
    rsrc=np.asarray([int(i) for i in src if int(i) in reach],dtype=int); preach=len(rsrc)/max(len(src),1)
    if not len(rsrc): continue
    cfg=dict(direction='X',temperature_K=300.,gap_eV=0.,xi0_A=3.,alpha_gap=0.,beta_T=0.,nu0_Hz=1e13,max_time_s=5.12e-9,max_steps=60000,n_trajectories=ntraj,electrode_width_A=2.,gap_activation='none',pair_scales={},layer_normal=[0,0,1],interlayer_scale=eta,interlayer_threshold_A=1.,intralayer_cutoff_A=1.75,interlayer_cutoff_A=3.75,source_indices_override=rsrc.tolist())
    out=simulate_detailed(atoms,graph,.05,cfg,{},rng_seed=1300000+100*r+int(round(eta*10000)),record_paths=0)
    rows.append(dict(realization=r,interlayer_scale=eta,source_reachable_fraction=preach,T_conditional=out['T_eff'],reduced_tau=1e13*out['mean_first_passage_s'],mean_tortuosity=out['mean_tortuosity'],censored_fraction=out['censored_fraction']))
raw=pd.DataFrame(rows); raw.to_csv(OUT/'graphene_kinetic_selected_raw.csv',index=False); qs=[]
for eta,d in raw.groupby('interlayer_scale'):
  q=lambda c,p: float(np.nanquantile(d[c],p))
  qs.append(dict(interlayer_scale=eta,n_realizations=len(d),preach_mean=d.source_reachable_fraction.mean(),T_conditional_mean=d.T_conditional.mean(),censored_mean=d.censored_fraction.mean(),tau_median=q('reduced_tau',.5),tau_q16=q('reduced_tau',.16),tau_q84=q('reduced_tau',.84),tort_median=q('mean_tortuosity',.5),tort_q16=q('mean_tortuosity',.16),tort_q84=q('mean_tortuosity',.84)))
s=pd.DataFrame(qs); s.to_csv(OUT/'graphene_kinetic_selected_summary.csv',index=False)
fig,axes=plt.subplots(1,2,figsize=(7.5,3.4)); x=s.interlayer_scale.to_numpy(); axes[0].plot(x,s.tau_median,'o-'); axes[0].fill_between(x,s.tau_q16,s.tau_q84,alpha=.2); axes[1].plot(x,s.tort_median,'o-'); axes[1].fill_between(x,s.tort_q16,s.tort_q84,alpha=.2)
for ax in axes: ax.set_xscale('log'); ax.grid(alpha=.25); ax.set_xlabel(r'$\eta_\perp$')
axes[0].set_yscale('log'); axes[0].set_ylabel(r'median $\nu_0\langle\tau_{FP}\rangle$'); axes[1].set_ylabel('median tortuosity')
fig.tight_layout(); fig.savefig(OUT/'figure_graphene_kinetic_selected.png',dpi=300,bbox_inches='tight'); plt.close(fig); print(s.to_string(index=False))
