#!/usr/bin/env python3
"""Controlled four-layer carbon topology statistics for dimensional rescue."""
import os, math
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from hopping3d.core import load_structure,apply_random_vacancies,build_finite_graph,direction_vector,electrode_indices
from hopping3d.graph import active_adjacency,reachable_from_targets
HERE=Path(__file__).resolve().parent; ROOT=HERE.parent; OUT=ROOT/'generated'; OUT.mkdir(exist_ok=True)
full=os.environ.get('FULL','0')=='1'; base=load_structure(ROOT/'inputs/stacked_graphene_4layers.cif')
pvals=[0.0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40]; reps=600 if full else 250; rows=[]
for pv in pvals:
  for r in range(reps):
    atoms=apply_random_vacancies(base,pv,seed=510000+10000*int(round(pv*100))+r); graph=build_finite_graph(atoms,3.75)
    u=np.asarray(atoms.positions)@direction_vector('X',atoms); src,umin,umax,w=electrode_indices(u,2.0,None); drain=np.where(u>=umax-w)[0]
    for label,eta in [('strict_2D',0.0),('3D_enabled',1.0)]:
      adj=active_adjacency(graph,[0,0,1],eta,1.75,3.75,1.0); reach=reachable_from_targets(adj,drain)
      frac=float(np.mean([int(i) in reach for i in src])) if len(src) else 0.0
      rows.append(dict(vacancy_fraction=pv,realization=r,topology=label,n_atoms=len(atoms),n_source=len(src),n_drain=len(drain),source_reachable_fraction=frac,all_nodes_reachable_fraction=len(reach)/max(len(atoms),1),fully_source_connected=float(frac>=1-1e-12)))
raw=pd.DataFrame(rows); raw.to_csv(OUT/'graphene_topology_raw.csv',index=False); summary=[]
for (pv,top),d in raw.groupby(['vacancy_fraction','topology']):
  x=d.source_reachable_fraction.to_numpy(); f=d.fully_source_connected.to_numpy(); n=len(d)
  summary.append(dict(vacancy_fraction=pv,topology=top,n_realizations=n,source_reach_mean=x.mean(),source_reach_std=x.std(ddof=1),source_reach_sem=x.std(ddof=1)/math.sqrt(n),full_connect_probability=f.mean(),full_connect_sem=f.std(ddof=1)/math.sqrt(n)))
s=pd.DataFrame(summary); s.to_csv(OUT/'graphene_topology_summary.csv',index=False)
fig,axes=plt.subplots(1,2,figsize=(8.3,3.5))
for top,d in s.groupby('topology'):
  axes[0].errorbar(d.vacancy_fraction,d.source_reach_mean,yerr=1.96*d.source_reach_sem,marker='o',capsize=2,label=top)
  axes[1].errorbar(d.vacancy_fraction,d.full_connect_probability,yerr=1.96*d.full_connect_sem,marker='o',capsize=2,label=top)
axes[0].set(xlabel='Vacancy fraction',ylabel='Reachable source fraction',ylim=(-.03,1.03)); axes[1].set(xlabel='Vacancy fraction',ylabel='P(all source sites reach drain)',ylim=(-.03,1.03))
for ax in axes: ax.grid(alpha=.25); ax.legend(frameon=False,fontsize=8)
fig.tight_layout(); fig.savefig(OUT/'figure_graphene_topology_ensemble.png',dpi=300,bbox_inches='tight'); plt.close(fig)
print(s.to_string(index=False))
