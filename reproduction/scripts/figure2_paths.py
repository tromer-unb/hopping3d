#!/usr/bin/env python3
"""Representative x-z trajectories showing layer switching after 3D edges are enabled."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from hopping3d.core import load_structure,apply_random_vacancies,build_finite_graph,direction_vector,electrode_indices
from hopping3d.graph import active_adjacency,reachable_from_targets
from hopping3d.survival import simulate_detailed
HERE=Path(__file__).resolve().parent; ROOT=HERE.parent; OUT=ROOT/'generated'; OUT.mkdir(exist_ok=True)
raw=pd.read_csv(OUT/'graphene_topology_raw.csv'); piv=raw[raw.vacancy_fraction==0.25].pivot(index='realization',columns='topology',values='source_reachable_fraction').dropna(); target=float(piv.strict_2D.mean()); cand=piv[piv['3D_enabled']>=1.0-1e-12].copy(); rr=int((cand.strict_2D-target).abs().idxmin()); seed=510000+10000*25+rr
base=load_structure(ROOT/'inputs/stacked_graphene_4layers.cif'); atoms=apply_random_vacancies(base,.25,seed=seed); graph=build_finite_graph(atoms,3.75)
u=np.asarray(atoms.positions)@direction_vector('X',atoms); src,umin,umax,w=electrode_indices(u,2.0,None); drain=np.where(u>=umax-w)[0]
etas=[0.0,0.01,0.10,1.0]; fig,axes=plt.subplots(1,4,figsize=(14.5,3.5),sharex=True,sharey=True); rows=[]
for ax,eta in zip(axes,etas):
  adj=active_adjacency(graph,[0,0,1],eta,1.75,3.75,1.0); reach=reachable_from_targets(adj,drain); rsrc=np.asarray([int(i) for i in src if int(i) in reach],dtype=int); preach=len(rsrc)/len(src)
  cfg=dict(direction='X',temperature_K=300.,gap_eV=0.,xi0_A=3.,alpha_gap=0.,beta_T=0.,nu0_Hz=1e13,max_time_s=5.12e-9,max_steps=50000,n_trajectories=120,electrode_width_A=2.0,gap_activation='none',pair_scales={},layer_normal=[0,0,1],interlayer_scale=eta,interlayer_threshold_A=1.,intralayer_cutoff_A=1.75,interlayer_cutoff_A=3.75,source_indices_override=rsrc.tolist())
  out=simulate_detailed(atoms,graph,.05,cfg,{},rng_seed=1200000+int(eta*10000),record_paths=18); pos=np.asarray(atoms.positions); ax.scatter(pos[:,0],pos[:,2],s=8,alpha=.18)
  for pth in out['paths']:
    ids=np.asarray(pth['indices'],int); ax.plot(pos[ids,0],pos[ids,2],lw=.75,alpha=.6)
  ax.set_title(rf'$\eta_\perp={eta:g}$'+'\n'+rf'$P_{{reach}}={preach:.2f}$'); ax.set_xlabel('x (Å)'); ax.grid(alpha=.12)
  rows.append(dict(realization=rr,seed=seed,interlayer_scale=eta,source_reachable_fraction=preach,T_conditional=out['T_eff'],reduced_tau=1e13*out['mean_first_passage_s'],mean_tortuosity=out['mean_tortuosity']))
axes[0].set_ylabel('layer coordinate z (Å)'); fig.tight_layout(); fig.savefig(OUT/'figure_graphene_paths_xz.png',dpi=300,bbox_inches='tight'); plt.close(fig); pd.DataFrame(rows).to_csv(OUT/'graphene_representative_paths_metrics.csv',index=False)
