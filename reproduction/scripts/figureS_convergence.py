#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent; ROOT=HERE.parent; REF=ROOT/'reference_data'; OUT=ROOT/'generated'; OUT.mkdir(exist_ok=True)
df=pd.read_csv(REF/'convergence.csv')
fig,axs=plt.subplots(1,2,figsize=(8.3,3.4)); size=df[df.test=='size']
for protocol,d in size.groupby('protocol'):
    d=d.sort_values('n_atoms'); axs[0].errorbar(d.n_atoms,d.T_eff,yerr=1.96*d.T_stderr,marker='o',capsize=2,label=protocol)
axs[0].set(xlabel='System size (atoms)',ylabel='Effective transmittance',ylim=(0,1.05)); axs[0].legend(frameon=False,fontsize=8)
stat=df[df.test=='statistics']; g=stat.groupby('n_trajectories').agg(T_mean=('T_eff','mean'),T_std=('T_eff','std')).reset_index(); axs[1].errorbar(g.n_trajectories,g.T_mean,yerr=g.T_std,marker='o',capsize=3); axs[1].set(xlabel='Number of trajectories',ylabel='Mean effective transmittance',ylim=(0.85,1.02))
for ax,l in zip(axs,['(a)','(b)']): ax.grid(alpha=.22); ax.text(.02,.96,l,transform=ax.transAxes,va='top',fontweight='bold')
fig.tight_layout(); fig.savefig(OUT/'figureS_convergence.png',dpi=300,bbox_inches='tight'); plt.close(fig)
