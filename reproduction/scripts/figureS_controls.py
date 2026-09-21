#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent; ROOT=HERE.parent; REF=ROOT/'reference_data'; OUT=ROOT/'generated'; OUT.mkdir(exist_ok=True)
att=pd.read_csv(REF/'attempt_frequency_invariance_c4.csv'); eq=pd.read_csv(REF/'equilibrium_initialization_c5.csv'); gap=pd.read_csv(REF/'gap_ablation_study6.csv')
fig,axs=plt.subplots(1,3,figsize=(11.6,3.4))
axs[0].plot(att.nu0_Hz,att.tau_s,'o-',label=r'$\tau$'); ax2=axs[0].twinx(); ax2.plot(att.nu0_Hz,att.reduced_tau,'s--',label=r'$\nu_0\tau$')
axs[0].set_xscale('log'); axs[0].set_yscale('log'); axs[0].set(xlabel=r'$\nu_0$ (s$^{-1}$)',ylabel=r'$\tau$ (s)'); ax2.set_ylabel(r'$\nu_0\tau$')
for mode,d in eq.groupby('initial_distribution'): axs[1].errorbar(d.observation_time_s,d.D_parallel_mean,yerr=d.D_parallel_std,marker='o',capsize=2,label=mode)
axs[1].set_xscale('log'); axs[1].set(xlabel='observation time (s)',ylabel=r'$D_{\parallel}$ (m$^2$/s)'); axs[1].legend(frameon=False)
labels={'gap_0_no_activation':'Eg=0','gap_1.3_localization_only':'Eg=1.3, localization only','gap_1.3_legacy_global':'legacy global activation'}
for model,d in gap.groupby('model'): axs[2].plot(d.V,d.T_eff,'o-',label=labels.get(model,model))
axs[2].set(xlabel='Bias (V)',ylabel='Effective transmittance',ylim=(-.03,1.03)); axs[2].legend(frameon=False,fontsize=7)
for ax,l in zip(axs,['(a)','(b)','(c)']): ax.grid(alpha=.22); ax.text(.02,.96,l,transform=ax.transAxes,va='top',fontweight='bold')
fig.tight_layout(); fig.savefig(OUT/'figureS_controls.png',dpi=300,bbox_inches='tight'); plt.close(fig)
