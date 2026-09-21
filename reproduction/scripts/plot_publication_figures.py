from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent.parent; GEN=ROOT/'generated'; REF=ROOT/'reference_data'; GEN.mkdir(exist_ok=True)
def data(name):
    p=GEN/name
    return pd.read_csv(p if p.exists() else REF/name)

top=data('graphene_topology_summary.csv'); kin=data('graphene_kinetic_selected_summary.csv')
fig,ax=plt.subplots(2,2,figsize=(10,8))
for label,d in top.groupby('topology'):
    ax[0,0].errorbar(d.vacancy_fraction,d.source_reach_mean,yerr=1.96*d.source_reach_sem,marker='o',label=label)
    ax[0,1].errorbar(d.vacancy_fraction,d.full_connect_probability,yerr=1.96*d.full_connect_sem,marker='o',label=label)
ax[1,0].plot(kin.interlayer_scale,kin.tau_median,'o-'); ax[1,0].fill_between(kin.interlayer_scale,kin.tau_q16,kin.tau_q84,alpha=.2)
ax[1,1].plot(kin.interlayer_scale,kin.tort_median,'o-'); ax[1,1].fill_between(kin.interlayer_scale,kin.tort_q16,kin.tort_q84,alpha=.2)
for a in ax.flat:a.grid(alpha=.25)
ax[0,0].legend(frameon=False); ax[0,1].legend(frameon=False); ax[1,0].set_xscale('log'); ax[1,0].set_yscale('log'); ax[1,1].set_xscale('log')
ax[0,0].set(xlabel='Vacancy fraction',ylabel='Reachable source fraction'); ax[0,1].set(xlabel='Vacancy fraction',ylabel='P(all source sites reach drain)'); ax[1,0].set(xlabel=r'$\eta_{inter}$',ylabel=r'$\nu_0\tau_{FP}$'); ax[1,1].set(xlabel=r'$\eta_{inter}$',ylabel='Tortuosity')
fig.tight_layout(); fig.savefig(GEN/'figure2_dimensional_rescue.png',dpi=300); plt.close(fig)

bn=data('BN_tensor_production_summary.csv'); fig,ax=plt.subplots(1,2,figsize=(9,3.6))
for c in ['D1','D2','D3']: ax[0].errorbar(bn.interlayer_scale,bn[c+'_mean'],yerr=1.96*bn[c+'_sem'],marker='o',label=c)
ax[1].errorbar(bn.interlayer_scale,bn.anisotropy_principal_mean,yerr=1.96*bn.anisotropy_principal_sem,marker='o',label='principal'); ax[1].errorbar(bn.interlayer_scale,bn.anisotropy_plane_mean,yerr=1.96*bn.anisotropy_plane_sem,marker='s',label='crystal planes')
for a in ax:a.set_xscale('log');a.grid(alpha=.25);a.legend(frameon=False)
ax[0].set_yscale('log'); ax[0].set(xlabel=r'$\eta_{inter}$',ylabel=r'$D$ (m$^2$/s)'); ax[1].set(xlabel=r'$\eta_{inter}$',ylabel='anisotropy')
fig.tight_layout(); fig.savefig(GEN/'figure3_bn_tensor.png',dpi=300); plt.close(fig)

ch=data('W2O6_chemistry_production_summary.csv'); fig,ax=plt.subplots(1,2,figsize=(9,3.6))
for h,d in ch.groupby('homobond_scale'):
    ax[0].errorbar(d.delta_site_eV,d.reduced_tau_mean,yerr=1.96*d.reduced_tau_sem,marker='o',label=f'w={h:g}'); ax[1].errorbar(d.delta_site_eV,d.W_visit_fraction_mean,yerr=1.96*d.W_visit_fraction_sem,marker='o',label=f'w={h:g}')
ax[0].set_yscale('log'); ax[0].set(xlabel=r'$\Delta\epsilon_{W-O}$ (eV)',ylabel=r'$\nu_0\tau_{FP}$'); ax[1].set(xlabel=r'$\Delta\epsilon_{W-O}$ (eV)',ylabel='W-site visit fraction')
for a in ax:a.grid(alpha=.25)
ax[0].legend(frameon=False,fontsize=8); fig.tight_layout(); fig.savefig(GEN/'figure4_chemistry.png',dpi=300); plt.close(fig)
print('Wrote publication figures to',GEN)
