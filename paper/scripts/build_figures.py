from pathlib import Path
import json, numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.patches import Circle, RegularPolygon, FancyArrowPatch
from mpl_toolkits.mplot3d import Axes3D  # noqa
H=Path(__file__).resolve().parents[1]; D=H/'data'; F=H/'figures'; F.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.labelsize':10,'legend.fontsize':8,'savefig.dpi':600,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
def save(fig,n):
    fig.savefig(F/f'{n}.png',dpi=600,bbox_inches='tight'); fig.savefig(F/f'{n}.svg',bbox_inches='tight'); plt.close(fig)
def lab(ax,s):
    (ax.text2D if hasattr(ax,'text2D') else ax.text)(.02,.98,s,transform=ax.transAxes,va='top',ha='left',fontweight='bold',fontsize=11)
S=json.load(open(D/'EHT_article_summary.json'))
# Fig1: compact workflow, no prose panels
fig,axs=plt.subplots(1,3,figsize=(11.2,3.2))
a=axs[0]; a.axis('off'); lab(a,'(a)')
xs=[.10,.36,.62,.88]; names=['structure','EHT','rates','transport']
for x,n in zip(xs,names):
    a.add_patch(Circle((x,.55),.075,fill=False,lw=1.5)); a.text(x,.33,n,ha='center')
for x1,x2 in zip(xs[:-1],xs[1:]): a.add_patch(FancyArrowPatch((x1+.08,.55),(x2-.08,.55),arrowstyle='->',mutation_scale=12))
a.set_xlim(0,1); a.set_ylim(0,1)
a=axs[1]; a.axis('off'); lab(a,'(b)')
# atomic, molecular, hybrid icons
for x,y in [(.12,.64),(.21,.72),(.30,.62),(.21,.52)]: a.add_patch(Circle((x,y),.022,fill=False,lw=1.2))
a.plot([.12,.21,.30,.21,.12],[.64,.72,.62,.52,.64],lw=1)
for x in [.56,.68,.80]: a.add_patch(RegularPolygon((x,.63),6,radius=.05,fill=False,lw=1.2))
a.plot([.48,.90],[.31,.31],lw=2); a.add_patch(RegularPolygon((.69,.43),6,radius=.055,fill=False,lw=1.2))
a.text(.21,.22,'atomic',ha='center'); a.text(.68,.22,'molecular / hybrid',ha='center'); a.set_xlim(0,1); a.set_ylim(0,1)
a=axs[2]; lab(a,'(c)')
a.set_xscale('log'); a.set_xlim(1e-4,1); a.set_ylim(0,1); a.set_yticks([]); a.set_xlabel('model parameter')
a.axvspan(1e-4,2e-3,alpha=.12); a.axvspan(2e-3,5e-2,alpha=.12); a.axvspan(5e-2,1,alpha=.12)
a.scatter([1.1e-4,2.9e-4,3.7e-3],[.75,.5,.25],marker='*',s=90)
a.set_ylabel('material-informed anchors')
fig.tight_layout(); save(fig,'Fig1_framework_eht')
# Reuse compact plot logic for Figs 2-4 from prior script
top=pd.read_csv(D/'graphene_topology_summary.csv'); kin=pd.read_csv(D/'graphene_kinetic_selected_summary.csv'); eta=S['carbon_stack']['EHT_interlayer_rate_scale']
fig,ax=plt.subplots(2,2,figsize=(9.3,7.0))
for name,g in top.groupby('topology'):
    L='3D' if name=='3D_enabled' else '2D'; ax[0,0].errorbar(g.vacancy_fraction,g.source_reach_mean,yerr=1.96*g.source_reach_sem,marker='o',capsize=2,label=L); ax[0,1].errorbar(g.vacancy_fraction,g.full_connect_probability,yerr=1.96*g.full_connect_sem,marker='o',capsize=2,label=L)
ax[0,0].set(xlabel='Vacancy fraction',ylabel='Reachable source fraction',ylim=(-.03,1.03)); ax[0,1].set(xlabel='Vacancy fraction',ylabel='P(all source sites connected)',ylim=(-.03,1.03)); ax[0,0].legend(frameon=False); ax[0,1].legend(frameon=False)
x=kin.interlayer_scale.to_numpy(); ax[1,0].plot(x,kin.tau_median,'o-'); ax[1,0].fill_between(x,kin.tau_q16,kin.tau_q84,alpha=.18); ax[1,1].plot(x,kin.tort_median,'o-'); ax[1,1].fill_between(x,kin.tort_q16,kin.tort_q84,alpha=.18)
for a in ax[1]: a.set_xscale('log'); a.axvline(eta,ls='--',lw=1.2); a.set_xlim(5e-5,1.5)
ax[1,0].set_yscale('log'); ax[1,0].set(xlabel=r'$\eta_\perp$',ylabel=r'median $\nu_0\tau_{FP}$'); ax[1,1].set(xlabel=r'$\eta_\perp$',ylabel='Median tortuosity')
for a,s in zip(ax.flat,['(a)','(b)','(c)','(d)']): lab(a,s); a.grid(alpha=.18)
fig.tight_layout(); save(fig,'Fig2_dimensional_rescue_eht')
bn=pd.read_csv(D/'BN_empirical_sweep.csv'); e=S['BN_tensor']['EHT_interlayer_rate_scale']; Ds=np.array(S['BN_tensor']['ensemble_D_principal_m2_s'])
fig,ax=plt.subplots(2,2,figsize=(9.3,7.0))
for c,m in [('D1','o'),('D2','s'),('D3','^')]: ax[0,0].errorbar(bn.interlayer_scale,bn[c+'_mean'],yerr=bn[c+'_sem'],marker=m,capsize=2,label=c)
for v in Ds: ax[0,0].scatter([e],[v],marker='*',s=95,zorder=5)
ax[0,0].set(xscale='log',yscale='log',xlabel=r'$\eta_\perp$',ylabel=r'$D$ (m$^2$s$^{-1}$)'); ax[0,0].set_xlim(1e-4,1.5); ax[0,0].legend(frameon=False)
ax[0,1].errorbar(bn.interlayer_scale,bn.anisotropy_principal_mean,yerr=bn.anisotropy_principal_sem,marker='o',capsize=2); ax[0,1].scatter([e],[S['BN_tensor']['ensemble_anisotropy_D1_D3']],marker='*',s=110)
ax[0,1].set(xscale='log',yscale='log',xlabel=r'$\eta_\perp$',ylabel=r'$D_1/D_3$'); ax[0,1].set_xlim(1e-4,1.5)
ax[1,0].bar(['intra','inter'],[S['BN_tensor']['J_ref_intralayer_eV'],S['BN_tensor']['median_interlayer_J_eV']]); ax[1,0].set_yscale('log'); ax[1,0].set_ylabel('|J| (eV)')
ax[1,1].bar(['$D_1$','$D_2$','$D_3$'],Ds); ax[1,1].set_yscale('log'); ax[1,1].set_ylabel(r'$D$ (m$^2$s$^{-1}$)')
for a,s in zip(ax.flat,['(a)','(b)','(c)','(d)']): lab(a,s); a.grid(alpha=.18)
fig.tight_layout(); save(fig,'Fig3_BN_tensor_eht')
w=pd.read_csv(D/'W2O6_empirical_map.csv'); ew=S['W2O6_chemistry']
fig,ax=plt.subplots(2,2,figsize=(9.3,7.0))
for h,g in w.groupby('homobond_scale'):
    ax[0,0].plot(g.delta_site_eV,g.reduced_tau_mean,'o-',label=rf'{h:g}'); ax[0,1].plot(g.delta_site_eV,g.W_visit_fraction_mean,'o-')
ax[0,0].set(yscale='log',xlabel=r'$\Delta\epsilon$ (eV)',ylabel=r'$\nu_0\langle\tau_{FP}\rangle$'); ax[0,0].legend(title=r'$w_{OO}$',frameon=False,ncol=2)
ax[0,1].set(xlabel=r'$\Delta\epsilon$ (eV)',ylabel='W visit fraction',ylim=(0,.52))
ax[1,0].scatter(w.delta_site_eV,w.homobond_scale,s=25); ax[1,0].scatter([ew['raw_W5d_minus_O2p_eV']],[ew['median_OO_rate_weight']],marker='*',s=150); ax[1,0].set(yscale='log',xlabel=r'$\Delta\epsilon$ (eV)',ylabel=r'$w_{OO}$'); ax[1,0].set_xlim(-.1,4.7)
vals=[ew['transport']['zero_site_contrast']['reduced_tau_mean'],ew['transport']['EHT_raw_site_contrast']['reduced_tau_mean']]; ax[1,1].bar(['zero','EHT'],vals); ax[1,1].set_yscale('log'); ax[1,1].set_ylabel(r'$\nu_0\langle\tau_{FP}\rangle$'); a2=ax[1,1].twinx(); a2.plot([0,1],[ew['transport']['zero_site_contrast']['W_visit_fraction_mean'],ew['transport']['EHT_raw_site_contrast']['W_visit_fraction_mean']],'D--'); a2.set_ylabel('W visit fraction'); a2.set_ylim(0,.55)
for a,s in zip(ax.flat,['(a)','(b)','(c)','(d)']): lab(a,s); a.grid(alpha=.18)
fig.tight_layout(); save(fig,'Fig4_W2O6_chemistry_eht')
# Fig5: explicit 100-molecule periodic box
B=json.load(open(D/'benzene_box_100_summary.json')); HJ=json.load(open(D/'benzene_box_pair_histograms.json'))
fig=plt.figure(figsize=(10.6,7.0)); ax1=fig.add_subplot(2,2,1,projection='3d'); ax2=fig.add_subplot(2,2,2,projection='3d'); ax3=fig.add_subplot(2,2,3); ax4=fig.add_subplot(2,2,4)
cent=np.array([[i*7,j*7,k*7] for i in range(5) for j in range(5) for k in range(4)],float)
ax1.scatter(*cent.T,s=7); q=cent[::4]; ax1.quiver(q[:,0],q[:,1],q[:,2],0,0,2,length=1,normalize=False,linewidth=.6); ax1.set_box_aspect((5,5,4)); ax1.set_xticks([]); ax1.set_yticks([]); ax1.set_zticks([]); lab(ax1,'(a)')
rng=np.random.default_rng(1234); u=rng.normal(size=(len(q),3)); u/=np.linalg.norm(u,axis=1)[:,None]; ax2.scatter(*cent.T,s=7); ax2.quiver(q[:,0],q[:,1],q[:,2],u[:,0],u[:,1],u[:,2],length=2,normalize=False,linewidth=.6); ax2.set_box_aspect((5,5,4)); ax2.set_xticks([]); ax2.set_yticks([]); ax2.set_zticks([]); lab(ax2,'(b)')

for label in ['aligned','isotropic']:
    edges=np.asarray(HJ[label]['bin_edges']); counts=np.asarray(HJ[label]['counts']); centers=.5*(edges[:-1]+edges[1:]); ax3.step(centers,counts,where='mid',label=label)
ax3.set(xlabel=r'$J_{MM}$ (meV)',ylabel='Pair count'); ax3.set_yscale('log'); ax3.legend(frameon=False); lab(ax3,'(c)'); ax3.grid(alpha=.15)
Da=np.asarray(B['aligned']['D_principal_m2_s']); Di=np.asarray(B['isotropic']['ensemble_D_principal_m2_s']); x=np.arange(3); ax4.plot(x,Da/Da[0],'o-',label='aligned'); ax4.plot(x,Di/Di[0],'s-',label='isotropic'); ax4.set_xticks(x,['$D_1$','$D_2$','$D_3$']); ax4.set_yscale('log'); ax4.set_ylabel(r'$D_i/D_1$'); ax4.legend(frameon=False); lab(ax4,'(d)'); ax4.grid(alpha=.15)
fig.tight_layout(); save(fig,'Fig5_benzene_box_eht')
# Fig6: multiple benzene molecules on graphene/slab
T=pd.read_csv(D/'graphene_benzene_tilt_eht.csv'); G=json.load(open(D/'graphene_4benzene_summary.json'))
fig,ax=plt.subplots(2,2,figsize=(9.3,7.0))
# simple top-view hybrid schematic
for i in range(8):
 for j in range(6):
  x=i+.5*(j%2); y=j*.86; ax[0,0].scatter(x,y,s=7,facecolors='none',edgecolors='0.55')
for x,y in [(1.8,1.4),(5.4,1.4),(1.8,3.8),(5.4,3.8)]: ax[0,0].add_patch(RegularPolygon((x,y),6,radius=.42,fill=False,lw=1.5))
ax[0,0].set_aspect('equal'); ax[0,0].set_xticks([]); ax[0,0].set_yticks([]); lab(ax[0,0],'(a)')
ax[0,1].plot(T.tilt_deg,T.Jhole_total_eht_meV,'o-',label='hole'); ax[0,1].plot(T.tilt_deg,T.Jelectron_total_eht_meV,'s-',label='electron'); ax[0,1].set(xlabel='Molecular tilt (deg)',ylabel=r'$J_{G-M}$ (meV)'); ax[0,1].set_yscale('log'); ax[0,1].legend(frameon=False); lab(ax[0,1],'(b)'); ax[0,1].grid(alpha=.15)
vg=np.asarray(G['electron_J_eff_meV']); ax[1,0].bar(np.arange(1,len(vg)+1),vg); ax[1,0].set(xlabel='Adsorbate',ylabel=r'$J_{G-M}$ (meV)'); lab(ax[1,0],'(c)'); ax[1,0].grid(axis='y',alpha=.15)
vm=np.sort(np.asarray(G['pair_J_electron_meV']))[::-1]; ax[1,1].plot(np.arange(1,len(vm)+1),vm,'o'); ax[1,1].set(xlabel='Molecule pair',ylabel=r'$J_{M-M}$ (meV)'); ax[1,1].set_yscale('log'); lab(ax[1,1],'(d)'); ax[1,1].grid(alpha=.15)
fig.tight_layout(); save(fig,'Fig6_graphene_benzene_eht')
print('wrote revised figures to',F)
