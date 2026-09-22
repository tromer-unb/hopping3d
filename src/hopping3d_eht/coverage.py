"""Exact coarse-grained coverage dynamics for molecule/graphene interfaces."""
from __future__ import annotations
import math
import numpy as np
from .reservoir import HBAR_EV_S
KB_EV=8.617333262145e-5

def molecular_marcus_rate(J_eff_eV,deltaG_eV,lambda_eV,T_K=300.0):
    J=abs(float(J_eff_eV)); d=float(deltaG_eV); lam=float(lambda_eV); T=float(T_K)
    if lam<=0 or T<=0: raise ValueError("lambda_eV and T_K must be positive")
    pref=(2.0*math.pi/HBAR_EV_S)*J*J/math.sqrt(4.0*math.pi*lam*KB_EV*T)
    return float(pref*math.exp(-((d+lam)**2)/(4.0*lam*KB_EV*T)))

def _stationary_distribution(generator):
    G=np.asarray(generator,float); n=len(G); A=G.T.copy(); b=np.zeros(n); A[-1,:]=1.0; b[-1]=1.0
    p=np.linalg.solve(A,b); p=np.maximum(p,0.0); return p/p.sum()

def molecular_excursion_metrics(k_in_s,k_out_s,k_mm_s,centers_A):
    kin=np.asarray(k_in_s,float); kout=np.asarray(k_out_s,float); K=np.asarray(k_mm_s,float); pos=np.asarray(centers_A,float); n=len(kin)
    if kout.shape!=(n,) or K.shape!=(n,n) or pos.shape!=(n,3): raise ValueError("coverage arrays have inconsistent shapes")
    if np.any(kin<0) or np.any(kout<0) or np.any(K<0): raise ValueError("rates must be non-negative")
    total_in=float(kin.sum())
    if total_in<=0: raise ValueError("total reservoir capture rate must be positive")
    pentry=kin/total_in; rmm=K.sum(axis=1); Q=K.copy(); np.fill_diagonal(Q,-(kout+rmm)); F=np.linalg.inv(-Q); ones=np.ones(n)
    mean_exc=float(pentry@F@ones); B=F@np.diag(kout); expected_mm=float(pentry@F@rmm)
    p_any=float(pentry@np.divide(rmm,rmm+kout,out=np.zeros(n),where=(rmm+kout)>0)); p_diff=float(sum(pentry[i]*(B[i].sum()-B[i,i]) for i in range(n)))
    msd=0.0
    for i in range(n):
        dr=pos-pos[i]; r2=np.sum(dr*dr,axis=1); msd+=float(pentry[i]*np.dot(B[i],r2))
    G=np.zeros((n+1,n+1),float); G[0,1:]=kin; G[0,0]=-total_in
    for i in range(n):
        G[i+1,0]=kout[i]
        for j in range(n):
            if i!=j: G[i+1,j+1]=K[i,j]
        G[i+1,i+1]=-G[i+1].sum()
    stat=_stationary_distribution(G); empty_wait=1.0/total_in; renewal_fraction=mean_exc/(empty_wait+mean_exc)
    return {"total_capture_rate_s-1":total_in,"mean_wait_to_capture_s":empty_wait,"mean_molecular_excursion_s":mean_exc,
            "probability_any_inter_molecular_hop":p_any,"expected_inter_molecular_hops_per_excursion":expected_mm,
            "probability_return_via_different_molecule":p_diff,"rms_lateral_migration_A":float(math.sqrt(max(msd,0.0))),
            "molecular_stationary_probability":float(stat[1:].sum()),"renewal_molecular_time_fraction":float(renewal_fraction),
            "stationary_distribution":[float(x) for x in stat],"entry_probabilities":[float(x) for x in pentry],"exit_probability_matrix":B.tolist()}
