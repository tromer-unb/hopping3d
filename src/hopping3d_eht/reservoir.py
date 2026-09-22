"""Low-cost molecule/electrode Marcus reservoir models for Hopping3D."""
from __future__ import annotations
import math
import numpy as np
KB_EV=8.617333262145e-5
from .eht import graphene_pi_eht,_p_params
HBAR_EV_S=6.582119569e-16

def fermi(E_eV,T_K=300.0):
    x=np.asarray(E_eV,float)/(KB_EV*float(T_K)); return 1.0/(1.0+np.exp(np.clip(x,-700,700)))

def marcus_prefactor_per_J2(deltaG_eV,lambda_eV,T_K=300.0):
    lam=float(lambda_eV); T=float(T_K); d=np.asarray(deltaG_eV,float); pref=(2.0*math.pi/HBAR_EV_S)/math.sqrt(4.0*math.pi*lam*KB_EV*T)
    return pref*np.exp(-((d+lam)**2)/(4.0*lam*KB_EV*T))

def graphene_spectral_strengths(graphene_atoms,V_site_micro_eV,mode='coherent',local_strengths_eV2=None):
    host=graphene_pi_eht(graphene_atoms,cutoff_A=1.85); eps,U=np.linalg.eigh(np.asarray(host['Heff'],float)); mu=float(_p_params('C')[0]); rel=eps-mu; V=np.asarray(V_site_micro_eV,float)
    if V.shape[0]!=len(graphene_atoms): raise ValueError('V_site_micro_eV first dimension must match graphene sites')
    if mode=='coherent': W=U.T@V; S=np.sum(W*W,axis=1)
    elif mode=='incoherent':
        if local_strengths_eV2 is None: local_strengths_eV2=np.sum(V*V,axis=1)
        s=np.asarray(local_strengths_eV2,float)
        if len(s)!=len(graphene_atoms): raise ValueError('local strength length mismatch')
        S=(U*U).T@s
    else: raise ValueError('mode must be coherent or incoherent')
    return {'energies_relative_mu_eV':rel,'spectral_strength_eV2':np.asarray(S,float),'eigenvectors':U,'mu_EHT_eV':mu,'host_eigenvalues_eV':eps}

def electron_reservoir_rates(spectral,molecular_level_eV,lambda_eV,degeneracy=2.0,T_K=300.0):
    en=np.asarray(spectral['energies_relative_mu_eV'],float); S=np.asarray(spectral['spectral_strength_eV2'],float); f=fermi(en,T_K); de=float(molecular_level_eV)-en
    kin_terms=S*f*marcus_prefactor_per_J2(de,lambda_eV,T_K); kout_terms=(S/float(degeneracy))*(1.0-f)*marcus_prefactor_per_J2(-de,lambda_eV,T_K)
    kin=float(np.sum(kin_terms)); kout=float(np.sum(kout_terms)); ratio=float(kin/kout) if kout>0 else float('inf'); expected=float(degeneracy)*math.exp(-float(molecular_level_eV)/(KB_EV*float(T_K))); peq=float(kin/(kin+kout)) if kin+kout>0 else float('nan')
    return {'k_reservoir_to_molecule_s-1':kin,'k_molecule_to_reservoir_s-1':kout,'rate_ratio':ratio,'expected_detailed_balance_ratio':expected,
            'relative_detailed_balance_error':float(abs(ratio-expected)/expected) if expected>0 and np.isfinite(ratio) else float('nan'),
            'equilibrium_molecular_probability':peq,'mean_empty_wait_s':float(1.0/kin) if kin>0 else float('inf'),'mean_occupied_dwell_s':float(1.0/kout) if kout>0 else float('inf')}
