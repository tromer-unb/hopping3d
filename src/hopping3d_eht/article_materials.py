"""EHT-only parametrization for the inorganic benchmarks used in the article.

Active targets:
- layered BN: one local p_z transport orbital per atom;
- W2O6: atom-site reduction of O 2p and W 5d frontier manifolds.
"""
from __future__ import annotations
from functools import lru_cache
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from .eht import K_WH,_p_params,p_overlap,local_normals

BOHR_A=0.529177210903
W_EHT={"6s":{"Hii_eV":-8.26,"zeta":2.341},"6p":{"Hii_eV":-5.17,"zeta":2.309},
       "5d":{"Hii_eV":-10.37,"zeta1":4.982,"c1":0.6940,"zeta2":2.068,"c2":0.5631},
       "source":"published Hoffmann-style extended-Huckel parameter table"}
O_2P_HII_EV=-14.80; O_2P_ZETA=2.275

def bn_pi_eht(atoms,cutoff_A=3.8,intralayer_cutoff_A=1.9,K=K_WH):
    sy=np.asarray(atoms.get_chemical_symbols(),dtype=object)
    if set(sy)!={"B","N"}: raise ValueError("bn_pi_eht requires B/N only")
    n=len(atoms); normals=local_normals(atoms); onsite=np.array([_p_params(s)[0] for s in sy]); H=np.diag(onsite); S=np.eye(n)
    for i in range(n):
        for j in range(i+1,n):
            v=np.asarray(atoms.get_distance(i,j,mic=True,vector=True),float); r=float(np.linalg.norm(v))
            if r>cutoff_A or r<1e-10: continue
            ov=p_overlap(str(sy[i]),str(sy[j]),v,normals[i],normals[j]); S[i,j]=S[j,i]=ov; H[i,j]=H[j,i]=K*ov*.5*(onsite[i]+onsite[j])
    w,U=np.linalg.eigh(S)
    if w.min()<=1e-6: raise ValueError(f"BN EHT overlap not positive definite: {w.min():.3e}")
    Sm=(U*(w**-0.5))@U.T; He=Sm@H@Sm; edges=[]; Jin=[]
    for i in range(n):
        for j in range(i+1,n):
            r=float(atoms.get_distance(i,j,mic=True))
            if r<=cutoff_A:
                J=abs(float(He[i,j])); intra=bool(r<=intralayer_cutoff_A); edges.append((i,j,r,str(sy[i]),str(sy[j]),J,intra))
                if intra and set((sy[i],sy[j]))=={"B","N"}: Jin.append(J)
    return _finish_bn(sy,He,w,edges,float(np.median(Jin)))

def _finish_bn(sy,He,w,edges,Jref):
    emap={}; rows=[]; inter=[]
    for i,j,r,si,sj,J,intra in edges:
        wt=(J/Jref)**2; emap[f"{i}-{j}"]=J; rows.append({"i":i,"j":j,"r_A":r,"pair":"-".join(sorted((si,sj))),"J_eV":J,"rate_weight":wt,"intralayer":intra})
        if not intra: inter.append(J)
    diag=np.diag(He); b=float(np.mean(diag[sy=="B"])); n=float(np.mean(diag[sy=="N"])); jm=float(np.median(inter)) if inter else 0.0
    return {"method":"layered B/N local-pz EHT","J_ref_eV":Jref,"edge_transfer_integrals_eV":emap,"edges":rows,
            "overlap_min_eigenvalue":float(np.min(w)),"species_level_mean_eV":{"B":b,"N":n},"raw_B_minus_N_level_eV":b-n,
            "median_interlayer_J_eV":jm,"median_interlayer_to_intralayer_J_ratio":jm/Jref,"median_interlayer_rate_scale":(jm/Jref)**2,
            "transport_note":"use common site energy for the article tensor benchmark; raw B/N orbital offset is diagnostic"}

def _Nr(n,z): return (2.0*z)**(n+0.5)/math.sqrt(math.factorial(2*n))
def _p_z(rho,z,zc,zeta):
    zz=z-zc; r=np.sqrt(rho*rho+zz*zz); return _Nr(2,zeta)*math.sqrt(3/(4*math.pi))*zz*np.exp(-zeta*r)
def _p_x_amp(rho,z,zc,zeta):
    zz=z-zc; r=np.sqrt(rho*rho+zz*zz); return _Nr(2,zeta)*math.sqrt(3/(4*math.pi))*rho*np.exp(-zeta*r)
def _d_z2(rho,z,zc,zeta):
    zz=z-zc; r2=rho*rho+zz*zz; r=np.sqrt(r2); return _Nr(5,zeta)*math.sqrt(5/(16*math.pi))*r2*(3*zz*zz-r2)*np.exp(-zeta*r)
def _d_xz_amp(rho,z,zc,zeta):
    zz=z-zc; r2=rho*rho+zz*zz; r=np.sqrt(r2); return _Nr(5,zeta)*math.sqrt(15/(4*math.pi))*rho*zz*r2*np.exp(-zeta*r)
def _d_combo(fun,rho,z,zc):
    d=W_EHT["5d"]; return d["c1"]*fun(rho,z,zc,d["zeta1"])+d["c2"]*fun(rho,z,zc,d["zeta2"])
def _integrate_cyl(fun,L,N=90):
    x,w=leggauss(int(N)); rho=.5*(x+1.0)*L; wr=.5*L*w; z=L*x; wz=L*w; RR,ZZ=np.meshgrid(rho,z,indexing="ij")
    return float(np.sum(fun(RR,ZZ)*RR*wr[:,None]*wz[None,:]))
@lru_cache(maxsize=1)
def _w_d_norm():
    L=14.0; q=_integrate_cyl(lambda r,z:2*math.pi*_d_combo(_d_z2,r,z,0.0)**2,L,120); return math.sqrt(q)
@lru_cache(maxsize=256)
def sto_channel_overlaps(distance_milli_A):
    R_A=float(distance_milli_A)/1000.0; R=R_A/BOHR_A; L=max(12.0,R+9.0); nd=_w_d_norm()
    pp_s=_integrate_cyl(lambda r,z:2*math.pi*_p_z(r,z,0,O_2P_ZETA)*_p_z(r,z,R,O_2P_ZETA),L)
    pp_p=_integrate_cyl(lambda r,z:math.pi*_p_x_amp(r,z,0,O_2P_ZETA)*_p_x_amp(r,z,R,O_2P_ZETA),L)
    pd_s=_integrate_cyl(lambda r,z:2*math.pi*_p_z(r,z,0,O_2P_ZETA)*_d_combo(_d_z2,r,z,R)/nd,L)
    pd_p=_integrate_cyl(lambda r,z:math.pi*_p_x_amp(r,z,0,O_2P_ZETA)*_d_combo(_d_xz_amp,r,z,R)/nd,L)
    return pp_s,pp_p,pd_s,pd_p

def _channel_J(distance_A,pair,K=K_WH):
    pp_s,pp_p,pd_s,pd_p=sto_channel_overlaps(int(round(float(distance_A)*1000)))
    if pair=="O-O": hs=K*pp_s*O_2P_HII_EV; hp=K*pp_p*O_2P_HII_EV
    elif pair=="O-W": hav=.5*(O_2P_HII_EV+W_EHT["5d"]["Hii_eV"]); hs=K*pd_s*hav; hp=K*pd_p*hav
    else: raise ValueError(f"No active W2O6 atom-site EHT channel for {pair}")
    return float(math.sqrt(hs*hs+2.0*hp*hp)),float(hs),float(hp)

def w2o6_atom_site_eht(atoms,cutoff_A=3.4,periodic=False,K=K_WH):
    sy=np.asarray(atoms.get_chemical_symbols(),dtype=object)
    if set(sy)!={"O","W"}: raise ValueError("w2o6_atom_site_eht requires O/W only")
    rows=[]; by={"O-W":[],"O-O":[],"W-W":[]}
    for i in range(len(atoms)):
        for j in range(i+1,len(atoms)):
            r=float(atoms.get_distance(i,j,mic=bool(periodic)))
            if r>cutoff_A: continue
            pair="-".join(sorted((str(sy[i]),str(sy[j]))))
            if pair=="W-W": by[pair].append((i,j,r,None)); continue
            J,Js,Jp=_channel_J(r,pair,K); by[pair].append((i,j,r,J)); rows.append({"i":i,"j":j,"r_A":r,"pair":pair,"J_block_eV":J,"J_sigma_eV":Js,"J_pi_eV":Jp})
    if not by["O-W"]: raise ValueError("No O-W edges inside cutoff")
    return _finish_w2o6(rows,by,float(np.median([x[3] for x in by["O-W"]])),cutoff_A,periodic)

def _finish_w2o6(rows,by,Jref,cutoff_A,periodic):
    emap={}; oo=[]
    for row in rows:
        key=f"{min(row['i'],row['j'])}-{max(row['i'],row['j'])}"; emap[key]=row["J_block_eV"]; row["rate_weight"]=(row["J_block_eV"]/Jref)**2
        if row["pair"]=="O-O": oo.append(row["rate_weight"])
    raw_delta=W_EHT["5d"]["Hii_eV"]-O_2P_HII_EV; offsets={"O":-raw_delta/2.0,"W":raw_delta/2.0}
    return {"method":"O(2p)-W(5d) atom-site EHT reduction","cutoff_A":float(cutoff_A),"periodic":bool(periodic),"K":float(K_WH),
            "J_ref_WO_eV":Jref,"edge_transfer_integrals_eV":emap,"edges":rows,"edge_counts":{k:len(v) for k,v in by.items()},
            "median_OO_rate_weight":float(np.median(oo)) if oo else None,"OO_rate_weight_range":[float(np.min(oo)),float(np.max(oo))] if oo else None,
            "raw_W5d_minus_O2p_eV":float(raw_delta),"screening_site_offsets_eV":offsets,"W_EHT_parameters":W_EHT,
            "transport_note":"raw orbital-level contrast is a screening scale; do not equate it automatically with a defect/carrier site energy",
            "pair_note":"no W-W edge exists at 3.4 A in the article geometry; homonuclear kinetics there are O-O only"}
