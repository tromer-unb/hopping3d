"""Low-cost extended-Huckel parametrization.

The engine uses the Wolfsberg-Helmholtz form H_ij = K S_ij (H_ii+H_jj)/2.
For graphene validation, the transport subspace is one locally oriented p_z orbital
per carbon. Slater-Koster angular factors and exponentially decaying overlap proxies
are used. This is an EHT parametrization layer, not a substitute for converged DFT.
"""
from __future__ import annotations
import numpy as np
from ase.data import atomic_numbers, covalent_radii

K_WH = 1.75

ELEMENTS = {
    "H":  {"valence":1, "orbitals":{"1s":(-13.60,1.30)}},
    "Li": {"valence":1, "orbitals":{"2s":(-5.39,0.65), "2p":(-3.50,0.65)}},
    "B":  {"valence":3, "orbitals":{"2s":(-15.20,1.30), "2p":(-8.30,1.30)}},
    "C":  {"valence":4, "orbitals":{"2s":(-21.40,1.625), "2p":(-11.40,1.625)}},
    "N":  {"valence":5, "orbitals":{"2s":(-26.00,1.95), "2p":(-13.40,1.95)}},
    "O":  {"valence":6, "orbitals":{"2s":(-32.30,2.275), "2p":(-14.80,2.275)}},
    "F":  {"valence":7, "orbitals":{"2s":(-40.00,2.425), "2p":(-18.10,2.425)}},
    "Na": {"valence":1, "orbitals":{"3s":(-5.14,0.73), "3p":(-3.00,0.73)}},
    "Si": {"valence":4, "orbitals":{"3s":(-17.30,1.38), "3p":(-9.20,1.38)}},
    "P":  {"valence":5, "orbitals":{"3s":(-18.60,1.60), "3p":(-10.10,1.60)}},
    "S":  {"valence":6, "orbitals":{"3s":(-20.00,1.80), "3p":(-11.00,1.80)}},
    "Cl": {"valence":7, "orbitals":{"3s":(-26.30,2.00), "3p":(-14.20,2.00)}},
}

OVERLAP_REF = {"ss_sigma":0.20, "sp_sigma":0.16, "pp_sigma":0.24, "pp_pi":0.136}
PAIR_RREF = {("C","C"):1.42}

def _p_params(symbol):
    d=ELEMENTS[symbol]["orbitals"]
    key=next(k for k in d if k.endswith("p"))
    return d[key]

def local_normals(atoms, n_neighbors=3):
    n=len(atoms); normals=np.zeros((n,3))
    for i in range(n):
        items=[]
        for j in range(n):
            if i==j: continue
            v=np.asarray(atoms.get_distance(i,j,mic=True,vector=True),float)
            items.append((np.linalg.norm(v),v))
        items.sort(key=lambda x:x[0])
        V=np.array([v for _,v in items[:n_neighbors]])
        cov=V.T@V
        _,vec=np.linalg.eigh(cov)
        normal=vec[:,0]
        if normal[2] < 0: normal=-normal
        normals[i]=normal/np.linalg.norm(normal)
    return normals

def _pair_rref(si,sj):
    key=tuple(sorted((si,sj)))
    if key in PAIR_RREF: return PAIR_RREF[key]
    return float(covalent_radii[atomic_numbers[si]] + covalent_radii[atomic_numbers[sj]])

def p_overlap(si,sj,rvec,ui,uj):
    r=float(np.linalg.norm(rvec)); rh=np.asarray(rvec)/r
    _,zi=_p_params(si); _,zj=_p_params(sj)
    beta=0.75*(zi+zj)
    dr=r-_pair_rref(si,sj)
    f=np.exp(-beta*dr)
    spp_sigma=OVERLAP_REF["pp_sigma"]*f
    spp_pi=OVERLAP_REF["pp_pi"]*f
    ci=float(np.dot(ui,rh)); cj=float(np.dot(uj,rh)); uu=float(np.dot(ui,uj))
    return ci*cj*spp_sigma + (uu-ci*cj)*spp_pi

def graphene_pi_eht(atoms, cutoff_A=3.0, K=K_WH):
    syms=atoms.get_chemical_symbols()
    if set(syms)!={"C"}: raise ValueError("graphene_pi_eht currently validates carbon-only systems")
    n=len(atoms); normals=local_normals(atoms)
    onsite=_p_params("C")[0]
    H=np.eye(n)*onsite; S=np.eye(n); raw=[]
    for i in range(n):
        for j in range(i+1,n):
            v=np.asarray(atoms.get_distance(i,j,mic=True,vector=True),float); r=np.linalg.norm(v)
            if r>cutoff_A: continue
            s=p_overlap("C","C",v,normals[i],normals[j])
            hij=K*s*onsite
            S[i,j]=S[j,i]=s; H[i,j]=H[j,i]=hij
            raw.append((i,j,r,s,hij))
    w,U=np.linalg.eigh(S)
    if w.min() <= 1e-6:
        raise ValueError(f"Overlap matrix not positive definite; min eigenvalue={w.min():.3e}")
    Smhalf=(U*(w**-0.5))@U.T
    Heff=Smhalf@H@Smhalf
    eps=np.diag(Heff).copy(); eps-=np.median(eps)
    return {"H":H,"S":S,"Heff":Heff,"site_energy_eV":eps,
            "normals":normals,"raw_pairs":raw,"overlap_eig_min":float(w.min())}

def hopping_pairs(atoms, eht_result, cutoff_A=1.85):
    Heff=eht_result["Heff"]; rows=[]
    for i in range(len(atoms)):
        for j in range(i+1,len(atoms)):
            r=float(atoms.get_distance(i,j,mic=True))
            if r<=cutoff_A:
                rows.append((i,j,r,float(Heff[i,j])))
    return rows

def fit_coupling_decay(pair_rows):
    arr=np.array([(r,abs(j)) for _,_,r,j in pair_rows if abs(j)>1e-10])
    x=arr[:,0]; y=np.log(arr[:,1]); slope,intercept=np.polyfit(x,y,1)
    beta=max(0.0,float(-slope))
    return {"beta_Ainv":beta,"xi_rate_A":float(1/beta) if beta>0 else np.inf,
            "J_at_zero_eV":float(np.exp(intercept))}

_P_AXES={"px":np.array([1.0,0.0,0.0]),"py":np.array([0.0,1.0,0.0]),"pz":np.array([0.0,0.0,1.0])}

def _expand_orbitals(symbol):
    if symbol not in ELEMENTS:
        raise KeyError(f"No EHT parameters for {symbol}; add an element parameter record before use")
    d=ELEMENTS[symbol]["orbitals"]; out=[]
    for key,(Hii,zeta) in d.items():
        if key.endswith('s'):
            out.append((key,float(Hii),float(zeta),None))
        elif key.endswith('p'):
            shell=key[:-1]
            for q,axis in _P_AXES.items(): out.append((shell+q,float(Hii),float(zeta),axis))
    return out

def _overlap_radial(kind, si, sj, zeta_i, zeta_j, r):
    beta=0.75*(float(zeta_i)+float(zeta_j)); f=np.exp(-beta*(float(r)-_pair_rref(si,sj)))
    return OVERLAP_REF[kind]*f

def build_valence_eht(atoms, cutoff_A=3.5, K=K_WH):
    syms=atoms.get_chemical_symbols(); basis=[]
    for ia,sym in enumerate(syms):
        for label,Hii,zeta,axis in _expand_orbitals(sym):
            basis.append((ia,sym,label,Hii,zeta,axis))
    nb=len(basis); H=np.zeros((nb,nb)); S=np.eye(nb)
    for a,(_,_,_,Hii,_,_) in enumerate(basis): H[a,a]=Hii
    atom_to_basis={i:[] for i in range(len(atoms))}
    for a,b in enumerate(basis): atom_to_basis[b[0]].append(a)
    for i in range(len(atoms)):
        for j in range(i+1,len(atoms)):
            v=np.asarray(atoms.get_distance(i,j,mic=True,vector=True),float); r=float(np.linalg.norm(v))
            if r>cutoff_A or r<1e-10: continue
            rh=v/r
            for a in atom_to_basis[i]:
                _,si,li,Hi,zi,ai=basis[a]
                for b in atom_to_basis[j]:
                    _,sj,lj,Hj,zj,aj=basis[b]
                    ispi=ai is not None; ispj=aj is not None
                    if not ispi and not ispj:
                        ov=_overlap_radial('ss_sigma',si,sj,zi,zj,r)
                    elif not ispi and ispj:
                        ov=-float(np.dot(aj,rh))*_overlap_radial('sp_sigma',si,sj,zi,zj,r)
                    elif ispi and not ispj:
                        ov=float(np.dot(ai,rh))*_overlap_radial('sp_sigma',si,sj,zi,zj,r)
                    else:
                        ci=float(np.dot(ai,rh)); cj=float(np.dot(aj,rh)); uu=float(np.dot(ai,aj))
                        sig=_overlap_radial('pp_sigma',si,sj,zi,zj,r); pi=_overlap_radial('pp_pi',si,sj,zi,zj,r)
                        ov=ci*cj*sig+(uu-ci*cj)*pi
                    S[a,b]=S[b,a]=ov
                    H[a,b]=H[b,a]=K*ov*0.5*(Hi+Hj)
    w,U=np.linalg.eigh(S)
    if w.min()<=1e-6: raise ValueError(f"Generic EHT overlap matrix not positive definite: {w.min():.3e}")
    Smhalf=(U*(w**-0.5))@U.T; Heff=Smhalf@H@Smhalf
    return {"basis":basis,"H":H,"S":S,"Heff":Heff,"overlap_eig_min":float(w.min())}

def elemental_p_site_energies(atoms):
    vals=[]
    for s in atoms.get_chemical_symbols(): vals.append(_p_params(s)[0])
    vals=np.asarray(vals,float); return vals-np.median(vals)

def hybrid_graphene_li_eht(atoms, cutoff_A=3.5, K=K_WH):
    syms=np.asarray(atoms.get_chemical_symbols(),dtype=object)
    if not set(syms).issubset({'C','Li'}) or 'C' not in set(syms):
        raise ValueError('hybrid_graphene_li_eht requires C/Li with at least one carbon')
    n=len(atoms); H=np.zeros((n,n)); S=np.eye(n)
    normals=np.zeros((n,3),float); cidx=np.where(syms=='C')[0]
    for i in cidx:
        items=[]
        for j in cidx:
            if i==j: continue
            v=np.asarray(atoms.get_distance(int(i),int(j),mic=True,vector=True),float)
            items.append((float(np.linalg.norm(v)),v))
        items.sort(key=lambda x:x[0]); V=np.array([v for _,v in items[:3]])
        _,vec=np.linalg.eigh(V.T@V); u=vec[:,0]
        if u[2]<0: u=-u
        normals[i]=u/max(np.linalg.norm(u),1e-15)
    Hii=np.array([_p_params('C')[0] if s=='C' else ELEMENTS['Li']['orbitals']['2s'][0] for s in syms],float)
    np.fill_diagonal(H,Hii)
    raw=[]
    for i in range(n):
        for j in range(i+1,n):
            v=np.asarray(atoms.get_distance(i,j,mic=True,vector=True),float); r=float(np.linalg.norm(v))
            if r>cutoff_A or r<1e-10: continue
            si,sj=syms[i],syms[j]; rh=v/r
            if si=='C' and sj=='C':
                s=p_overlap('C','C',v,normals[i],normals[j])
            elif si=='Li' and sj=='Li':
                zi=ELEMENTS['Li']['orbitals']['2s'][1]; zj=zi
                s=_overlap_radial('ss_sigma','Li','Li',zi,zj,r)
            else:
                if si=='C':
                    c=i; sign=1.0; zc=_p_params('C')[1]; zl=ELEMENTS['Li']['orbitals']['2s'][1]
                else:
                    c=j; sign=-1.0; zc=_p_params('C')[1]; zl=ELEMENTS['Li']['orbitals']['2s'][1]
                s=sign*float(np.dot(normals[c],rh))*_overlap_radial('sp_sigma','C','Li',zc,zl,r)
            hij=K*s*0.5*(Hii[i]+Hii[j])
            S[i,j]=S[j,i]=s; H[i,j]=H[j,i]=hij
            raw.append((i,j,r,str(si),str(sj),float(s),float(hij)))
    w,U=np.linalg.eigh(S)
    if w.min()<=1e-6: raise ValueError(f'Graphene-Li EHT overlap not positive definite: {w.min():.3e}')
    Smhalf=(U*(w**-0.5))@U.T; Heff=Smhalf@H@Smhalf
    ref=float(np.median(np.diag(Heff)[cidx])); eps=np.diag(Heff).copy()-ref
    return {'H':H,'S':S,'Heff':Heff,'site_energy_eV':eps,'normals':normals,
            'raw_pairs':raw,'overlap_eig_min':float(w.min()),'carbon_reference_eV':ref}

def graphene_li_hopping_pairs(atoms, eht_result, cc_cutoff_A=1.90, cli_cutoff_A=3.20, lili_cutoff_A=3.50):
    syms=atoms.get_chemical_symbols(); H=eht_result['Heff']; rows=[]
    for i in range(len(atoms)):
        for j in range(i+1,len(atoms)):
            r=float(atoms.get_distance(i,j,mic=True)); pair='-'.join(sorted((syms[i],syms[j])))
            cutoff={'C-C':cc_cutoff_A,'C-Li':cli_cutoff_A,'Li-Li':lili_cutoff_A}.get(pair,0.0)
            if r<=cutoff: rows.append((i,j,r,pair,float(H[i,j])))
    return rows
