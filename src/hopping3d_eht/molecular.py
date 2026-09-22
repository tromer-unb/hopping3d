"""Molecular EHT helpers for Hopping3D.

Frontier-state transfer between degenerate molecular subspaces is represented by
basis-invariant singular values of the projected intermolecular Hamiltonian block.
This is a low-cost semiempirical bridge, not a quantitative quantum-chemistry method.
"""
from __future__ import annotations
import numpy as np
from ase import Atoms
from .eht import ELEMENTS, build_valence_eht

def valence_electron_count(atoms: Atoms) -> int:
    return int(sum(ELEMENTS[s]["valence"] for s in atoms.get_chemical_symbols()))

def frontier_subspaces(atoms: Atoms, n_homo: int = 2, n_lumo: int = 2, cutoff_A: float = 3.2):
    eht = build_valence_eht(atoms, cutoff_A=cutoff_A)
    evals, evecs = np.linalg.eigh(eht["Heff"])
    ne = valence_electron_count(atoms)
    if ne % 2: raise ValueError("Current frontier helper assumes a closed-shell even-electron molecule")
    nocc = ne // 2
    if nocc < n_homo or nocc + n_lumo > len(evals): raise ValueError("Requested frontier subspace is outside the molecular spectrum")
    return {"eht":eht,"energies_eV":evals,"vectors":evecs,"n_electrons":ne,"n_occupied":nocc,
            "homo_indices":list(range(nocc-n_homo,nocc)),"lumo_indices":list(range(nocc,nocc+n_lumo))}

def _effective_subspace_coupling(V):
    V=np.asarray(V,float); s=np.linalg.svd(V,compute_uv=False); m=max(V.shape[0],1)
    return float(np.linalg.norm(V,"fro")/np.sqrt(m)),np.asarray(s,float)

def molecular_dimer_frontier_couplings(atoms,group_a,group_b,intra_cutoff_A=3.2,dimer_cutoff_A=6.0,n_homo=2,n_lumo=2):
    ia=np.asarray(group_a,int); ib=np.asarray(group_b,int); A=atoms[ia]; B=atoms[ib]; A.set_pbc(False); B.set_pbc(False)
    fa=frontier_subspaces(A,n_homo,n_lumo,intra_cutoff_A); fb=frontier_subspaces(B,n_homo,n_lumo,intra_cutoff_A)
    de=build_valence_eht(atoms,cutoff_A=dimer_cutoff_A)
    ba=[k for k,b in enumerate(de["basis"]) if int(b[0]) in set(ia.tolist())]; bb=[k for k,b in enumerate(de["basis"]) if int(b[0]) in set(ib.tolist())]
    Hab=de["Heff"][np.ix_(ba,bb)]
    ah=fa["vectors"][:,fa["homo_indices"]]; bh=fb["vectors"][:,fb["homo_indices"]]
    al=fa["vectors"][:,fa["lumo_indices"]]; bl=fb["vectors"][:,fb["lumo_indices"]]
    Vh=ah.T@Hab@bh; Vl=al.T@Hab@bl; Jh,sh=_effective_subspace_coupling(Vh); Jl,sl=_effective_subspace_coupling(Vl)
    return {"J_hole_eV":Jh,"J_electron_eV":Jl,"hole_singular_values_eV":sh,"electron_singular_values_eV":sl,
            "V_hole_eV":Vh,"V_electron_eV":Vl,"monomer_homo_energies_eV":fa["energies_eV"][fa["homo_indices"]],
            "monomer_lumo_energies_eV":fa["energies_eV"][fa["lumo_indices"]],"dimer_overlap_min_eigenvalue":float(de["overlap_eig_min"])}

def _plane_normal(atoms,carbon_only=True):
    idx=[i for i,s in enumerate(atoms.get_chemical_symbols()) if s=='C'] if carbon_only else list(range(len(atoms)))
    X=np.asarray(atoms.positions[idx],float); X-=X.mean(axis=0); _,v=np.linalg.eigh(X.T@X); n=v[:,0]
    return n/max(np.linalg.norm(n),1e-15)

def conjugated_carbon_pi_eht(atoms: Atoms,groups=None,K=1.75,intra_cutoff_A=1.8,inter_cutoff_A=6.0):
    from .eht import _p_params,p_overlap
    cidx=[i for i,s in enumerate(atoms.get_chemical_symbols()) if s=='C']
    if not cidx: raise ValueError('No carbon atoms found for pi subspace')
    if groups is None: groups=[list(range(len(atoms)))]
    gid={}; normals={}
    for g,idx in enumerate(groups):
        aa=atoms[list(idx)]; normals[g]=_plane_normal(aa,carbon_only=True)
        for i in idx:
            if atoms[int(i)].symbol=='C': gid[int(i)]=g
    onsite=float(_p_params('C')[0]); n=len(cidx); H=np.eye(n)*onsite; S=np.eye(n)
    for a,i in enumerate(cidx):
        for b,j in enumerate(cidx[a+1:],a+1):
            same=gid[int(i)]==gid[int(j)]; cutoff=intra_cutoff_A if same else inter_cutoff_A
            v=np.asarray(atoms.positions[j]-atoms.positions[i],float); r=float(np.linalg.norm(v))
            if r>cutoff or r<1e-12: continue
            s=p_overlap('C','C',v,normals[gid[int(i)]],normals[gid[int(j)]])
            S[a,b]=S[b,a]=s; H[a,b]=H[b,a]=float(K)*s*onsite
    w,U=np.linalg.eigh(S)
    if w.min()<=1e-6: raise ValueError(f'Pi-EHT overlap matrix not positive definite: {w.min():.3e}')
    Sm=(U*(w**-0.5))@U.T; Heff=Sm@H@Sm; evals,evecs=np.linalg.eigh(Heff)
    return {'H':H,'S':S,'Heff':Heff,'energies_eV':evals,'vectors':evecs,'carbon_indices':cidx,'overlap_eig_min':float(w.min()),'K':float(K)}

def calibrate_pi_K_to_gap(atoms,target_gap_eV,lo=1.05,hi=4.0):
    target=float(target_gap_eV)
    def gap(K):
        e=conjugated_carbon_pi_eht(atoms,K=K)['energies_eV']; nocc=len([i for i,s in enumerate(atoms.get_chemical_symbols()) if s=='C'])//2
        return float(e[nocc]-e[nocc-1])
    glo,ghi=gap(lo),gap(hi)
    if not (glo<=target<=ghi): raise ValueError(f'target gap {target} eV outside calibration bracket [{glo},{ghi}]')
    for _ in range(60):
        mid=.5*(lo+hi)
        if gap(mid)<target: lo=mid
        else: hi=mid
    K=.5*(lo+hi); return {'K':float(K),'target_gap_eV':target,'achieved_gap_eV':gap(K)}

def conjugated_pi_dimer_couplings(atoms,group_a,group_b,K,intra_cutoff_A=1.8,inter_cutoff_A=6.0):
    ia=np.asarray(group_a,int); ib=np.asarray(group_b,int); A=atoms[ia]; B=atoms[ib]; A.set_pbc(False); B.set_pbc(False)
    ea=conjugated_carbon_pi_eht(A,K=K,intra_cutoff_A=intra_cutoff_A,inter_cutoff_A=inter_cutoff_A)
    eb=conjugated_carbon_pi_eht(B,K=K,intra_cutoff_A=intra_cutoff_A,inter_cutoff_A=inter_cutoff_A)
    de=conjugated_carbon_pi_eht(atoms,groups=[ia,ib],K=K,intra_cutoff_A=intra_cutoff_A,inter_cutoff_A=inter_cutoff_A)
    na=len(ea['carbon_indices']); nb=len(eb['carbon_indices'])
    if na!=6 or nb!=6: raise ValueError('Current conjugated dimer helper is validated for benzene-like C6 pi systems')
    Ha=ea['vectors'][:,1:3]; Hb=eb['vectors'][:,1:3]; La=ea['vectors'][:,3:5]; Lb=eb['vectors'][:,3:5]
    Hab=de['Heff'][:na,na:na+nb]; Vh=Ha.T@Hab@Hb; Vl=La.T@Hab@Lb
    Jh,sh=_effective_subspace_coupling(Vh); Jl,sl=_effective_subspace_coupling(Vl)
    return {'J_hole_eV':Jh,'J_electron_eV':Jl,'hole_singular_values_eV':sh,'electron_singular_values_eV':sl,
            'V_hole_eV':Vh,'V_electron_eV':Vl,'monomer_homo_eV':ea['energies_eV'][1:3],
            'monomer_lumo_eV':ea['energies_eV'][3:5],'monomer_gap_eV':float(ea['energies_eV'][3]-ea['energies_eV'][2]),
            'dimer_overlap_min_eigenvalue':de['overlap_eig_min'],'K':float(K)}

def benzene_monomer():
    from ase.build import molecule
    b=molecule('C6H6'); b.positions-=b.positions.mean(axis=0); b.pbc=False; return b

def _rotation_axis(axis,angle_rad):
    a=np.asarray(axis,float); a/=max(np.linalg.norm(a),1e-15); x,y,z=a; c=np.cos(angle_rad); s=np.sin(angle_rad); C=1.0-c
    return np.array([[c+x*x*C,x*y*C-z*s,x*z*C+y*s],[y*x*C+z*s,c+y*y*C,y*z*C-x*s],[z*x*C-y*s,z*y*C+x*s,c+z*z*C]])

def make_benzene_dimer(theta_deg,phi_deg,axis_deg,slip_A,slip_angle_deg,target_cc_A):
    A=benzene_monomer(); B=benzene_monomer(); axis=np.array([np.cos(np.deg2rad(axis_deg)),np.sin(np.deg2rad(axis_deg)),0.0])
    R=_rotation_axis(axis,np.deg2rad(theta_deg))@_rotation_axis([0,0,1],np.deg2rad(phi_deg)); B.positions=B.positions@R.T
    B.positions+=np.array([slip_A*np.cos(np.deg2rad(slip_angle_deg)),slip_A*np.sin(np.deg2rad(slip_angle_deg)),0.0])
    ca=[i for i,s in enumerate(A.get_chemical_symbols()) if s=='C']; cb=[i for i,s in enumerate(B.get_chemical_symbols()) if s=='C']
    def mincc(z):
        P=A.positions[ca]; Q=B.positions[cb]+np.array([0.,0.,z]); return float(np.linalg.norm(P[:,None,:]-Q[None,:,:],axis=2).min())
    lo,hi=0.,10.
    for _ in range(70):
        mid=.5*(lo+hi)
        if mincc(mid)<target_cc_A: lo=mid
        else: hi=mid
    z=.5*(lo+hi); B.positions+=np.array([0.,0.,z]); d=A+B; d.pbc=False
    return d,{'theta_deg':float(theta_deg),'phi_deg':float(phi_deg),'axis_deg':float(axis_deg),'slip_A':float(slip_A),
              'slip_angle_deg':float(slip_angle_deg),'target_min_cc_A':float(target_cc_A),'z_shift_A':float(z),
              'center_distance_A':float(np.linalg.norm(B.positions.mean(0)-A.positions.mean(0)))}

def graphene_benzene_pi_interface(atoms,graphene_indices,benzene_indices,K_benzene,K_graphene=1.75,K_cross=None,intra_cutoff_A=1.8,inter_cutoff_A=6.0):
    from .eht import _p_params,p_overlap
    gi=np.asarray(graphene_indices,int); mi=np.asarray(benzene_indices,int); gC=[int(i) for i in gi if atoms[int(i)].symbol=='C']; mC=[int(i) for i in mi if atoms[int(i)].symbol=='C']
    if len(mC)!=6: raise ValueError('Current molecular interface helper expects benzene C6')
    if K_cross is None: K_cross=float(np.sqrt(float(K_graphene)*float(K_benzene)))
    ng=_plane_normal(atoms[gi],carbon_only=True); nm=_plane_normal(atoms[mi],carbon_only=True); cidx=gC+mC; n=len(cidx)
    onsite=float(_p_params('C')[0]); H=np.eye(n)*onsite; S=np.eye(n); gset=set(gC); mset=set(mC)
    for aa,i in enumerate(cidx):
        for bb,j in enumerate(cidx[aa+1:],aa+1):
            gg=i in gset and j in gset; mm=i in mset and j in mset
            if gg: cutoff=intra_cutoff_A; Kij=float(K_graphene); ui=uj=ng
            elif mm: cutoff=intra_cutoff_A; Kij=float(K_benzene); ui=uj=nm
            else: cutoff=inter_cutoff_A; Kij=float(K_cross); ui=ng if i in gset else nm; uj=ng if j in gset else nm
            v=np.asarray(atoms.positions[j]-atoms.positions[i],float); r=float(np.linalg.norm(v))
            if r>cutoff or r<1e-12: continue
            s=p_overlap('C','C',v,ui,uj); S[aa,bb]=S[bb,aa]=s; H[aa,bb]=H[bb,aa]=Kij*s*onsite
    w,U=np.linalg.eigh(S)
    if w.min()<=1e-6: raise ValueError(f'Graphene-benzene pi-EHT overlap not positive definite: {w.min():.3e}')
    Sm=(U*(w**-0.5))@U.T; Heff=Sm@H@Sm; b=atoms[mi]; b.set_pbc(False)
    eb=conjugated_carbon_pi_eht(b,K=K_benzene,intra_cutoff_A=intra_cutoff_A); C_h=eb['vectors'][:,1:3]; C_l=eb['vectors'][:,3:5]
    Hab=Heff[:len(gC),len(gC):]; Vh=Hab@C_h; Vl=Hab@C_l; Jh=np.linalg.norm(Vh,axis=1)/np.sqrt(2.0); Jl=np.linalg.norm(Vl,axis=1)/np.sqrt(2.0)
    return {'H':H,'S':S,'Heff':Heff,'graphene_carbon_indices':gC,'benzene_carbon_indices':mC,'J_hole_rms_eV':Jh,
            'J_electron_rms_eV':Jl,'hole_coupling_matrix_eV':Vh,'electron_coupling_matrix_eV':Vl,
            'graphene_diagonal_shift_eV':np.diag(Heff)[:len(gC)]-np.median(np.diag(Heff)[:len(gC)]),
            'benzene_homo_eV':eb['energies_eV'][1:3],'benzene_lumo_eV':eb['energies_eV'][3:5],
            'overlap_min_eigenvalue':float(w.min()),'K_graphene':float(K_graphene),'K_benzene':float(K_benzene),'K_cross':float(K_cross),
            'coupling_definition':'J_rms=sqrt(sum_alpha |V_site,alpha|^2 / 2) for HOMO/LUMO doublet'}
