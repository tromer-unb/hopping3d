"""Automatic low-cost interface parametrization for Hopping3D.

Validated target: graphene + one separated conjugated C/H adsorbate,
with benzene receiving the calibrated parameter set. The default runtime
uses direct EHT for interface parametrization.
"""
from __future__ import annotations
from collections import Counter
import numpy as np
from ase import Atoms
from ase.neighborlist import neighbor_list,natural_cutoffs
from .alignment import molecule_graphene_alignment_eht,K_BENZENE_PI
from .molecular import graphene_benzene_pi_interface

BENZENE_REORGANIZATION_EV={"hole":0.057962630462082734,"electron":0.2046250207569526,
                           "source":"offline GFN1-xTB four-point calibration"}

def _formula(symbols):
    c=Counter(symbols); order=[]
    for s in ("C","H"):
        if s in c: order.append(s)
    order+=sorted(k for k in c if k not in order)
    return "".join(s+(str(c[s]) if c[s]!=1 else "") for s in order)

def bonded_components(atoms:Atoms,mult=1.15):
    cut=natural_cutoffs(atoms,mult=float(mult)); ii,jj=neighbor_list("ij",atoms,cut); parent=list(range(len(atoms)))
    def find(x):
        while parent[x]!=x: parent[x]=parent[parent[x]]; x=parent[x]
        return x
    def union(a,b):
        a,b=find(int(a)),find(int(b))
        if a!=b: parent[b]=a
    for i,j in zip(ii,jj): union(i,j)
    groups={}
    for i in range(len(atoms)): groups.setdefault(find(i),[]).append(i)
    return sorted((np.asarray(v,int) for v in groups.values()),key=len,reverse=True)

def _plane_normal(points):
    X=np.asarray(points,float); X=X-X.mean(axis=0); w,v=np.linalg.eigh(X.T@X); n=v[:,0]
    return n/max(np.linalg.norm(n),1e-15),float(np.sqrt(max(w[0],0.0)/max(len(X),1)))

def _unwrapped_component_positions(atoms,idx):
    idx=np.asarray(idx,int); allowed=set(idx.tolist()); cut=natural_cutoffs(atoms,mult=1.15); ii,jj,DD=neighbor_list("ijD",atoms,cut); adj={int(i):[] for i in idx}
    for i,j,d in zip(ii,jj,DD):
        if int(i) in allowed and int(j) in allowed: adj[int(i)].append((int(j),np.asarray(d,float)))
    root=int(idx[0]); pos={root:np.asarray(atoms.positions[root],float)}; stack=[root]
    while stack:
        i=stack.pop()
        for j,d in adj[i]:
            if j not in pos: pos[j]=pos[i]+d; stack.append(j)
    if len(pos)!=len(idx): raise ValueError("Could not unwrap molecular component from covalent graph")
    return np.asarray([pos[int(i)] for i in idx],float)

def _component_record(atoms,idx):
    sy=[atoms[int(i)].symbol for i in idx]
    return {"n_atoms":int(len(idx)),"formula":_formula(sy),"symbols":dict(Counter(sy)),"indices":[int(i) for i in idx]}

def detect_graphene_molecule(atoms:Atoms):
    comps=bonded_components(atoms); records=[_component_record(atoms,c) for c in comps]; host_pos=None
    for k,c in enumerate(comps):
        sy={atoms[int(i)].symbol for i in c}
        if sy=={"C"} and len(c)>=12: host_pos=k; break
    if host_pos is None: raise ValueError("No graphene-like all-carbon component with >=12 atoms detected")
    host_idx=comps[host_pos]; ads=[(k,c) for k,c in enumerate(comps) if k!=host_pos and len(c)>=2]
    if not ads: raise ValueError("No separated molecular adsorbate component detected")
    if len(ads)>1: raise ValueError(f"Detected {len(ads)} molecular adsorbates; use `hopping3d-eht coverage` for multi-adsorbate systems")
    mol_pos,mol_idx=ads[0]; mol_sy={atoms[int(i)].symbol for i in mol_idx}
    if not mol_sy.issubset({"C","H"}): raise ValueError("First interface mode supports a separated conjugated C/H adsorbate")
    host=atoms[host_idx]; host.set_pbc(atoms.pbc); cut=natural_cutoffs(host,mult=1.15); ii,_=neighbor_list("ij",host,cut); deg=np.bincount(ii,minlength=len(host)); n,rms=_plane_normal(host.positions)
    if not (np.median(deg)>=2.5 and rms<0.25): raise ValueError("Largest carbon component does not pass the current graphene-like geometry checks")
    return host_idx,mol_idx,{"components":records,"host_component":int(host_pos),"molecule_component":int(mol_pos),
        "host_type":"graphene-like carbon pi sheet","molecule_formula":_formula([atoms[int(i)].symbol for i in mol_idx]),
        "host_median_coordination":float(np.median(deg)),"host_plane_rms_A":rms}

def _interface_geometry(atoms,host_idx,mol_idx):
    hC=[int(i) for i in host_idx if atoms[int(i)].symbol=='C']; mC=[int(i) for i in mol_idx if atoms[int(i)].symbol=='C']; hn,hrms=_plane_normal(atoms.positions[hC])
    mpos=_unwrapped_component_positions(atoms,mol_idx); mC_local=[k for k,i in enumerate(mol_idx) if atoms[int(i)].symbol=='C']; mn,mrms=_plane_normal(mpos[mC_local])
    tilt=float(np.degrees(np.arccos(np.clip(abs(np.dot(hn,mn)),0.0,1.0)))); hc=np.mean(atoms.positions[hC],axis=0); mc=np.mean(mpos[mC_local],axis=0)
    height=float(abs(np.dot(mc-hc,hn))); dmin=min(float(atoms.get_distance(i,j,mic=True)) for i in hC for j in mC)
    return {"adsorption_height_A":height,"tilt_deg":tilt,"min_host_molecule_CC_A":dmin,"host_plane_rms_A":hrms,"molecule_plane_rms_A":mrms,
            "host_normal":[float(x) for x in hn],"molecule_normal":[float(x) for x in mn]}

def _local_interface_cluster(atoms,host_idx,mol_idx,radius_A=8.0):
    host=atoms[host_idx].copy(); mol=atoms[mol_idx].copy(); mol.positions=_unwrapped_component_positions(atoms,mol_idx); hn,_=_plane_normal(host.positions)
    cell=np.asarray(atoms.cell.array,float); reps=[1,1,1]; shift=np.zeros(3)
    for a in range(3):
        L=np.linalg.norm(cell[a]); inplane=L>1e-8 and abs(np.dot(cell[a]/L,hn))<0.5
        if inplane and bool(atoms.pbc[a]): reps[a]=3; shift+=cell[a]
    hrep=host.repeat(tuple(reps)); mol.positions+=shift; mc=np.mean(mol.positions[[i for i,s in enumerate(mol.get_chemical_symbols()) if s=='C']],axis=0); hp=np.asarray(hrep.positions,float)
    patch=hrep[np.where(np.linalg.norm(hp-mc,axis=1)<=float(radius_A))[0]]; patch.set_pbc(False); mol.set_pbc(False); combo=patch+mol; combo.set_pbc(False)
    return combo,np.arange(len(patch),dtype=int),np.arange(len(patch),len(combo),dtype=int),{"n_host_patch_atoms":int(len(patch)),"repeat":reps,"radius_A":float(radius_A)}

def _parameter(value,unit,source,confidence,recommended_use,caveat=None):
    out={"value":value,"unit":unit,"source":source,"confidence":confidence,"recommended_use":recommended_use}
    if caveat: out["caveat"]=caveat
    return out

def parameterize_interface(atoms:Atoms,patch_radius_A=8.0):
    host_idx,mol_idx,det=detect_graphene_molecule(atoms); host=atoms[host_idx].copy(); mol=atoms[mol_idx].copy(); mol.set_pbc(False); geom=_interface_geometry(atoms,host_idx,mol_idx); align=molecule_graphene_alignment_eht(host,mol); formula=det["molecule_formula"]
    params={"graphene_work_function_proxy":_parameter(align["graphene_reference"]["work_function_proxy_eV"],"eV","portable EHT with offline PAH calibration","medium","screening/reference alignment"),
            "delta_mu_molecule_minus_graphene":_parameter(align["delta_mu_molecule_minus_graphene_eV"],"eV","portable EHT with offline PAH calibration","medium-low","qualitative charge-transfer screening"),
            "hole_localization_offset":_parameter(align["hole_localization_offset_eV"],"eV","isolated-fragment EHT/redox calibration","low","regime screening only"),
            "electron_localization_offset":_parameter(align["electron_localization_offset_eV"],"eV","isolated-fragment EHT/redox calibration","low","regime screening only")}
    result={"schema":"hopping3d-interface-0.1","detection":det,"geometry":geom,"parameters":params,
            "warnings":["Level alignment omits image-charge polarization, dielectric screening and interface dipoles."],"external_executables_required":[]}
    if formula=="C6H6":
        result["detection"]["molecule_type"]="benzene"; result["parameters"]["lambda_hole"]=_parameter(BENZENE_REORGANIZATION_EV["hole"],"eV",BENZENE_REORGANIZATION_EV["source"],"medium","Marcus hole-transfer kinetics"); result["parameters"]["lambda_electron"]=_parameter(BENZENE_REORGANIZATION_EV["electron"],"eV",BENZENE_REORGANIZATION_EV["source"],"medium","Marcus electron-transfer kinetics")
        combo,gi,mi,patch=_local_interface_cluster(atoms,host_idx,mol_idx,patch_radius_A); eint=graphene_benzene_pi_interface(combo,gi,mi,K_BENZENE_PI); Vh=np.asarray(eint["hole_coupling_matrix_eV"],float); Ve=np.asarray(eint["electron_coupling_matrix_eV"],float); Sh=float(np.sum(Vh*Vh)); Se=float(np.sum(Ve*Ve))
        result["interface_patch"]=patch; result["parameters"]["hole_interface_strength"]=_parameter(Sh,"eV^2","direct coherent pi-EHT interface block","medium","single-geometry interface kinetics"); result["parameters"]["electron_interface_strength"]=_parameter(Se,"eV^2","direct coherent pi-EHT interface block","medium","single-geometry interface kinetics"); result["parameters"]["hole_J_eff"]=_parameter(float(np.sqrt(Sh/2.0)),"eV","sqrt(sum |V|^2 / HOMO-doublet degeneracy)","medium","coarse molecular-site coupling"); result["parameters"]["electron_J_eff"]=_parameter(float(np.sqrt(Se/2.0)),"eV","sqrt(sum |V|^2 / LUMO-doublet degeneracy)","medium","coarse molecular-site coupling")
        result["diagnostics"]={"interface_overlap_min_eigenvalue":float(eint["overlap_min_eigenvalue"]),"K_graphene":float(eint["K_graphene"]),"K_benzene":float(eint["K_benzene"]),"K_cross":float(eint["K_cross"])}
        result["transport_templates"]={"hole":{"rate_model":"marcus_gerischer","molecular_degeneracy":2,"lambda_eV":BENZENE_REORGANIZATION_EV["hole"],"interface_strength_eV2":Sh,"energy_offset_eV":None,"screening_only_offset_eV":float(align["hole_localization_offset_eV"])},
            "electron":{"rate_model":"marcus_gerischer","molecular_degeneracy":2,"lambda_eV":BENZENE_REORGANIZATION_EV["electron"],"interface_strength_eV2":Se,"energy_offset_eV":None,"screening_only_offset_eV":float(align["electron_localization_offset_eV"])}}
    else:
        result["detection"]["molecule_type"]="conjugated C/H molecule (unclassified)"; result["warnings"].append("No shipped reorganization energy or interface-coupling calibration for this molecule.")
    nperiodic=int(np.sum(np.asarray(atoms.pbc,dtype=bool)))
    if nperiodic>=2: recommendation={"model":"electronic_reservoir","confidence":"high","reason":"extended graphene-like host is better treated as a delocalized pi reservoir than atom-by-atom hopping","kinetics":"Marcus-Gerischer for molecule<->host; direct EHT retains coherent multi-site interference"}
    else: recommendation={"model":"finite_pi_host_or_localized_benchmark","confidence":"medium","reason":"host is not clearly an extended periodic sheet; choose reservoir vs localized CTMC from the physical problem"}
    result["solver_recommendation"]=recommendation; result["provenance"]={"runtime":"Python/NumPy/SciPy/ASE + Hopping3D EHT code","xTB_required":False,"electronic_parameterization":"direct EHT","scope":"screening/semiempirical; not DFT-equivalent"}
    return result

def render_interface_markdown(result):
    d=result["detection"]; g=result["geometry"]; lines=["# Hopping3D EHT interface report","",f"Detected host: **{d['host_type']}**",f"Detected adsorbate: **{d.get('molecule_type',d['molecule_formula'])} ({d['molecule_formula']})**","","## Geometry","",f"- adsorption height: {g['adsorption_height_A']:.3f} Å",f"- molecular tilt: {g['tilt_deg']:.2f}°",f"- minimum host–molecule C–C distance: {g['min_host_molecule_CC_A']:.3f} Å","","## Estimated parameters","","| parameter | value | source | confidence | recommended use |","|---|---:|---|---|---|"]
    for name,p in result["parameters"].items():
        v=p["value"]; txt=f"{v:.6g} {p['unit']}" if isinstance(v,(int,float)) else f"{v} {p['unit']}"; lines.append(f"| {name} | {txt} | {p['source']} | {p['confidence']} | {p['recommended_use']} |")
    r=result["solver_recommendation"]; lines += ["","## Recommended transport representation","",f"**{r['model']}** ({r['confidence']} confidence): {r['reason']}","",f"Kinetics: {r.get('kinetics','not assigned')}","","## Warnings",""]+[f"- {w}" for w in result.get("warnings",[])]
    return "\n".join(lines)+"\n"
