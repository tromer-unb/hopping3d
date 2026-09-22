"""Hopping3D direct-EHT command-line parametrizer (graphene implementation)."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from ase.io import read
from .eht import graphene_pi_eht,hopping_pairs

def suggest_graphene_cutoff(atoms):
    ds=[]
    for i in range(len(atoms)):
        for j in range(i+1,len(atoms)):
            d=float(atoms.get_distance(i,j,mic=True))
            if 0.8<d<3.0: ds.append(d)
    vals=np.sort(np.unique(np.round(ds,3)))
    if len(vals)<2: return 1.9
    gaps=np.diff(vals); candidates=np.where(gaps>0.25)[0]
    if len(candidates):
        k=int(candidates[0]); return float(0.5*(vals[k]+vals[k+1]))
    return float(vals[0]+0.35)

def parameterize_graphene(atoms):
    if set(atoms.get_chemical_symbols())!={"C"}: raise ValueError("This parametrizer is validated only for carbon graphene")
    eht=graphene_pi_eht(atoms,cutoff_A=3.0); cutoff=suggest_graphene_cutoff(atoms); pairs=hopping_pairs(atoms,eht,cutoff_A=cutoff); J=np.asarray([x[3] for x in pairs]); Jref=float(np.median(np.abs(J)))
    return {"method":"extended_huckel_hopping3d_graphene","graph":{"cutoff_A":cutoff},"chemistry":{"site_offsets_array_eV":[0.0]*len(atoms)},
      "transport":{"edge_transfer_integrals_eV":{f"{i}-{j}":float(v) for (i,j,_,v) in pairs},"J_ref_eV":Jref,
        "explicit_J_mapping":"rates use |J_ij/J_ref|^2; do not double-count a separate distance exponential","nu0_Hz":None,"reorganization_energy_eV":None},
      "diagnostics":{"mean_abs_J_eV":float(np.mean(np.abs(J))),"std_abs_J_eV":float(np.std(np.abs(J))),
        "lowdin_site_shift_std_eV":float(np.std(eht['site_energy_eV'])),"overlap_min_eigenvalue":float(eht['overlap_eig_min'])},
      "provenance":{"geometry":"input structure","site_energies":"elemental carbon p onsite; relative graphene onsite set to zero",
        "couplings":"EHT Wolfsberg-Helmholtz + Lowdin orthogonalization","confidence":"screening/relative trends; not DFT-equivalent"}}

def main(argv=None):
    ap=argparse.ArgumentParser(description="Hopping3D semiempirical EHT parameter estimator"); ap.add_argument('structure'); ap.add_argument('-o','--output',default='hopping3d_parameters.json'); ns=ap.parse_args(argv); atoms=read(ns.structure); out=parameterize_graphene(atoms); Path(ns.output).write_text(json.dumps(out,indent=2)); print(json.dumps(out,indent=2))

if __name__=='__main__': main()
