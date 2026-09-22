"""Low-cost EHT-only level-alignment helpers for conjugated C/H systems.

The runtime path does not require xTB. The affine corrections below were fitted
offline to GFN1-xTB VIP/VEA values for a small PAH calibration set and are
therefore explicitly a rough carbon-pi calibration, not a universal model.
"""
from __future__ import annotations
import numpy as np
from .eht import graphene_pi_eht,_p_params
from .molecular import conjugated_carbon_pi_eht

K_BENZENE_PI=2.417879719037213
PAH_CALIBRATION={
    "ip_from_minus_homo":(1.98000850,-16.46516978),
    "ea_from_minus_lumo":(1.57692245,-14.58591048),
    "mu_from_eht_midpoint":(0.84615174,4.73532757),
    "loo_mae_ip_eV":0.05645,"loo_mae_ea_eV":0.10998,"loo_mae_mu_eV":0.00959,
    "reference":"offline GFN1-xTB --vip/--vea PAH calibration",
}

def _affine(key,x):
    a,b=PAH_CALIBRATION[key]; return float(a*float(x)+b)

def conjugated_frontier_eht(atoms,K=K_BENZENE_PI):
    syms=set(atoms.get_chemical_symbols())
    if not syms.issubset({"C","H"}): raise ValueError("Current frontier calibration is validated only for conjugated C/H systems")
    nC=sum(s=="C" for s in atoms.get_chemical_symbols())
    if nC<2 or nC%2: raise ValueError("Current pi frontier helper expects an even number of carbon pz electrons")
    e=conjugated_carbon_pi_eht(atoms,K=K); vals=np.asarray(e["energies_eV"],float); nocc=nC//2
    homo,lumo=float(vals[nocc-1]),float(vals[nocc])
    return {"homo_eV":homo,"lumo_eV":lumo,"midpoint_eV":0.5*(homo+lumo),"gap_eV":lumo-homo,
            "n_carbon":int(nC),"overlap_min_eigenvalue":e["overlap_eig_min"]}

def calibrated_pah_redox_from_eht(homo_eV,lumo_eV):
    ip=_affine("ip_from_minus_homo",-float(homo_eV)); ea=_affine("ea_from_minus_lumo",-float(lumo_eV)); midpoint=.5*(float(homo_eV)+float(lumo_eV)); mu=_affine("mu_from_eht_midpoint",midpoint)
    return {"IP_vertical_eV":ip,"EA_vertical_eV":ea,"chemical_potential_proxy_eV":mu,"calibration":dict(PAH_CALIBRATION),
            "confidence":"rough; validated only on small conjugated C/H PAHs"}

def graphene_fermi_proxy_eht(atoms,K=1.75):
    e=graphene_pi_eht(atoms,cutoff_A=1.85,K=K); vals=np.linalg.eigvalsh(e["Heff"]); n=len(vals)
    if n%2: raise ValueError("Graphene pi supercell must contain an even number of carbon sites")
    homo,lumo=float(vals[n//2-1]),float(vals[n//2]); finite_midpoint=.5*(homo+lumo); dirac_eht=float(_p_params("C")[0]); mu_cal=_affine("mu_from_eht_midpoint",dirac_eht)
    return {"homo_proxy_eV":homo,"lumo_proxy_eV":lumo,"finite_supercell_midpoint_eV":finite_midpoint,"mu_EHT_eV":dirac_eht,
            "mu_calibrated_eV":mu_cal,"work_function_proxy_eV":-mu_cal,"finite_supercell_gap_eV":lumo-homo,
            "reference_definition":"carbon p onsite / Dirac level; supercell midpoint is diagnostic only","overlap_min_eigenvalue":e["overlap_eig_min"]}

def molecule_graphene_alignment_eht(graphene_atoms,molecule_atoms,K_molecule=K_BENZENE_PI):
    mf=conjugated_frontier_eht(molecule_atoms,K=K_molecule); mr=calibrated_pah_redox_from_eht(mf["homo_eV"],mf["lumo_eV"]); gr=graphene_fermi_proxy_eht(graphene_atoms)
    W=float(gr["work_function_proxy_eV"]); hole_offset=float(mr["IP_vertical_eV"]-W); electron_offset=float(W-mr["EA_vertical_eV"]); delta_mu=float(mr["chemical_potential_proxy_eV"]-gr["mu_calibrated_eV"])
    return {"molecule_frontier_EHT":mf,"molecule_redox_proxy":mr,"graphene_reference":gr,"delta_mu_molecule_minus_graphene_eV":delta_mu,
            "hole_localization_offset_eV":hole_offset,"electron_localization_offset_eV":electron_offset,
            "interpretation":{"positive_offset":"molecular localization unfavorable relative to graphene",
                              "missing_physics":["image-charge polarization","environment screening","interface dipole","coverage-dependent electrostatics"]},
            "runtime_dependencies":"Python/NumPy/ASE only; xTB was used offline to fit calibration constants"}
