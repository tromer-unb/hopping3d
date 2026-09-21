#!/usr/bin/env python3
"""
hopping3d/core.py

Trajectory-resolved kinetic Monte Carlo (KMC) transport on finite 3D atomic graphs.
Designed as a reusable atomistic hopping-transport engine.

Key differences implemented in this release:
  * true 3D Cartesian graph with arbitrary transport direction vector;
  * cKDTree neighbor construction (instead of O(N^2));
  * chemistry-aware site offsets and pair-dependent coupling scales;
  * inter-layer coupling control for 2D -> quasi-2D -> 3D crossover;
  * legacy global gap activation kept only as an explicit ablation option;
  * trajectory observables: first-passage time, hops, path length, tortuosity;
  * periodic unbiased KMC for a diffusion/mobility tensor;
  * Landauer-like current retained only as a backwards-compatibility diagnostic.

The code intentionally separates geometry/topology from electronic parameters so that
DFT/Wannier-derived site energies and transfer integrals can replace the provisional
chemistry parameters without changing the KMC engine.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ase import Atoms
from ase.io import read, write
from ase.neighborlist import neighbor_list
from scipy.spatial import cKDTree

KB_EV = 8.617333262145e-5
KB_J = 1.380649e-23
E_CHARGE = 1.602176634e-19
PLANCK = 6.62607015e-34
G0 = 2.0 * E_CHARGE**2 / PLANCK
ANGSTROM = 1.0e-10

IONIZATION_POTENTIAL_EV = {
    "H": 13.598, "Li": 5.392, "Na": 5.139, "K": 4.341,
    "Rb": 4.177, "Cs": 3.894, "Be": 9.323, "Mg": 7.646,
    "Ca": 6.113, "B": 8.298, "C": 11.260, "N": 14.534,
    "O": 13.618, "F": 17.423, "Al": 5.986, "Si": 8.152,
    "P": 10.487, "S": 10.360, "Cl": 12.968, "Ti": 6.828,
    "V": 6.746, "Cr": 6.767, "Mn": 7.434, "Fe": 7.902,
    "Co": 7.881, "Ni": 7.640, "Cu": 7.726, "Zn": 9.394,
    "W": 7.864, "Tl": 6.108,
}


@dataclass
class Edge:
    j: int
    r: float
    dr: np.ndarray


def load_json(path: str | Path) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def save_json(path: str | Path, obj: Any) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, sort_keys=True)


def unit_vector(v: Sequence[float]) -> np.ndarray:
    a = np.asarray(v, dtype=float)
    n = float(np.linalg.norm(a))
    if n <= 0.0:
        raise ValueError("Direction vector must be non-zero")
    return a / n


def direction_vector(spec: Any, atoms: Optional[Atoms] = None) -> np.ndarray:
    if isinstance(spec, str):
        s = spec.strip().upper()
        if s == "X": return np.array([1.0, 0.0, 0.0])
        if s == "Y": return np.array([0.0, 1.0, 0.0])
        if s == "Z": return np.array([0.0, 0.0, 1.0])
        if s in ("A", "B", "C"):
            if atoms is None:
                raise ValueError("A/B/C direction requires atoms")
            ax = ["A", "B", "C"].index(s)
            return unit_vector(np.asarray(atoms.cell.reciprocal()[ax], dtype=float))
        if s in ("A_DIRECT", "B_DIRECT", "C_DIRECT"):
            if atoms is None:
                raise ValueError("*_DIRECT direction requires atoms")
            ax = ["A_DIRECT", "B_DIRECT", "C_DIRECT"].index(s)
            return unit_vector(np.asarray(atoms.cell[ax], dtype=float))
        vals = [float(x) for x in s.replace(",", " ").split()]
        if len(vals) != 3:
            raise ValueError(f"Cannot parse direction: {spec}")
        return unit_vector(vals)
    if isinstance(spec, (list, tuple, np.ndarray)) and len(spec) == 3:
        return unit_vector(spec)
    raise ValueError(f"Unsupported direction specification: {spec!r}")


def load_structure(cif: str | Path, repeat: Sequence[int] = (1, 1, 1)) -> Atoms:
    atoms = read(str(cif))
    rep = tuple(int(x) for x in repeat)
    if rep != (1, 1, 1):
        atoms = atoms.repeat(rep)
    return atoms


def apply_random_vacancies(atoms: Atoms, fraction: float, seed: int = 0,
                           allowed_symbols: Optional[Sequence[str]] = None) -> Atoms:
    fraction = float(fraction)
    if fraction <= 0.0:
        return atoms.copy()
    if fraction >= 1.0:
        raise ValueError("vacancy fraction must be < 1")
    symbols = np.asarray(atoms.get_chemical_symbols())
    if allowed_symbols:
        mask = np.isin(symbols, np.asarray(list(allowed_symbols), dtype=str))
        candidates = np.where(mask)[0]
    else:
        candidates = np.arange(len(atoms), dtype=int)
    nremove = int(round(fraction * len(candidates)))
    if nremove <= 0:
        return atoms.copy()
    rng = np.random.default_rng(int(seed))
    remove = set(rng.choice(candidates, size=nremove, replace=False).tolist())
    keep = [i for i in range(len(atoms)) if i not in remove]
    return atoms[keep]


def build_finite_graph(atoms: Atoms, cutoff_A: float) -> List[List[Edge]]:
    pos = np.asarray(atoms.positions, dtype=float)
    tree = cKDTree(pos)
    pairs = tree.query_pairs(float(cutoff_A), output_type="ndarray")
    neigh: List[List[Edge]] = [[] for _ in range(len(atoms))]
    for i, j in pairs:
        dr = pos[j] - pos[i]
        r = float(np.linalg.norm(dr))
        if r <= 1e-14:
            continue
        neigh[int(i)].append(Edge(int(j), r, dr.copy()))
        neigh[int(j)].append(Edge(int(i), r, -dr.copy()))
    return neigh


def build_crystal_device_graph(atoms: Atoms, cutoff_A: float, open_axis: int | str) -> List[List[Edge]]:
    if isinstance(open_axis, str):
        ax = {"A": 0, "B": 1, "C": 2}[open_axis.strip().upper()]
    else:
        ax = int(open_axis)
    work = atoms.copy()
    pbc = np.array([True, True, True], dtype=bool)
    pbc[ax] = False
    work.set_pbc(pbc)
    i_arr, j_arr, d_arr, D_arr = neighbor_list("ijdD", work, float(cutoff_A))
    neigh: List[List[Edge]] = [[] for _ in range(len(work))]
    for i, j, r, dr in zip(i_arr, j_arr, d_arr, D_arr):
        if float(r) <= 1e-14:
            continue
        neigh[int(i)].append(Edge(int(j), float(r), np.asarray(dr, float)))
    return neigh


def build_periodic_graph(atoms: Atoms, cutoff_A: float) -> List[List[Edge]]:
    i_arr, j_arr, d_arr, D_arr = neighbor_list("ijdD", atoms, float(cutoff_A))
    neigh: List[List[Edge]] = [[] for _ in range(len(atoms))]
    for i, j, r, dr in zip(i_arr, j_arr, d_arr, D_arr):
        if float(r) <= 1e-14:
            continue
        neigh[int(i)].append(Edge(int(j), float(r), np.asarray(dr, float)))
    return neigh


def pair_key(a: str, b: str) -> str:
    return "-".join(sorted((str(a), str(b))))


def species_site_offsets(atoms: Atoms, table: Dict[str, float]) -> np.ndarray:
    symbols = atoms.get_chemical_symbols()
    return np.asarray([float(table.get(s, 0.0)) for s in symbols], dtype=float)


def chemistry_site_offsets(atoms: Atoms, chemistry_cfg: Dict[str, Any]) -> np.ndarray:
    arr = species_site_offsets(atoms, chemistry_cfg.get("site_offsets_eV", {}) or {})
    by_index = chemistry_cfg.get("site_offsets_by_index_eV", {}) or {}
    for k, v in by_index.items():
        idx = int(k)
        if idx < 0 or idx >= len(arr):
            raise IndexError(f"site offset index {idx} outside 0..{len(arr)-1}")
        arr[idx] += float(v)
    full = chemistry_cfg.get("site_offsets_array_eV", None)
    if full is not None:
        full = np.asarray(full, dtype=float)
        if len(full) != len(arr):
            raise ValueError("site_offsets_array_eV length must equal number of atoms")
        arr += full
    return arr


def legacy_impurity_potential(atoms: Atoms, cfg: Dict[str, Any]) -> np.ndarray:
    pot = np.zeros(len(atoms), dtype=float)
    if not cfg or not bool(cfg.get("enabled", False)):
        return pot
    symbols = np.asarray(atoms.get_chemical_symbols())
    pos = np.asarray(atoms.positions)
    ref = str(cfg.get("reference", "O"))
    if ref not in IONIZATION_POTENTIAL_EV:
        raise ValueError(f"No IP value for reference element {ref}")
    ip_ref = IONIZATION_POTENTIAL_EV[ref]
    scale = float(cfg.get("scale", 0.2))
    species_cfg = cfg.get("species", {})
    for sp, scfg in species_cfg.items():
        if sp not in IONIZATION_POTENTIAL_EV:
            raise ValueError(f"No IP value for impurity element {sp}")
        centers = pos[symbols == sp]
        if len(centers) == 0:
            continue
        delta = -scale * (IONIZATION_POTENTIAL_EV[sp] - ip_ref)
        eta = float(scfg.get("eta_A", 2.0))
        radius = float(scfg.get("radius_A", 4.0))
        for c in centers:
            r = np.linalg.norm(pos - c, axis=1)
            pot += np.where(r <= radius, delta * np.exp(-r / max(eta, 1e-12)), 0.0)
    return pot


def bias_energy(atoms: Atoms, voltage_V: float, direction: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]:
    pos = np.asarray(atoms.positions)
    u = pos @ direction
    umin, umax = float(u.min()), float(u.max())
    L = max(umax - umin, 1e-12)
    e = -float(voltage_V) * (u - umin) / L
    return e, u, L


def localization_length_A(gap_eV: float, T_K: float, xi0_A: float,
                          alpha_gap: float, beta_T: float,
                          B_T: float = 0.0, B0_T: float = 10.0) -> float:
    xi = float(xi0_A) / (1.0 + float(alpha_gap) * max(float(gap_eV), 0.0))
    xi *= 1.0 + float(beta_T) * (float(T_K) / 300.0 - 1.0)
    if float(B0_T) > 0.0:
        xi /= math.sqrt(1.0 + (float(B_T) / float(B0_T))**2)
    return max(xi, 1e-12)


def edge_coupling_weight(i: int, edge: Edge, sym_i: str, sym_j: str, cfg: Dict[str, Any]) -> float:
    w = 1.0
    pair_scales = cfg.get("pair_scales", {}) or {}
    w *= float(pair_scales.get(pair_key(sym_i, sym_j), 1.0))
    jmap = cfg.get("edge_transfer_integrals_eV", {}) or {}
    if jmap:
        key = f"{min(int(i), int(edge.j))}-{max(int(i), int(edge.j))}"
        if key in jmap:
            J = abs(float(jmap[key]))
            Jref = max(abs(float(cfg.get("J_ref_eV", 1.0))), 1e-15)
            w *= (J / Jref) ** 2
    emap = cfg.get("edge_scales_by_index", {}) or {}
    if emap:
        key = f"{min(int(i), int(edge.j))}-{max(int(i), int(edge.j))}"
        if key in emap:
            w *= max(float(emap[key]), 0.0)
    if "interlayer_scale" in cfg or "intralayer_cutoff_A" in cfg or "interlayer_cutoff_A" in cfg:
        scale = max(float(cfg.get("interlayer_scale", 1.0)), 0.0)
        normal = unit_vector(cfg.get("layer_normal", [0.0, 0.0, 1.0]))
        threshold = float(cfg.get("interlayer_threshold_A", 1.0))
        dz = abs(float(np.dot(edge.dr, normal)))
        is_inter = dz > threshold
        if is_inter:
            if edge.r > float(cfg.get("interlayer_cutoff_A", 1.0e9)):
                return 0.0
            w *= scale
        else:
            if edge.r > float(cfg.get("intralayer_cutoff_A", 1.0e9)):
                return 0.0
    tensor = cfg.get("orientation_tensor", None)
    if tensor is not None:
        M = np.asarray(tensor, dtype=float).reshape(3, 3)
        u = edge.dr / max(edge.r, 1e-12)
        orient = float(u @ M @ u)
        w *= max(orient, 0.0)
    return max(float(w), 0.0)


def transition_rates(i: int, neighbors_i: List[Edge], energies_eV: np.ndarray,
                     symbols: Sequence[str], cfg: Dict[str, Any], voltage_V: float) -> Tuple[np.ndarray, np.ndarray]:
    T = float(cfg.get("temperature_K", 300.0))
    gap = float(cfg.get("gap_eV", 0.0))
    xi = localization_length_A(gap, T, float(cfg.get("xi0_A", 3.0)),
        float(cfg.get("alpha_gap", 0.0)), float(cfg.get("beta_T", 0.0)),
        float(cfg.get("B_T", 0.0)), float(cfg.get("B0_T", 10.0)))
    nu0 = float(cfg.get("nu0_Hz", 1.0))
    gap_mode = str(cfg.get("gap_activation", "none")).lower()
    f_gap = 1.0
    if gap_mode == "legacy_global" and T > 0.0:
        gamma = float(cfg.get("gap_gamma", 1.0))
        barrier = max(0.0, gap - gamma * abs(float(voltage_V)))
        f_gap = math.exp(-barrier / (KB_EV * T))
    elif gap_mode not in ("none", "off", "false"):
        raise ValueError(f"Unknown gap_activation mode: {gap_mode}")
    js: List[int] = []
    rs: List[float] = []
    for edge in neighbors_i:
        dE = float(energies_eV[edge.j] - energies_eV[i])
        if T > 0.0:
            thermal = math.exp(-dE / (KB_EV * T)) if dE > 0.0 else 1.0
        else:
            thermal = 1.0 if dE <= 0.0 else 0.0
        coupling = edge_coupling_weight(i, edge, symbols[i], symbols[edge.j], cfg)
        rate = nu0 * coupling * math.exp(-2.0 * edge.r / xi) * thermal * f_gap
        if rate > 0.0 and np.isfinite(rate):
            js.append(edge.j)
            rs.append(rate)
    return np.asarray(js, dtype=int), np.asarray(rs, dtype=float)


def prepare_rate_table(neighbors: List[List[Edge]], energies_eV: np.ndarray,
                       symbols: Sequence[str], cfg: Dict[str, Any], voltage_V: float):
    table = []
    for i, edges in enumerate(neighbors):
        js, rates = transition_rates(i, edges, energies_eV, symbols, cfg, voltage_V)
        dmap = {e.j: e.r for e in edges}
        lengths = np.asarray([dmap[int(j)] for j in js], dtype=float)
        table.append((js, rates, lengths))
    return table


def electrode_indices(u: np.ndarray, width_A: Optional[float] = None,
                      width_fraction: Optional[float] = None) -> Tuple[np.ndarray, float, float, float]:
    umin, umax = float(u.min()), float(u.max())
    L = max(umax - umin, 1e-12)
    if width_A is None:
        width_A = L * float(width_fraction if width_fraction is not None else 0.05)
    width = min(max(float(width_A), 1e-8), 0.45 * L)
    source = np.where(u <= umin + width)[0]
    if len(source) == 0:
        source = np.asarray([int(np.argmin(u))])
    return source, umin, umax, width


def simulate_first_passage(atoms: Atoms, neighbors: List[List[Edge]], voltage_V: float,
                           transport_cfg: Dict[str, Any], chemistry_cfg: Optional[Dict[str, Any]] = None,
                           rng_seed: int = 0, record_paths: int = 0) -> Dict[str, Any]:
    chemistry_cfg = chemistry_cfg or {}
    direction = direction_vector(transport_cfg.get("direction", "X"), atoms)
    ebias, u, L_A = bias_energy(atoms, voltage_V, direction)
    offsets = chemistry_site_offsets(atoms, chemistry_cfg)
    imp = legacy_impurity_potential(atoms, chemistry_cfg.get("legacy_impurity", {}) or {})
    energies = ebias + offsets + imp
    symbols = atoms.get_chemical_symbols()
    rate_table = prepare_rate_table(neighbors, energies, symbols, transport_cfg, voltage_V)
    source, umin, umax, width = electrode_indices(u, transport_cfg.get("electrode_width_A", None),
        transport_cfg.get("electrode_width_fraction", 0.05))
    drain_threshold = umax - width
    ntraj = int(transport_cfg.get("n_trajectories", 500))
    max_steps = int(transport_cfg.get("max_steps", 5000))
    max_time = float(transport_cfg.get("max_time_s", 1.0))
    rng = np.random.default_rng(int(rng_seed))
    passed = 0
    success_times: List[float] = []
    success_hops: List[int] = []
    success_lengths: List[float] = []
    success_tortuosity: List[float] = []
    saved_paths: List[Dict[str, Any]] = []
    pos = np.asarray(atoms.positions)
    for itraj in range(ntraj):
        i = int(rng.choice(source)); start = i; t = 0.0; hops = 0; path_length = 0.0
        path = [i] if len(saved_paths) < int(record_paths) else None
        times = [0.0] if path is not None else None
        ok = bool(u[i] >= drain_threshold)
        while (not ok) and hops < max_steps and t < max_time:
            js, rates, lengths = rate_table[i]
            if len(js) == 0:
                break
            total = float(rates.sum())
            dt = -math.log(max(rng.random(), 1e-15)) / total
            if not np.isfinite(dt) or t + dt > max_time:
                break
            xsel = rng.random() * total
            ksel = int(np.searchsorted(np.cumsum(rates), xsel, side="right"))
            ksel = min(ksel, len(js) - 1)
            j = int(js[ksel]); path_length += float(lengths[ksel]); t += dt; i = j; hops += 1
            if path is not None:
                path.append(i); times.append(t)
            if u[i] >= drain_threshold:
                ok = True
        if ok:
            passed += 1
            success_times.append(t); success_hops.append(hops); success_lengths.append(path_length)
            straight = float(np.linalg.norm(pos[i] - pos[start]))
            success_tortuosity.append(path_length / max(straight, 1e-12))
            if path is not None:
                saved_paths.append({"indices": path, "times_s": times})
    T_eff = passed / max(ntraj, 1)
    stderr = math.sqrt(max(T_eff * (1.0 - T_eff) / max(ntraj, 1), 0.0))
    mt = float(np.mean(success_times)) if success_times else math.nan
    mh = float(np.mean(success_hops)) if success_hops else math.nan
    mpl = float(np.mean(success_lengths)) if success_lengths else math.nan
    mtor = float(np.mean(success_tortuosity)) if success_tortuosity else math.nan
    drift_v = math.nan; mobility = math.nan
    if success_times and voltage_V != 0.0 and L_A > 0.0:
        L_m = L_A * ANGSTROM; drift_v = L_m / mt; field = abs(float(voltage_V)) / L_m
        mobility = drift_v / field if field > 0.0 else math.nan
    G = G0 * T_eff; I = G * float(voltage_V)
    return {"V": float(voltage_V), "T_eff": T_eff, "T_stderr": stderr,
        "n_pass": int(passed), "n_total": int(ntraj), "mean_first_passage_s": mt,
        "mean_hops": mh, "mean_path_length_A": mpl, "mean_tortuosity": mtor,
        "drift_velocity_m_s": drift_v, "mobility_drift_m2_Vs": mobility,
        "G_landauer_compat_S": G, "I_landauer_compat_A": I,
        "device_length_A": L_A, "electrode_width_A": width,
        "direction_vector": direction.tolist(), "paths": saved_paths}


def run_voltage_sweep(atoms: Atoms, neighbors: List[List[Edge]], voltages: Iterable[float],
                      transport_cfg: Dict[str, Any], chemistry_cfg: Optional[Dict[str, Any]] = None,
                      seed: int = 0) -> pd.DataFrame:
    rows = []
    for k, V in enumerate(voltages):
        out = simulate_first_passage(atoms, neighbors, float(V), transport_cfg,
                                     chemistry_cfg=chemistry_cfg,
                                     rng_seed=int(seed) + 1009 * k, record_paths=0)
        out.pop("paths", None); rows.append(out)
    return pd.DataFrame(rows)


def plot_sweep(df: pd.DataFrame, out_png: str | Path, title: str = "") -> None:
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.5))
    ax = axes[0, 0]; ax.errorbar(df["V"], df["T_eff"], yerr=1.96 * df["T_stderr"], marker="o", capsize=3)
    ax.set(xlabel="Bias V (V)", ylabel="Effective transmittance", ylim=(-0.03, 1.03)); ax.grid(alpha=0.25)
    ax = axes[0, 1]; ax.plot(df["V"], df["mean_first_passage_s"], "o-")
    ax.set(xlabel="Bias V (V)", ylabel="Mean first-passage time (s)"); ax.set_yscale("log"); ax.grid(alpha=0.25)
    ax = axes[1, 0]; ax.plot(df["V"], df["mean_tortuosity"], "o-")
    ax.set(xlabel="Bias V (V)", ylabel="Mean tortuosity"); ax.grid(alpha=0.25)
    ax = axes[1, 1]; ax.plot(df["V"], df["mobility_drift_m2_Vs"], "o-")
    ax.set(xlabel="Bias V (V)", ylabel="Drift mobility proxy (m$^2$/V s)"); ax.grid(alpha=0.25)
    if title: fig.suptitle(title)
    fig.tight_layout(); fig.savefig(out_png, dpi=180); plt.close(fig)


def plot_structure_3d(atoms: Atoms, out_png: str | Path, title: str = "Structure") -> None:
    from ase.data import atomic_numbers, covalent_radii
    pos = np.asarray(atoms.positions); syms = np.asarray(atoms.get_chemical_symbols())
    fig = plt.figure(figsize=(7, 6)); ax = fig.add_subplot(111, projection="3d")
    for sp in sorted(set(syms.tolist())):
        idx = np.where(syms == sp)[0]; z = atomic_numbers[sp]
        size = max(float(covalent_radii[z]) * 35.0, 18.0)
        ax.scatter(pos[idx, 0], pos[idx, 1], pos[idx, 2], s=size, label=sp, alpha=0.82)
    ax.set(xlabel="x (Å)", ylabel="y (Å)", zlabel="z (Å)", title=title); ax.legend(loc="best", fontsize=8)
    fig.tight_layout(); fig.savefig(out_png, dpi=180); plt.close(fig)


def plot_paths_3d(atoms: Atoms, paths: List[Dict[str, Any]], out_png: str | Path,
                  title: str = "Successful KMC trajectories") -> None:
    pos = np.asarray(atoms.positions); fig = plt.figure(figsize=(7, 6)); ax = fig.add_subplot(111, projection="3d")
    ax.scatter(pos[:, 0], pos[:, 1], pos[:, 2], s=5, alpha=0.15)
    for p in paths:
        idx = np.asarray(p["indices"], dtype=int); ax.plot(pos[idx, 0], pos[idx, 1], pos[idx, 2], alpha=0.35, lw=0.8)
    ax.set(xlabel="x (Å)", ylabel="y (Å)", zlabel="z (Å)", title=title)
    fig.tight_layout(); fig.savefig(out_png, dpi=180); plt.close(fig)


def simulate_periodic_diffusion(atoms: Atoms, cutoff_A: float, transport_cfg: Dict[str, Any],
                                chemistry_cfg: Optional[Dict[str, Any]] = None,
                                n_walkers: int = 1000, observation_time_s: float = 1e-10,
                                seed: int = 0) -> Dict[str, Any]:
    chemistry_cfg = chemistry_cfg or {}
    neigh = build_periodic_graph(atoms, cutoff_A)
    offsets = chemistry_site_offsets(atoms, chemistry_cfg)
    offsets += legacy_impurity_potential(atoms, chemistry_cfg.get("legacy_impurity", {}) or {})
    symbols = atoms.get_chemical_symbols(); cfg = dict(transport_cfg); cfg["gap_activation"] = "none"
    ptable = []
    for i, edges in enumerate(neigh):
        kept_edges = []; kept_rates = []
        for edge in edges:
            _js, _rates = transition_rates(i, [edge], offsets, symbols, cfg, voltage_V=0.0)
            if len(_rates): kept_edges.append(edge); kept_rates.append(float(_rates[0]))
        ptable.append((kept_edges, np.asarray(kept_rates, dtype=float)))
    T = float(cfg.get("temperature_K", 300.0)); rng = np.random.default_rng(int(seed))
    disps = np.zeros((int(n_walkers), 3), dtype=float); hop_counts = np.zeros(int(n_walkers), dtype=int)
    init_mode = str(cfg.get("initial_distribution", "equilibrium")).lower()
    if init_mode in ("equilibrium", "boltzmann") and T > 0.0:
        e0 = offsets - float(np.min(offsets)); weights = np.exp(-np.clip(e0 / (KB_EV * T), 0.0, 700.0)); start_prob = weights / weights.sum()
    elif init_mode in ("uniform", "legacy"):
        start_prob = None
    else:
        raise ValueError(f"Unknown initial_distribution: {init_mode}")
    for w in range(int(n_walkers)):
        i = int(rng.integers(0, len(atoms))) if start_prob is None else int(rng.choice(len(atoms), p=start_prob))
        t = 0.0; d = np.zeros(3, dtype=float)
        while t < observation_time_s:
            edges, rates = ptable[i]
            if len(rates) == 0: break
            total = float(rates.sum()); dt = -math.log(max(rng.random(), 1e-15)) / total
            if t + dt > observation_time_s: break
            x = rng.random() * total; k = int(np.searchsorted(np.cumsum(rates), x, side="right")); k = min(k, len(edges) - 1)
            edge = edges[k]; d += edge.dr; i = int(edge.j); t += dt; hop_counts[w] += 1
        disps[w] = d
    outer_mean = np.mean(np.einsum("ni,nj->nij", disps, disps), axis=0)
    D_A2_s = outer_mean / (2.0 * float(observation_time_s)); D = D_A2_s * ANGSTROM**2; mu = E_CHARGE * D / (KB_J * T)
    return {"D_m2_s": D.tolist(), "mu_m2_Vs": mu.tolist(), "mean_hops": float(np.mean(hop_counts)),
        "median_hops": float(np.median(hop_counts)), "zero_hop_fraction": float(np.mean(hop_counts == 0)),
        "initial_distribution": init_mode, "n_walkers": int(n_walkers), "observation_time_s": float(observation_time_s),
        "mean_squared_displacement_A2": float(np.mean(np.sum(disps**2, axis=1)))}


def save_tensor_csv(path: str | Path, tensor: Sequence[Sequence[float]], prefix: str) -> None:
    arr = np.asarray(tensor, dtype=float); labels = [f"{prefix}_x", f"{prefix}_y", f"{prefix}_z"]
    pd.DataFrame(arr, index=labels, columns=labels).to_csv(path)


def prepare_from_config(config_path: str | Path) -> Tuple[Dict[str, Any], Path, Atoms, List[List[Edge]]]:
    config_path = Path(config_path).resolve(); cfg = load_json(config_path); base = config_path.parent
    scfg = cfg.get("structure", {}); cif = Path(scfg["cif"])
    if not cif.is_absolute(): cif = base / cif
    atoms = load_structure(cif, scfg.get("repeat", [1, 1, 1]))
    atoms = apply_random_vacancies(atoms, scfg.get("vacancy_fraction", 0.0), scfg.get("vacancy_seed", 0), scfg.get("vacancy_symbols", None))
    cutoff = float(cfg.get("graph", {}).get("cutoff_A", 3.0)); graph = build_finite_graph(atoms, cutoff)
    return cfg, base, atoms, graph


def run_from_config(config_path: str | Path) -> pd.DataFrame:
    cfg, base, atoms, graph = prepare_from_config(config_path)
    outdir = base / cfg.get("output", {}).get("directory", "results"); outdir.mkdir(parents=True, exist_ok=True)
    write(outdir / "structure_used.cif", atoms); plot_structure_3d(atoms, outdir / "structure.png", cfg.get("output", {}).get("title", "Structure"))
    tcfg = dict(cfg.get("transport", {})); tcfg.update(cfg.get("graph", {})); chem = cfg.get("chemistry", {})
    volts = tcfg.get("voltages_V", [0.0, 0.5, 1.0]); seed = int(cfg.get("random_seed", 1234))
    df = run_voltage_sweep(atoms, graph, volts, tcfg, chem, seed); df.to_csv(outdir / "sweep.csv", index=False)
    plot_sweep(df, outdir / "sweep_summary.png", cfg.get("output", {}).get("title", ""))
    path_V = float(cfg.get("output", {}).get("path_voltage_V", max(volts))); npaths = int(cfg.get("output", {}).get("record_paths", 40))
    pout = simulate_first_passage(atoms, graph, path_V, tcfg, chem, rng_seed=seed + 99991, record_paths=npaths)
    plot_paths_3d(atoms, pout["paths"], outdir / "paths_3d.png", f"Successful trajectories at V={path_V:g} V")
    meta = {k: v for k, v in pout.items() if k != "paths"}; meta["n_atoms"] = len(atoms); meta["mean_degree"] = float(np.mean([len(n) for n in graph]))
    save_json(outdir / "path_metrics.json", meta); return df
