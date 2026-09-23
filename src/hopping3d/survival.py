#!/usr/bin/env python3
"""First-passage diagnostics with administrative censoring and adaptive Tmax.

This module is intentionally separate from the production engine.  It records one row
per injected trajectory, distinguishes arrival from numerical censoring, and can increase
the observation window until the end-of-window censored fraction is small and T stabilizes.
"""
from __future__ import annotations
import math
from typing import Any, Dict, Optional, List
import numpy as np
import pandas as pd

from .core import (
    ANGSTROM, G0, bias_energy, chemistry_site_offsets, direction_vector,
    electrode_indices, legacy_impurity_potential, prepare_rate_table,
)


def simulate_detailed(atoms, neighbors, voltage_V: float, transport_cfg: Dict[str, Any],
                      chemistry_cfg: Optional[Dict[str, Any]] = None,
                      rng_seed: int = 0, record_paths: int = 0) -> Dict[str, Any]:
    """KMC first-passage simulation with per-trajectory termination records.

    Each trajectory receives its own deterministic RNG stream. Therefore rerunning with
    the same rng_seed but a larger Tmax preserves each trajectory up to the old censor time.
    """
    chemistry_cfg = chemistry_cfg or {}
    direction = direction_vector(transport_cfg.get('direction', 'X'), atoms)
    ebias, u, L_A = bias_energy(atoms, voltage_V, direction)
    offsets = chemistry_site_offsets(atoms, chemistry_cfg)
    imp = legacy_impurity_potential(atoms, chemistry_cfg.get('legacy_impurity', {}) or {})
    energies = ebias + offsets + imp
    symbols = atoms.get_chemical_symbols()
    rate_table = prepare_rate_table(neighbors, energies, symbols, transport_cfg, voltage_V)

    source, umin, umax, width = electrode_indices(
        u, transport_cfg.get('electrode_width_A', None),
        transport_cfg.get('electrode_width_fraction', 0.05))
    override = transport_cfg.get('source_indices_override', None)
    if override is not None:
        allowed = set(int(x) for x in override)
        source = np.asarray([int(i) for i in source if int(i) in allowed], dtype=int)
        if len(source) == 0:
            raise ValueError('source_indices_override leaves no source-electrode sites')
    drain_threshold = umax - width
    ntraj = int(transport_cfg.get('n_trajectories', 500))
    max_steps = int(transport_cfg.get('max_steps', 5000))
    max_time = float(transport_cfg.get('max_time_s', 1.0))
    pos = np.asarray(atoms.positions, dtype=float)

    records: List[Dict[str, Any]] = []
    saved_paths = []
    for itraj in range(ntraj):
        rng = np.random.default_rng(int(rng_seed) + 104729 * itraj)
        i = int(rng.choice(source)); start = i
        t = 0.0; hops = 0; path_length = 0.0
        keep_path = len(saved_paths) < int(record_paths)
        path = [i] if keep_path else None
        times = [0.0] if keep_path else None
        arrived = bool(u[i] >= drain_threshold)
        reason = 'arrived' if arrived else None

        while (not arrived) and hops < max_steps and t < max_time:
            js, rates, lengths = rate_table[i]
            if len(js) == 0 or float(rates.sum()) <= 0.0:
                reason = 'no_rates'; break
            total = float(rates.sum())
            dt = -math.log(max(rng.random(), 1e-15)) / total
            if (not np.isfinite(dt)) or t + dt > max_time:
                t = max_time; reason = 'max_time'; break
            xsel = rng.random() * total
            ksel = int(np.searchsorted(np.cumsum(rates), xsel, side='right'))
            ksel = min(ksel, len(js)-1)
            i = int(js[ksel]); path_length += float(lengths[ksel]); t += dt; hops += 1
            if keep_path:
                path.append(i); times.append(t)
            if u[i] >= drain_threshold:
                arrived = True; reason = 'arrived'

        if not arrived and reason is None:
            reason = 'max_steps' if hops >= max_steps else 'max_time'
        straight = float(np.linalg.norm(pos[i] - pos[start]))
        records.append(dict(
            trajectory=itraj, event_observed=int(arrived), arrival_time_s=t if arrived else np.nan,
            terminal_time_s=t, termination=reason, hops=hops, path_length_A=path_length,
            tortuosity=path_length/max(straight,1e-12) if arrived else np.nan,
            start_index=start, end_index=i))
        if arrived and keep_path:
            saved_paths.append({'indices': path, 'times_s': times})

    rec = pd.DataFrame(records)
    T_eff = float(rec.event_observed.mean())
    stderr = math.sqrt(max(T_eff*(1-T_eff)/max(ntraj,1),0.0))
    succ = rec[rec.event_observed == 1]
    mt = float(succ.arrival_time_s.mean()) if len(succ) else math.nan
    L_m = L_A * ANGSTROM
    mobility = math.nan; drift_v = math.nan
    if len(succ) and voltage_V != 0.0 and L_m > 0.0:
        drift_v = L_m / mt
        field = abs(float(voltage_V)) / L_m
        mobility = drift_v / field if field > 0 else math.nan
    return dict(
        V=float(voltage_V), T_eff=T_eff, T_stderr=stderr, n_pass=int(rec.event_observed.sum()),
        n_total=ntraj, mean_first_passage_s=mt,
        mean_hops=float(succ.hops.mean()) if len(succ) else math.nan,
        mean_path_length_A=float(succ.path_length_A.mean()) if len(succ) else math.nan,
        mean_tortuosity=float(succ.tortuosity.mean()) if len(succ) else math.nan,
        drift_velocity_m_s=drift_v, mobility_drift_m2_Vs=mobility,
        G_landauer_compat_S=G0*T_eff, I_landauer_compat_A=G0*T_eff*float(voltage_V),
        device_length_A=L_A, electrode_width_A=width, direction_vector=direction.tolist(),
        n_source_sites=int(len(source)),
        censored_fraction=float((rec.termination == 'max_time').mean()),
        max_steps_fraction=float((rec.termination == 'max_steps').mean()),
        trapped_fraction=float((rec.termination == 'no_rates').mean()),
        records=rec, paths=saved_paths)


def arrival_curve(records: pd.DataFrame, n_grid: int = 160) -> pd.DataFrame:
    """Empirical first-arrival CDF and survival curve under common administrative censoring."""
    if len(records) == 0:
        return pd.DataFrame(columns=['time_s','F_arrival','S_not_arrived'])
    tmax = float(records.terminal_time_s.max())
    grid = np.linspace(0.0, tmax, int(n_grid))
    arrivals = records.loc[records.event_observed == 1, 'arrival_time_s'].to_numpy(float)
    F = np.asarray([(arrivals <= t).sum()/len(records) for t in grid], dtype=float)
    return pd.DataFrame({'time_s': grid, 'F_arrival': F, 'S_not_arrived': 1.0-F})


def adaptive_first_passage(atoms, neighbors, voltage_V: float, transport_cfg: Dict[str, Any],
                           chemistry_cfg: Optional[Dict[str, Any]] = None, rng_seed: int = 0,
                           initial_time_s: float = 2e-11, max_time_s: float = 2e-9,
                           growth: float = 2.0, delta_T_tol: float = 0.01,
                           censor_tol: float = 0.01) -> Dict[str, Any]:
    """Increase Tmax until numerical right-censoring is small and T has stabilized."""
    history=[]; previous=None; t=float(initial_time_s); last=None
    while True:
        cfg=dict(transport_cfg); cfg['max_time_s']=t
        out=simulate_detailed(atoms,neighbors,voltage_V,cfg,chemistry_cfg,rng_seed,record_paths=0)
        delta = math.nan if previous is None else abs(out['T_eff']-previous)
        history.append(dict(max_time_s=t,T_eff=out['T_eff'],T_stderr=out['T_stderr'],
                            censored_fraction=out['censored_fraction'],trapped_fraction=out['trapped_fraction'],
                            max_steps_fraction=out['max_steps_fraction'],delta_T=delta,
                            mean_first_passage_s=out['mean_first_passage_s']))
        last=out
        stable = previous is not None and delta <= float(delta_T_tol)
        low_censor = out['censored_fraction'] <= float(censor_tol)
        if stable and low_censor: break
        if t >= float(max_time_s): break
        previous=out['T_eff']; t=min(float(max_time_s),t*float(growth))
    return {'history': pd.DataFrame(history), 'final': last}
