from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .core import (
    apply_random_vacancies,
    build_finite_graph,
    load_json,
    load_structure,
    run_from_config,
    save_json,
    save_tensor_csv,
    simulate_periodic_diffusion,
)
from .graph import source_drain_reachability
from .tensor import principal_decomposition


def _structure_path(base: Path, scfg: dict) -> Path:
    """Return a structure path, accepting the preferred `file` key and legacy `cif`."""
    raw = scfg.get("file", scfg.get("cif"))
    if raw is None:
        raise KeyError("structure.file is required (legacy structure.cif is also accepted)")
    path = Path(raw)
    return path if path.is_absolute() else base / path


def _resolve_atoms(config_path: Path, cfg: dict):
    base = config_path.resolve().parent
    scfg = cfg.get("structure", {})
    structure = _structure_path(base, scfg)
    atoms = load_structure(structure, scfg.get("repeat", [1, 1, 1]))
    atoms = apply_random_vacancies(
        atoms,
        scfg.get("vacancy_fraction", 0.0),
        scfg.get("vacancy_seed", 0),
        scfg.get("vacancy_symbols", None),
    )
    return base, atoms


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="hopping3d",
        description=(
            "Trajectory-resolved hopping KMC from atomistic structure files. "
            "Use calculation=first_passage, diffusion_tensor, or reachability."
        ),
    )
    p.add_argument("config", help="JSON configuration file")
    args = p.parse_args(argv)

    cfg = load_json(args.config)
    mode = str(cfg.get("calculation", "first_passage")).lower()

    if mode in {"first_passage", "device"}:
        df = run_from_config(args.config)
        print(df.to_string(index=False))
        return 0

    if mode in {"diffusion", "diffusion_tensor", "periodic_diffusion"}:
        config_path = Path(args.config)
        base, atoms = _resolve_atoms(config_path, cfg)
        graph_cfg = dict(cfg.get("graph", {}))
        tcfg = dict(cfg.get("transport", {}))
        tcfg.update(graph_cfg)
        chem = cfg.get("chemistry", {})
        dcfg = cfg.get("diffusion", {})
        outdir = base / cfg.get("output", {}).get("directory", "results")
        outdir.mkdir(parents=True, exist_ok=True)
        cutoff = float(graph_cfg.get("cutoff_A", 3.8))
        result = simulate_periodic_diffusion(
            atoms,
            cutoff,
            tcfg,
            chem,
            n_walkers=int(dcfg.get("n_walkers", 500)),
            observation_time_s=float(dcfg.get("observation_time_s", 8e-11)),
            seed=int(cfg.get("random_seed", 1234)),
        )
        D = np.asarray(result["D_m2_s"], float)
        vals, vecs = principal_decomposition(D)
        save_tensor_csv(outdir / "diffusion_tensor.csv", D, "D")
        result["principal_diffusivities_m2_s"] = vals.tolist()
        result["principal_axes"] = vecs.tolist()
        save_json(outdir / "diffusion_summary.json", result)
        print("D tensor (m^2/s):\n", D)
        print("principal diffusivities:", vals)
        return 0

    if mode in {"reachability", "connectivity", "graph_reachability"}:
        config_path = Path(args.config)
        base, atoms = _resolve_atoms(config_path, cfg)
        graph_cfg = dict(cfg.get("graph", {}))
        tcfg = dict(cfg.get("transport", {}))
        cutoff = float(graph_cfg.get("cutoff_A", 3.8))
        graph = build_finite_graph(atoms, cutoff)

        eta = float(graph_cfg.get("interlayer_scale", 1.0))
        intralayer_cutoff = float(graph_cfg.get("intralayer_cutoff_A", cutoff))
        interlayer_cutoff = float(graph_cfg.get("interlayer_cutoff_A", cutoff))
        threshold = float(graph_cfg.get("interlayer_threshold_A", 1.0))
        layer_normal = graph_cfg.get("layer_normal", [0.0, 0.0, 1.0])

        result = source_drain_reachability(
            atoms,
            graph,
            direction=tcfg.get("direction", "X"),
            electrode_width_fraction=float(tcfg.get("electrode_width_fraction", 0.08)),
            layer_normal=layer_normal,
            eta=eta,
            intralayer_cutoff_A=intralayer_cutoff,
            interlayer_cutoff_A=interlayer_cutoff,
            interlayer_threshold_A=threshold,
        )
        summary = {
            "n_atoms": len(atoms),
            "n_source": int(len(result["source"])),
            "n_drain": int(len(result["drain"])),
            "n_reachable_source": int(len(result["reachable_source"])),
            "source_reachable_fraction": float(result["source_reachable_fraction"]),
            "all_nodes_reachable_fraction": float(result["all_nodes_reachable_fraction"]),
            "all_source_sites_reach_drain": bool(
                len(result["reachable_source"]) == len(result["source"])
            ),
        }
        outdir = base / cfg.get("output", {}).get("directory", "results")
        outdir.mkdir(parents=True, exist_ok=True)
        save_json(outdir / "reachability.json", summary)
        print(json.dumps(summary, indent=2))
        return 0

    raise ValueError(f"Unknown calculation mode: {mode}")


if __name__ == "__main__":
    raise SystemExit(main())
