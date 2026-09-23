from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from .core import (load_json, load_structure, apply_random_vacancies, run_from_config,
                   simulate_periodic_diffusion, save_tensor_csv, save_json,
                   prepare_from_config)
from .tensor import principal_decomposition
from .graph import source_drain_reachability, source_drain_reachability_generic


def _resolve_atoms(config_path: Path, cfg: dict):
    base=config_path.resolve().parent
    scfg=cfg.get("structure", {})
    structure_name=scfg.get("file", scfg.get("cif"))
    if structure_name is None:
        raise KeyError("structure.file (or legacy structure.cif) is required")
    structure_path=Path(structure_name)
    if not structure_path.is_absolute(): structure_path=base/structure_path
    atoms=load_structure(structure_path, scfg.get("repeat", [1,1,1]))
    atoms=apply_random_vacancies(atoms, scfg.get("vacancy_fraction",0.0),
                                 scfg.get("vacancy_seed",0), scfg.get("vacancy_symbols",None))
    return base, atoms


def main(argv=None):
    p=argparse.ArgumentParser(prog="hopping3d", description="Trajectory-resolved hopping KMC from CIF structures")
    p.add_argument("config", help="JSON configuration file")
    args=p.parse_args(argv)
    cfg=load_json(args.config); mode=str(cfg.get("calculation","first_passage")).lower()
    if mode in {"first_passage","device"}:
        df=run_from_config(args.config); print(df.to_string(index=False)); return 0
    if mode in {"reachability","connectivity"}:
        config_path=Path(args.config); full_cfg,base,atoms,graph=prepare_from_config(config_path)
        gcfg=dict(full_cfg.get("graph",{})); tcfg=dict(full_cfg.get("transport",{}))
        direction=tcfg.get("direction","X")
        width_A=tcfg.get("electrode_width_A",None)
        width_fraction=float(tcfg.get("electrode_width_fraction",0.05))
        layered=any(k in gcfg for k in ("interlayer_scale","intralayer_cutoff_A","interlayer_cutoff_A"))
        if layered:
            result=source_drain_reachability(
                atoms, graph, direction=direction, electrode_width_fraction=width_fraction,
                layer_normal=gcfg.get("layer_normal",[0,0,1]),
                eta=float(gcfg.get("interlayer_scale",1.0)),
                intralayer_cutoff_A=float(gcfg.get("intralayer_cutoff_A",gcfg.get("cutoff_A",3.0))),
                interlayer_cutoff_A=float(gcfg.get("interlayer_cutoff_A",gcfg.get("cutoff_A",3.0))),
                interlayer_threshold_A=float(gcfg.get("interlayer_threshold_A",1.0)))
        else:
            result=source_drain_reachability_generic(
                atoms, graph, direction=direction, electrode_width_A=width_A,
                electrode_width_fraction=width_fraction)
        serial={k:(sorted(v) if isinstance(v,set) else np.asarray(v).tolist() if hasattr(v,"shape") else v)
                for k,v in result.items()}
        outdir=base/full_cfg.get("output",{}).get("directory","results"); outdir.mkdir(parents=True,exist_ok=True)
        save_json(outdir/"reachability.json",serial)
        print("source_reachable_fraction =", result["source_reachable_fraction"])
        print("all_nodes_reachable_fraction =", result["all_nodes_reachable_fraction"])
        return 0
    if mode in {"diffusion","diffusion_tensor","periodic_diffusion"}:
        config_path=Path(args.config); base,atoms=_resolve_atoms(config_path,cfg)
        graph_cfg=dict(cfg.get("graph",{})); tcfg=dict(cfg.get("transport",{})); tcfg.update(graph_cfg)
        chem=cfg.get("chemistry",{}); dcfg=cfg.get("diffusion",{})
        outdir=base/cfg.get("output",{}).get("directory","results"); outdir.mkdir(parents=True,exist_ok=True)
        cutoff=float(graph_cfg.get("cutoff_A",3.8))
        result=simulate_periodic_diffusion(atoms, cutoff, tcfg, chem,
            n_walkers=int(dcfg.get("n_walkers",500)),
            observation_time_s=float(dcfg.get("observation_time_s",8e-11)),
            seed=int(cfg.get("random_seed",1234)))
        D=np.asarray(result["D_m2_s"],float); vals,vecs=principal_decomposition(D)
        save_tensor_csv(outdir/"diffusion_tensor.csv",D,"D")
        result["principal_diffusivities_m2_s"]=vals.tolist(); result["principal_axes"]=vecs.tolist()
        save_json(outdir/"diffusion_summary.json",result)
        print("D tensor (m^2/s):\n", D)
        print("principal diffusivities:", vals)
        return 0
    raise ValueError(f"Unknown calculation mode: {mode}")

if __name__ == "__main__":
    raise SystemExit(main())
