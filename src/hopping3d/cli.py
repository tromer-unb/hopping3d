from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from .core import (load_json, load_structure, apply_random_vacancies, run_from_config,
                   simulate_periodic_diffusion, save_tensor_csv, save_json)
from .tensor import principal_decomposition


def _resolve_atoms(config_path: Path, cfg: dict):
    base=config_path.resolve().parent
    scfg=cfg.get("structure", {})
    cif=Path(scfg["cif"])
    if not cif.is_absolute(): cif=base/cif
    atoms=load_structure(cif, scfg.get("repeat", [1,1,1]))
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
