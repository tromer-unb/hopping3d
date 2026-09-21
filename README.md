# hopping3d

[![CI](https://github.com/tromer-unb/hopping3d/actions/workflows/ci.yml/badge.svg)](https://github.com/tromer-unb/hopping3d/actions/workflows/ci.yml) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**Trajectory-resolved kinetic Monte Carlo for hopping transport in 2D, quasi-2D, and 3D atomic networks.**

`hopping3d` reads atomic structures directly from CIF files, builds sparse hopping graphs, and resolves transport through individual continuous-time kinetic Monte Carlo trajectories. The repository accompanies a published extension of an earlier planar-network framework and is organized so that a new user can run a minimal calculation in a few commands, while the full article statistics remain reproducible from explicit scripts and seeds.

## What the code computes

- finite source-to-drain first-passage transport: `T_eff`, first-passage time, reduced time `nu0*tau_FP`, hop count, path length, and tortuosity;
- exact graph reachability to separate topological disconnection from slow kinetics;
- periodic 3D diffusion and mobility tensors, including principal diffusivities and axes;
- dimensional crossover through the transverse hopping scale `eta_inter`;
- phenomenological chemical heterogeneity through site-energy offsets and pair-dependent hopping weights;
- reproducible vacancy ensembles and trajectory-level path statistics.

## Install

```bash
git clone https://github.com/tromer-unb/hopping3d.git
cd hopping3d
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

## 60-second examples

Finite first-passage calculation:

```bash
hopping3d examples/first_passage/params.json
```

Periodic diffusion tensor:

```bash
hopping3d examples/diffusion_tensor/params.json
```

Both examples use the same pattern as the article calculations: **CIF + JSON parameters -> atomic graph -> KMC -> CSV/JSON/figures**.

## Reproducing the article

```bash
# regenerate Figures 2--4 directly from the committed production summaries
python reproduction/scripts/plot_publication_figures.py

# rerun reduced-statistics simulations as a smoke test
bash reproduction/reproduce.sh quick

# rerun the production ensemble sizes used for the reported statistics
bash reproduction/reproduce.sh full
```

See [`docs/reproducibility.md`](docs/reproducibility.md) for the mapping between scripts, inputs, outputs, and reported observables.

## Repository map

- `src/hopping3d/` — reusable simulation engine and analysis utilities;
- `examples/` — minimal end-to-end calculations;
- `reproduction/` — article-specific simulations, input CIFs, reference summaries, and figure generation;
- `manuscript/` — LaTeX source associated with the article;
- `docs/` — model, parameters, and reproducibility notes;
- `tests/` — numerical/unit smoke tests;
- `.github/workflows/` — multi-version CI and reproducibility checks.

## Physical scope

The main rate is a Miller--Abrahams-type hopping rate. `eta_inter`, species/site energy offsets, and pair weights are intentionally exposed as phenomenological controls. They isolate dimensionality, connectivity, and chemistry without implying that every parameter is a first-principles transfer integral. The code also retains selected legacy diagnostics for controlled ablation, but they are not enabled by default.

## Citation and license

See [`CITATION.cff`](CITATION.cff) for citation metadata. Released under the MIT License.
