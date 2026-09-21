# hopping3d

[![CI](https://github.com/tromer-unb/hopping3d/actions/workflows/ci.yml/badge.svg)](https://github.com/tromer-unb/hopping3d/actions/workflows/ci.yml) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**Trajectory-resolved kinetic Monte Carlo for hopping transport in 2D, quasi-2D, and 3D networks.**

`hopping3d` reads an atomistic structure, builds a sparse hopping graph, and propagates continuous-time kinetic Monte Carlo trajectories. It is intended as a reusable extension of an earlier planar-network framework: the public interface is a structure file + one JSON parameter file + one command.

> **Important:** a normal user should not run the scripts inside `reproduction/` unless the goal is to reproduce the article. For a new material, use the installed `hopping3d` command.

## What do I actually download and run?

At present, `hopping3d` is distributed as a **Python package**, not as a single standalone binary downloaded from the Releases page.

After installation, `pip` creates the command-line program `hopping3d` in your Python environment. On Linux/macOS this is a console executable in the virtual environment; on Windows Python creates the corresponding `hopping3d.exe` launcher inside `.venv\Scripts`.

The normal workflow is therefore:

```text
structure.cif / structure.xyz / POSCAR / other ASE-readable file
                          +
                     params.json
                          |
                          v
                     hopping3d
                          |
       -----------------------------------------
       |                  |                    |
 first_passage      diffusion_tensor      reachability
       |                  |                    |
 finite device       periodic bulk         graph only
```

### Linux / macOS

```bash
git clone https://github.com/tromer-unb/hopping3d.git
cd hopping3d
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -e .

# check the installation
hopping3d --help
```

### Windows PowerShell

```powershell
git clone https://github.com/tromer-unb/hopping3d.git
cd hopping3d
py -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -e .

# check the installation
hopping3d --help
```

If `git` is unavailable, the repository can also be downloaded as a ZIP from GitHub, extracted, and installed with `pip install -e .` from the extracted directory.

## Which calculation should I use?

| Physical question | `calculation` value | Main output | Typical system |
|---|---|---|---|
| Bulk / periodic transport | `diffusion_tensor` | full diffusion tensor, principal diffusivities and axes | 3D crystal, layered solid, periodic molecular solid |
| Transport between source and drain | `first_passage` | `T_eff`, first-passage time, `nu0*tau_FP`, hops, path length, tortuosity | finite device, ribbon, flake, nanocone, finite 3D object |
| Pure connectivity, without KMC time scales | `reachability` | source-to-drain connected fraction | percolation / defect / dimensional-rescue analysis |
| Reproduce the published article | scripts under `reproduction/` | published CSV summaries and figures | article-specific only |

## Minimal examples

### 1. Finite source-to-drain transport

```bash
hopping3d examples/first_passage/params.json
```

Use this when the question is: **can a carrier cross this finite object, how long does it take, and how tortuous is the path?**

### 2. Periodic 3D diffusion tensor

```bash
hopping3d examples/diffusion_tensor/params.json
```

Use this for a periodic bulk material. The main result is

```text
Dxx Dxy Dxz
Dyx Dyy Dyz
Dzx Dzy Dzz
```

plus the principal diffusivities `D1 >= D2 >= D3` and their principal axes.

### 3. Pure graph reachability

Set

```json
{
  "calculation": "reachability"
}
```

and run

```bash
hopping3d params.json
```

No stochastic time scale is needed. This answers whether the source and drain are connected by the active graph.

## Generic input pattern

A minimal finite-device input is:

```json
{
  "calculation": "first_passage",
  "structure": {
    "file": "my_structure.cif",
    "repeat": [1, 1, 1]
  },
  "graph": {
    "cutoff_A": 3.5
  },
  "transport": {
    "direction": "X",
    "temperature_K": 300,
    "xi0_A": 3.0,
    "nu0_Hz": 1e13,
    "voltages_V": [0.02, 0.05, 0.10],
    "n_trajectories": 500,
    "max_time_s": 5e-10,
    "max_steps": 5000
  },
  "random_seed": 1234,
  "output": {
    "directory": "results"
  }
}
```

The preferred key is now `structure.file`. The older `structure.cif` key remains supported for backward compatibility. ASE is used internally, so CIF, XYZ, POSCAR and other ASE-readable structure formats can be used when they contain the geometry required by the selected calculation.

## Physical scope: what is a hopping site?

The present engine is **atomistic by default**:

```text
one atom  ->  one hopping site
```

This is appropriate for the atomic-network calculations reported with this code and for other models where localized states are intentionally attached to atomic positions.

For a molecular solid, that mapping may not be the physically appropriate coarse graining. For example, in a benzene box one often wants

```text
one benzene molecule  ->  one hopping site
```

with intermolecular rather than intramolecular charge-transfer rates. In that case the molecular center, site energy, transfer integral, reorganization energy, and molecular orientation can matter. The current public release does **not** automatically coarse-grain molecules or compute orientation-dependent molecular transfer integrals. See [`docs/use_cases.md`](docs/use_cases.md) for the recommended representation and extension path.

## Chemical heterogeneity and impurities

The current atomistic model can already assign:

- species-dependent site energies (`site_offsets_eV`);
- individual site-energy corrections (`site_offsets_by_index_eV` or `site_offsets_array_eV`);
- pair-dependent hopping weights (`pair_scales`);
- an optional local impurity-potential diagnostic.

This makes systems such as graphene + Li structurally possible within the atomistic model. The key physical decision is whether Li is itself treated as a hopping site, only as a perturbation to nearby carbon site energies, or both. The code does not decide that physics automatically; the model must be chosen explicitly.

## Finite nonplanar systems

`first_passage` does not require a planar structure. A finite nanocone, nanoparticle, curved sheet, cage, ribbon, or irregular 3D cluster can be treated as a finite graph. Source and drain are defined from the projection of atomic positions onto the selected transport direction.

For a nanocone, for example, one can choose the cone axis as `direction` and define source/drain near the two extrema of that projection. No periodicity is required for this calculation.

## What the code computes

- finite source-to-drain first-passage transport: `T_eff`, first-passage time, reduced time `nu0*tau_FP`, hop count, path length, and tortuosity;
- exact graph reachability to separate topological disconnection from slow kinetics;
- periodic 3D diffusion and mobility tensors, including principal diffusivities and axes;
- dimensional crossover through the transverse hopping scale `interlayer_scale`;
- phenomenological chemical heterogeneity through site-energy offsets and pair-dependent hopping weights;
- reproducible vacancy ensembles and trajectory-level path statistics.

## Reproducing the article

```bash
# regenerate Figures 2--4 directly from committed production summaries
python reproduction/scripts/plot_publication_figures.py

# rerun reduced-statistics simulations as a smoke test
bash reproduction/reproduce.sh quick

# rerun the production ensemble sizes used for the reported statistics
bash reproduction/reproduce.sh full
```

See [`docs/reproducibility.md`](docs/reproducibility.md) for the article-specific mapping between scripts, inputs, outputs, and reported observables.

## Repository map

- `src/hopping3d/` — reusable simulation engine and analysis utilities;
- `examples/` — minimal user-facing calculations;
- `reproduction/` — article-specific simulations, inputs, reference summaries, and figure generation;
- `manuscript/` — LaTeX source associated with the article;
- `docs/` — model, parameters, use cases, and reproducibility notes;
- `tests/` — numerical/unit smoke tests;
- `.github/workflows/` — multi-version CI and reproducibility checks.

## Physical interpretation

The default rate is Miller--Abrahams-like. `interlayer_scale`, site-energy offsets, and pair weights are phenomenological controls. They isolate dimensionality, connectivity, and chemistry without claiming that every parameter is a first-principles transfer integral.

For quantitative molecular charge transport, orientation-dependent electronic coupling and often Marcus-type rates are more appropriate. Such a molecular coarse-grained layer is a natural extension, but it should not be confused with the current atom-as-site model.

## Citation and license

See [`CITATION.cff`](CITATION.cff) for citation metadata. Released under the MIT License.
