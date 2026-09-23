# Hopping3D

Hopping3D is a research code for **localized-state transport on atomic, molecular, and hybrid networks**. The repository contains the code paths used in the accompanying paper: direct Extended-Huckel (EHT) parametrization, graph construction, CTMC/KMC transport, first-passage analysis, and diffusion tensors.

The repository contains the research code, example structures, tests, numerical reproduction scripts, and the revised manuscript package. The article source, frozen processed datasets, and reproducible figure-building workflow live under `paper/`; executable numerical workflows remain separated under `examples/paper/` and `reproduction/paper/`.

## What the code does

```text
CIF / XYZ / POSCAR / ASE structure
              |
              v
        hopping3d-eht
              |
      electronic parameters
      J_ij, screening levels
              |
              v
          hopping3d
              |
   reachability / first passage /
   diffusion tensor / anisotropy
```

EHT is used as a **low-cost screening parametrization**, not as a universal replacement for DFT. The intended use is to establish robust transport physics first (connectivity, anisotropy, chemical-path selection, orientation effects), then refine sensitive quantities with DFT/Wannier/experiment if needed.

## Install

```bash
git clone https://github.com/tromer-unb/hopping3d.git
cd hopping3d
python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -e .
```

Check the installation:

```bash
hopping3d-eht doctor
pytest
```

## Paper code examples

All paper-related examples are under `examples/paper/`:

| directory | physical problem | EHT subspace |
|---|---|---|
| `01_carbon_stack` | dimensional rescue in layered carbon | C 2pz |
| `02_bn_tensor` | anisotropic 3D diffusion tensor | B/N 2pz |
| `03_w2o6_chemistry` | chemical-path selection | O 2p + W 5d |
| `04_benzene_box` | 100-molecule periodic box | molecular HOMO/LUMO |
| `05_graphene_benzene` | molecule/2D-slab interface | graphene 2pz + molecular frontier states |

Run any case with its local `run.sh`, for example:

```bash
bash examples/paper/02_bn_tensor/run.sh
bash examples/paper/05_graphene_benzene/run.sh
```

## Manuscript package

The revised EHT-anchored article is versioned under `paper/`:

- `paper/main.tex` and `paper/SI.tex` — article and Supplemental Material sources;
- `paper/data/` — frozen processed datasets and example structures used by the figures;
- `paper/figures/` — generated PNG/SVG publication figures (ignored by Git);
- `paper/scripts/` — deterministic figure-building scripts.

To rebuild figures and compile the manuscript:

```bash
cd paper
python3 scripts/build_figures.py
bash build.sh
```

Generated figures, LaTeX build files, and compiled PDFs are intentionally not tracked.

## Reproduce all numerical paper datasets

Quick smoke-size reproduction:

```bash
bash reproduction/run_all.sh
```

Production statistics for the heavier ensemble calculations:

```bash
FULL=1 python3 reproduction/paper/atomic_cases.py
FULL=1 python3 reproduction/paper/benzene_box.py
python3 reproduction/paper/graphene_benzene.py
```

Generated files are written to `reproduction/paper/results/` and are ignored by Git.

## EHT first, transport second

For an atomic structure the normal workflow is:

```bash
hopping3d-eht bn examples/paper/02_bn_tensor/BN_bulk.cif \
  --repeat 2 2 2 -o bn_eht.json
```

The EHT output contains the edge-resolved couplings and the reference coupling used by the transport model. The reproduction scripts show the exact mapping from these values to Hopping3D rate parameters.

When explicit EHT couplings are supplied, Hopping3D uses weights proportional to

```text
|J_ij / J_ref|^2
```

and does **not** add a second empirical distance decay unless the user explicitly requests a different model.

## Electrons and holes

The EHT electronic problem is solved once for a given geometry. Molecular occupied and unoccupied frontier subspaces are then projected from the same Hamiltonian:

```text
one EHT Hamiltonian
  +-- HOMO-like / occupied  -> hole J_ij
  +-- LUMO-like / unoccupied -> electron J_ij
```

Electron and hole transport are **separate single-carrier simulations** because they may have different couplings, driving energies, and reorganization energies. See `docs/carriers.md`.

## Documentation

- `docs/getting_started.md` — first run from clone to result.
- `docs/eht_parameterization.md` — EHT basis, projections, outputs, and limitations.
- `docs/eht_to_transport.md` — detailed mapping from EHT output to hopping rates.
- `docs/carriers.md` — electron vs hole channels.
- `docs/model.md` — graph, CTMC/KMC, first passage, and diffusion tensor.
- `docs/parameters.md` — Hopping3D JSON fields.
- `docs/article_reproduction.md` — mapping from paper cases to scripts.
- `docs/calculation_lineage.md` — how the legacy, correction, extension, EHT, and revised-manuscript stages map into the curated public tree.

## Scope and limitations

EHT is most useful here for **relative coupling hierarchies and mechanism screening**. Quantitative quasiparticle level alignment, self-consistent charge transfer, dielectric/image-charge corrections, charged-defect energetics, and nuclear reorganization energies should be supplied by a higher-fidelity method when they control the conclusion.

For pristine delocalized conductors such as graphene, the repository does not claim that the substrate itself is universally described by incoherent C-C hopping. The graphene/benzene example is primarily an interface-coupling parametrization; a reservoir/interfacial-transfer description is more appropriate for quantitative charge transfer.

## License

MIT. See `LICENSE`.
