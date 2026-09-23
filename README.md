# Hopping3D

Hopping3D is a compact research code for **localized-state transport on atomic, molecular, and hybrid networks**. The current release couples a direct Extended-Huckel (EHT) parametrization layer to graph analysis, continuous-time kinetic Monte Carlo (CTMC/KMC), first-passage observables, and diffusion tensors.

The design is intentionally mechanism-first: use inexpensive EHT couplings to identify robust regimes (connectivity, anisotropy/isotropy, chemical-path selection, molecular-orientation effects, molecule-surface coupling), then replace sensitive parameters with DFT/Wannier/experiment when quantitative accuracy is required.

## Install

```bash
git clone https://github.com/tromer-unb/hopping3d.git
cd hopping3d
python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -e .
```

Check the EHT runtime:

```bash
hopping3d-eht doctor
```

## Two-stage workflow

```text
structure (CIF/XYZ/POSCAR/...)
        |
        v
  hopping3d-eht
        |
        +--> hole couplings / screening energies
        +--> electron couplings / screening energies
        |
        v
 transport input
        |
        v
    hopping3d
        |
        +--> reachability
        +--> first passage / trapping / tortuosity
        +--> diffusion tensor / principal axes
```

**One EHT electronic problem produces both occupied (hole) and unoccupied (electron) projections. Electron and hole transport are then separate single-carrier simulations.** See `docs/carriers.md`.

## EHT examples

Graphene/carbon:

```bash
hopping3d-eht graphene examples/article_systems/carbon_stack/stacked_graphene_4layers.cif -o carbon_eht.json
```

BN:

```bash
hopping3d-eht bn examples/article_systems/bn/BN_bulk.cif --repeat 2 2 2 -o bn_eht.json
```

W2O6:

```bash
hopping3d-eht w2o6 examples/article_systems/w2o6/W2O6.cif --cutoff 3.4 -o w2o6_eht.json
```

One molecule on graphene:

```bash
hopping3d-eht interface examples/article_systems/graphene_benzene/graphene_benzene.cif \
  -o interface.json --report interface.md
```

Several molecules on graphene:

```bash
hopping3d-eht coverage examples/article_systems/graphene_benzene/graphene_4benzene.cif \
  -o coverage.json --report coverage.md
```

## Reproduce article numerical datasets

No manuscript TeX or publication figures are stored in this repository. The numerical source datasets are reproduced with:

```bash
bash reproduction/article_cases/run_all.sh
```

Production statistics:

```bash
FULL=1 python3 reproduction/article_cases/run_atomic_cases.py
FULL=1 python3 reproduction/article_cases/run_benzene_box.py
python3 reproduction/article_cases/run_interface.py
```

See `docs/article_reproduction.md` for the mapping from scripts to article cases.

## What is EHT used for?

The same EHT engine is projected onto different physically motivated subspaces:

| case | subspace |
|---|---|
| carbon / graphene | C 2pz |
| BN | B/N 2pz |
| W-O | O 2p + W 5d |
| molecular solid | molecular HOMO/LUMO manifolds |
| molecule on graphene | graphene 2pz + molecular frontier manifolds |

EHT is a **screening parametrization**, not a universal DFT replacement. Relative coupling hierarchies are generally the intended use. Absolute quasiparticle level alignment, self-consistent charge transfer, polarization/image-charge corrections, defect energetics, and nuclear reorganization energies should be refined when they control the conclusion.

## Documentation

- `docs/workflow.md` — end-to-end workflow.
- `docs/eht_parameterization.md` — how EHT parameters are obtained and interpreted.
- `docs/eht_to_transport.md` — detailed structure -> EHT -> Hopping3D workflow.
- `docs/carriers.md` — electrons vs holes.
- `docs/parameters.md` — transport JSON parameters.
- `docs/model.md` — graph/KMC model.
- `docs/article_reproduction.md` — article code/data reproduction.
- `examples/README.md` — example systems.

## Scope

This code is intended for mechanism-resolved localized transport and inexpensive screening. It does not claim that pristine delocalized conductors are universally describable by incoherent hopping, and it does not infer high-accuracy electronic energetics from EHT.

## License

MIT. See `LICENSE`.
