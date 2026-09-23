# Changelog

## 2.2.0

- Reorganized the public repository around the five paper code cases.
- Added a versioned C192 layered-carbon input and per-case `run.sh` entry points.
- Consolidated numerical reproduction under `reproduction/paper/` with a single `reproduction/run_all.sh`.
- Rewrote the top-level workflow documentation around structure -> EHT -> transport.
- Kept manuscript TeX and publication figures out of the repository by design.
- Verified the clean-install test suite and quick numerical reproduction.

## 2.1.0 — 2026-09-22
- Synced the public transport core with the active EHT-aware engine, including explicit edge transfer integrals and Marcus rates.
- Added `hopping3d-eht bn` and `hopping3d-eht w2o6` commands.
- Added a detailed EHT-to-transport guide and paper-focused example READMEs.
- Simplified the example tree to the active article systems plus generic reachability/first-passage/diffusion examples.
- Made article reproduction scripts robust to environments that expose `python3` rather than `python`.

## 2.0.0 — 2026-09-22
- Added direct Extended-Huckel parametrization as `hopping3d-eht`.
- Added carrier-resolved molecular couplings: one EHT Hamiltonian, separate hole/electron frontier projections.
- Added BN and W-O article material parametrizers.
- Added molecular HOMO/LUMO coupling and graphene/benzene interface parametrization for one or multiple adsorbates.
- Added numerical reproduction workflows for layered carbon, BN, W2O6, the 100-benzene box, and graphene+benzene.
- Removed manuscript TeX and publication-figure assets from the repository; the public repository now focuses on code, structures, tests, documentation, and reproducible numerical datasets.

## 1.0.0 — 2026-09-21
- Public, package-based release of the 3D trajectory-resolved hopping framework.
- Finite first-passage KMC, periodic diffusion tensors, graph reachability diagnostics, and phenomenological chemistry controls.
