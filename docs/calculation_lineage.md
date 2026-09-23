# Calculation lineage

This repository is the curated public endpoint of a longer development sequence. The public tree is intentionally organized by **reproducible scientific role** rather than by every intermediate working-directory name.

## Public organization

| Public location | Role |
|---|---|
| `src/hopping3d/` | graph construction, transport, first-passage, diffusion-tensor machinery |
| `src/hopping3d_eht/` | low-cost Extended-Huckel parametrization used by the article workflows |
| `examples/paper/` | small, readable entry points for the five article systems |
| `reproduction/paper/` | executable scripts for regenerating the numerical paper datasets |
| `paper/data/` | frozen processed datasets and structures consumed by the manuscript figures |
| `paper/figures/` | generated publication figures reconstructed from the frozen article data (not versioned) |
| `paper/` | revised manuscript and Supplemental Material sources |

## Development sequence

The work evolved through several classes of calculations before reaching the current article-oriented layout:

1. **Legacy voltage-sweep calculations** established the original hopping-network workflow and supplied historical structure/input snapshots.
2. **3D mechanism studies** tested dimensional crossover, anisotropic diffusion, chemistry-aware paths, impurity effects, and finite-size/convergence behavior.
3. **Correction/validation stages** revised geometry handling, tensor extraction, electronic-parameter interfaces, transport-regime diagnostics, and production statistics.
4. **Extension studies** added molecular orientation, Marcus-rate tests, finite-size checks, lithium robustness, and related ensemble calculations.
5. **EHT-focused studies** introduced direct electronic coupling estimates for BN, W2O6, molecular benzene systems, and graphene/benzene interfaces.
6. **Manuscript revision 2** froze the EHT-anchored article datasets and figures now stored under `paper/`.

## What is intentionally not copied verbatim

Intermediate caches, virtual environments, LaTeX build products, repeated generated figures, temporary logs, and multi-gigabyte development run directories are not mirrored wholesale into Git. Their scientifically relevant outputs are represented by the curated examples, reproduction scripts, frozen processed datasets, tests, and repository history.

This separation keeps the repository cloneable while preserving the inputs and outputs needed to inspect the claims made in the revised manuscript.
