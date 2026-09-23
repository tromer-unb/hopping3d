# Choosing the physical representation

The central modeling choice is the definition of a **transport site**. Hopping3D does not assume that every system should use the same coarse graining.

| system | transport site | electronic parametrization | transport observable |
|---|---|---|---|
| layered carbon | atom / C 2pz | direct EHT C--C couplings | reachability + first passage |
| layered BN | atom / B,N 2pz | direct EHT B--N couplings | diffusion tensor |
| W2O6 | atom-site reduction | direct EHT O 2p + W 5d blocks | first passage + chemical visitation |
| molecular box | one molecule | HOMO/LUMO projected EHT couplings | Marcus diffusion tensor |
| molecule on graphene | graphene 2pz + molecular frontier site | interface EHT projection | interface/coverage kinetics |

## Atomic systems

For the atomic article cases, nodes are atomic localized-state proxies and the graph is constructed from a geometric cutoff. EHT supplies edge-specific transfer integrals. Site-energy differences from raw EHT orbitals are treated as screening scales unless independently validated.

## Molecular systems

For a molecular solid, **one molecule is one transport site**. A box containing 100 benzene molecules therefore becomes a 100-site graph. For every neighboring molecular pair the EHT dimer Hamiltonian is projected onto occupied and unoccupied frontier manifolds, yielding separate hole and electron couplings. Marcus transport additionally requires a reorganization energy `lambda`, which is not supplied by ordinary EHT and must be provided externally.

## Molecules on a slab

For graphene plus adsorbates, the structure is decomposed into the graphene host and molecular components. The electronic layer returns both graphene--molecule (`J_GM`) and molecule--molecule (`J_MM`) couplings. For a delocalized conductor such as pristine graphene, a reservoir/interface treatment is usually more physical than forcing the substrate into incoherent C--C hopping. Quantitative charge-transfer rates require a reliable level alignment and reorganization energy.

## Electrons and holes

The EHT Hamiltonian is built once. Occupied frontier projections generate hole parameters and unoccupied frontier projections generate electron parameters. The subsequent transport simulations are **separate single-carrier ensembles**. See `docs/carriers.md`.

## What EHT is for

Use EHT to identify coupling hierarchies, weak/strong directions, molecular-orientation effects, and chemical-path selection. If a conclusion depends sensitively on a site energy, level alignment, screening correction, or reorganization energy, replace that quantity with DFT/Wannier/constrained-DFT/experiment rather than tuning EHT to obtain a desired result.
