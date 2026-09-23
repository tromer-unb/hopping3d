# Benzene on graphene: hybrid atom--molecule interface

Two structures are provided: one benzene adsorbate and four benzene adsorbates on graphene.

Single adsorbate:

```bash
hopping3d-eht interface graphene_benzene.cif -o interface.json --report interface.md
```

Multiple adsorbates / coverage:

```bash
hopping3d-eht coverage graphene_4benzene.cif -o coverage.json --report coverage.md
```

A single EHT Hamiltonian provides occupied and unoccupied frontier projections, so the output contains separate hole/electron interface couplings. Transport of electrons and holes is subsequently performed as separate single-carrier calculations. Quantitative transfer rates also require a validated molecule--substrate level alignment and a reorganization energy.
