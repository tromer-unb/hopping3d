# 05 - Benzene on graphene: hybrid molecule/slab interface

Parameterize both the single-adsorbate and four-adsorbate examples:

```bash
bash examples/paper/05_graphene_benzene/run.sh
```

The parser identifies the graphene host and molecular components automatically. A single EHT electronic calculation provides occupied and unoccupied frontier projections, so outputs contain separate hole/electron interface couplings. With several molecules, the code also evaluates molecule-molecule couplings.

The coupling network is an EHT screening result. Quantitative charge-transfer rates additionally require a validated molecule-substrate level alignment and carrier-specific reorganization energy. For pristine graphene, a reservoir/interfacial-transfer model is more appropriate than forcing the entire substrate into incoherent C-C hopping.
