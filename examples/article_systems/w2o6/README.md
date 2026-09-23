# W2O6 chemical-path benchmark

This example demonstrates chemically heterogeneous transport with an O `2p` / W `5d` EHT reduction.

Generate the EHT parameters:

```bash
hopping3d-eht w2o6 W2O6.cif --cutoff 3.4 -o w2o6_eht.json
```

At the 3.4 A article cutoff the graph contains W--O and O--O edges but no W--W edges. The raw W(5d)--O(2p) level difference is a **screening scale**, not a DFT-quality trap energy.

The transport comparison between zero site contrast and the raw EHT screening contrast is implemented in:

```bash
python3 reproduction/article_cases/run_atomic_cases.py
```
