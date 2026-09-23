# BN tensor benchmark

The article uses a 2x2x2 replication of `BN_bulk.cif` to demonstrate anisotropic diffusion from direct EHT couplings.

Generate the EHT parameter file:

```bash
hopping3d-eht bn BN_bulk.cif --repeat 2 2 2 -o bn_eht.json
```

The active electronic subspace is local B/N `2pz`. In the tensor benchmark the site-energy contrast is intentionally set to zero so that the result isolates the directional coupling hierarchy. The raw B--N EHT orbital offset is diagnostic only.

Reproduce the transport calculation from the repository root with:

```bash
python3 reproduction/article_cases/run_atomic_cases.py
```
