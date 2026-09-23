# 02 - BN: anisotropic diffusion tensor

Generate direct-EHT parameters:

```bash
bash examples/paper/02_bn_tensor/run.sh
```

The calculation uses local B/N `2pz` orbitals and a 2x2x2 replication. The paper tensor benchmark intentionally keeps site energies equivalent and uses the EHT coupling hierarchy to isolate geometrical anisotropy. The raw B-N orbital-level difference is therefore diagnostic rather than a claimed quantitative trap energy.

The complete EHT-to-diffusion calculation is implemented in `reproduction/paper/atomic_cases.py`.
