# Minimal first-passage example

From the repository root:

```bash
hopping3d examples/first_passage/params.json
```

The example reads `structure.cif`, removes a reproducible 10% fraction of sites, builds a finite hopping graph, and runs source-to-drain KMC at three biases. Results are written to `examples/first_passage/results/`.
