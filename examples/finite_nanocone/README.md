# Finite nanocone template

A nanocone is a finite, curved atomic network and should normally be treated with `first_passage` or `reachability`, not with `diffusion_tensor`.

Choose a transport direction corresponding to the physical experiment. For apex-to-base transport, use the cone axis. Source and drain are selected from the extrema of the atomic positions projected onto that direction.

Example workflow:

```text
nanocone structure
      |
      v
finite atomic graph
      |
      +--> reachability: does an apex-to-base path exist?
      |
      +--> first_passage: how long does traversal take and which path is used?
```

Use `params_template.json` after replacing `nanocone.cif` by a valid ASE-readable structure file and setting the transport direction to the actual cone axis.

No periodic boundary conditions are required for this finite-device calculation.
