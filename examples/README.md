# Examples

The `paper/` directory contains the five physical cases discussed in the paper. Every case has its own README and an executable `run.sh`.

```text
paper/
  01_carbon_stack/       dimensional rescue / weak interlayer coupling
  02_bn_tensor/          anisotropic diffusion tensor
  03_w2o6_chemistry/     W-O / O-O chemical-path selection
  04_benzene_box/        periodic box of 100 benzene molecules
  05_graphene_benzene/   molecule-on-2D-slab interface
```

The three small generic examples (`reachability`, `first_passage`, and `diffusion_tensor`) are kept as minimal tests of the transport engine independent of EHT.

Generated JSON/CSV output files are not versioned. Re-run the scripts to create them locally.
