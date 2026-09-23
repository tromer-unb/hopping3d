# Layered carbon / dimensional rescue

This example corresponds to the controlled four-layer carbon network used to separate topology from kinetics.

Build the structure if needed:

```bash
python3 build_layered_carbon.py
```

Estimate direct-EHT carbon couplings:

```bash
hopping3d-eht graphene stacked_graphene_4layers.cif -o carbon_eht.json
```

The article reproduction code is `reproduction/article_cases/run_atomic_cases.py`. The reported EHT interlayer rate scale is derived from `(J_inter/J_intra)^2`; the broader interlayer sweep is a mechanism map, not a fitted EHT calculation.
