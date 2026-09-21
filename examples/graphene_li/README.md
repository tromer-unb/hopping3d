# Graphene + Li impurity template

This directory documents three physically distinct reduced models for Li on or near graphene.

1. **Li as a perturbation only:** carbon atoms are hopping sites; Li shifts nearby carbon site energies.
2. **Li as a hopping site:** C and Li are both graph nodes, with separate `C-C`, `C-Li`, and `Li-Li` pair weights.
3. **Combined model:** Li has its own site energy and also perturbs nearby C sites.

The JSON template below illustrates model 2. Replace `graphene_li.cif` by your structure and treat the numerical parameters as placeholders unless they have been calibrated.

```bash
cp params_template.json params.json
hopping3d params.json
```

For quantitative work, site-energy shifts and C--Li coupling weights should be obtained from an electronic-structure model or constrained by experiment rather than assumed universal.
