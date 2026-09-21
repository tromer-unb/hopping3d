# Molecular benzene template

For a benzene molecular solid, the physically natural coarse graining is usually

```text
one benzene molecule = one hopping site
```

rather than one atom = one site.

## Orientation matters

Two molecules at the same center-to-center distance should not generally be assigned the same transfer rate if their relative orientations differ. A face-to-face pi-stacked pair can have a very different electronic coupling from an edge-on or nearly perpendicular pair.

A molecular model should therefore use either:

- externally supplied pair transfer integrals `J_ij`; or
- an empirical orientation-dependent coupling based on distance, plane-normal alignment, and slip.

For quantitative molecular charge transport, Marcus-type rates are often more appropriate than the default atomic Miller--Abrahams rate because they can include `J_ij`, reorganization energy `lambda`, and the driving-force/site-energy difference.

## Current status

`hopping3d` v1.0 does **not** automatically identify molecules, replace them by centers, or compute orientation-dependent molecular transfer integrals. The JSON template in this directory is therefore a model-design template, not a ready-to-run quantitative benzene calculation.

A proper preprocessing stage should generate one transport site per molecule and provide the molecular site energies and pair couplings before KMC.
