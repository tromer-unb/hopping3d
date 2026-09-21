# Model summary

The code represents a crystal as an atomic graph built from a CIF structure. Directed Miller--Abrahams rates are evaluated on graph edges,

`Gamma_ij = nu0 * w_ij * exp(-2 r_ij / xi) * min[1, exp(-(E_j-E_i)/(k_B T))]`.

The framework deliberately separates (i) graph connectivity, (ii) kinetic accessibility, and (iii) local energetic/chemical contrast. Finite devices provide source-to-drain first-passage observables; periodic calculations provide a 3x3 diffusion tensor. The phenomenological parameters `interlayer_scale`, `site_offsets_eV`, and `pair_scales` are model controls and are not presented as first-principles transfer integrals.
