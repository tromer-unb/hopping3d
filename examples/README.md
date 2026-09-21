# Examples

These examples are for users applying `hopping3d` to their own systems. The `reproduction/` directory is separate and exists only to reproduce the published article.

| Example | Status | What it demonstrates |
|---|---|---|
| `first_passage/` | executable | finite source-to-drain KMC |
| `diffusion_tensor/` | executable | periodic bulk diffusion tensor |
| `reachability/` | executable | pure graph connectivity, no KMC time scale |
| `graphene_li/` | template | host + impurity modeling choices |
| `molecular_benzene/` | physical-model template | molecule-as-site coarse graining and orientation dependence |
| `finite_nanocone/` | template | finite curved nanostructure with an arbitrary transport axis |

For the first three examples, install the package and run `hopping3d params.json` from the example directory.

The last three are deliberately documented as templates because the physically correct site/rate model must be chosen by the user; they are not intended to imply that a single default parameterization is quantitatively valid for all impurities, molecular solids, or nanostructures.
