# Recommended workflow

## 1. Start from a structure

ASE-supported formats are accepted. Verify cell vectors, periodicity, bond lengths, vacuum, and molecular connectivity before interpreting transport.

## 2. Parameterize electronically

Use `hopping3d-eht` for low-cost screening. For atomic systems the code projects selected atomic orbitals. For molecular systems it identifies/uses molecular frontier subspaces. For molecule/slab systems it computes substrate-molecule and molecule-molecule coupling blocks.

## 3. Inspect confidence and missing physics

EHT couplings are useful for relative hierarchy and mechanism maps. Absolute level alignment, self-consistent charge transfer, dielectric/image-charge effects, and Marcus reorganization energies may require external values.

## 4. Build the transport input

Use explicit edge couplings whenever available. Hopping3D maps them to rate weights without applying a second distance-decay factor. Choose Miller-Abrahams for generic localized atomic hopping or Marcus for molecular electron transfer when `lambda` and driving forces are available.

## 5. Choose the observable

- `reachability`: topology only;
- `first_passage`: source-to-drain kinetics, trapping, tortuosity;
- `diffusion_tensor`: periodic bulk/molecular diffusion and anisotropy.

## 6. Converge statistics

Increase walkers/trajectories, observation time, structural realizations, and cell size until the observable of interest is stable. A finite observation window must not be confused with true topological disconnection.
