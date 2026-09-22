# Direct Extended-Huckel parametrization

`hopping3d-eht` is the low-cost electronic parametrization layer used by the current workflow. It reads structures through ASE and returns geometry-dependent electronic couplings; it is **not** a replacement for DFT/Wannier when quantitative level alignment, polarization, charge transfer, reorganization energies, or absolute mobilities are required.

## One electronic calculation, two carrier channels

For a fixed geometry the Extended-Huckel Hamiltonian is built once. Occupied and unoccupied frontier subspaces are then projected separately:

- **holes** use occupied/HOMO-like states and produce `J_hole`;
- **electrons** use unoccupied/LUMO-like states and produce `J_electron`.

Thus EHT supplies both channels from the same electronic problem. Transport is then run as two independent single-carrier calculations because electrons and holes can have different couplings, site energies, level offsets, and Marcus reorganization energies.

## Active orbital choices

| system | EHT transport subspace | main output |
|---|---|---|
| graphene / carbon stack | C 2pz | C-C, intra/inter-layer couplings |
| BN | B/N 2pz | in-plane/interlayer couplings |
| W-O benchmark | O 2p + W 5d | W-O and O-O coupling hierarchy |
| molecular box | molecular HOMO/LUMO subspaces | hole/electron J_MM |
| graphene + molecule | graphene C 2pz + molecular frontier subspace | hole/electron J_GM and J_MM |

Do not tune EHT until a desired transport result is obtained. Use standard parameters for mechanism screening. If the result lies near a crossover, refine the sensitive electronic quantity with DFT, constrained DFT, Wannier functions, higher-level quantum chemistry, or experiment.

## CLI examples

```bash
hopping3d-eht doctor
hopping3d-eht graphene examples/article_systems/carbon_stack/stacked_graphene_4layers.cif -o carbon_eht.json
hopping3d-eht interface examples/article_systems/graphene_benzene/graphene_benzene.cif -o interface.json --report interface.md
hopping3d-eht coverage examples/article_systems/graphene_benzene/graphene_4benzene.cif -o coverage.json --report coverage.md
```

The parameterization output is an electronic input to the transport model. Reorganization energies and quantitatively reliable molecule-substrate level alignments must be supplied separately when the chosen rate law needs them.
