# From a structure to EHT parameters and transport

This page is the practical bridge between the two executables shipped with the repository:

```text
structure -> hopping3d-eht -> electronic parameters -> hopping3d -> observables
```

The EHT layer is intentionally a low-cost screening model. It is useful for relative coupling hierarchies and mechanism studies; it is not a universal replacement for DFT/Wannier or experiment.

## 1. Check the installation

```bash
hopping3d-eht doctor
```

A normal installation requires NumPy, SciPy and ASE. No external electronic-structure executable is required for the direct-EHT workflows documented here.

## 2. Atomic carbon / graphene

```bash
hopping3d-eht graphene structure.cif -o carbon_eht.json
```

The output contains a suggested graph cutoff, a reference transfer integral `J_ref_eV`, and an edge map `edge_transfer_integrals_eV`. For explicit EHT couplings, Hopping3D uses rate weights proportional to `|J_ij/J_ref|^2`; do not multiply by a second empirical distance exponential unless that is intentionally part of a different model.

## 3. Layered BN

```bash
hopping3d-eht bn examples/article_systems/bn/BN_bulk.cif \
  --repeat 2 2 2 \
  -o bn_eht.json
```

The EHT subspace is B/N local `2pz`. The output reports the in-plane reference coupling, interlayer couplings, and the raw B--N orbital-level contrast. In the article tensor benchmark the B/N site energies are kept equal so that the calculation isolates the geometrical coupling hierarchy; the raw EHT level contrast is diagnostic, not asserted to be a quantitative carrier trap energy.

## 4. W2O6

```bash
hopping3d-eht w2o6 examples/article_systems/w2o6/W2O6.cif \
  --cutoff 3.4 \
  -o w2o6_eht.json
```

The active reduction uses O `2p` and W `5d` channels. The output includes `J_ref_WO_eV`, edge-specific couplings, the O--O rate-weight distribution, and the raw W(5d)--O(2p) level contrast. At the article cutoff there are no W--W edges. The raw level contrast is used as a screening scale only.

## 5. Molecular systems: electrons and holes

A single EHT Hamiltonian is built for a molecular pair. The occupied frontier manifold is projected to obtain hole coupling `J_hole`; the unoccupied frontier manifold is projected to obtain electron coupling `J_electron`.

```text
one EHT Hamiltonian
   +-> occupied/HOMO-like projection -> hole parameters
   +-> empty/LUMO-like projection    -> electron parameters
```

The transport calculations are then separate single-carrier simulations. They may use different `J_ij`, driving forces and reorganization energies.

## 6. Periodic molecular box

The article example is generated and run with:

```bash
python3 reproduction/article_cases/run_benzene_box.py
```

or, with production statistics:

```bash
FULL=1 python3 reproduction/article_cases/run_benzene_box.py
```

The script constructs 100 benzene molecules, identifies 300 nearest-neighbor molecular pairs, evaluates EHT hole/electron couplings for every pair, and propagates a Marcus diffusion tensor. The Marcus reorganization energy is an external parameter (`lambda = 0.20 eV` in the controlled benchmark), not an EHT prediction.

## 7. Benzene on graphene / slab

For one adsorbate:

```bash
hopping3d-eht interface examples/article_systems/graphene_benzene/graphene_benzene.cif \
  -o interface.json --report interface.md
```

For several adsorbates:

```bash
hopping3d-eht coverage examples/article_systems/graphene_benzene/graphene_4benzene.cif \
  -o coverage.json --report coverage.md
```

The parser identifies the graphene host and molecular components and returns both graphene--molecule (`J_GM`) and molecule--molecule (`J_MM`) couplings. A quantitative transfer rate also needs a validated level alignment and reorganization energy.

## 8. Feeding explicit EHT couplings to Hopping3D

A periodic atomic calculation can use a transport block of the form:

```json
{
  "calculation": "diffusion_tensor",
  "structure": {"file": "structure.cif"},
  "graph": {"cutoff_A": 3.8},
  "transport": {
    "temperature_K": 300,
    "rate_model": "miller_abrahams",
    "nu0_Hz": 1.0e13,
    "xi0_A": 1.0e9,
    "J_ref_eV": 1.0,
    "edge_transfer_integrals_eV": {"0-1": 0.1},
    "strict_edge_transfer_integrals": true
  },
  "chemistry": {"site_offsets_array_eV": [0.0]},
  "diffusion": {"n_walkers": 1000, "observation_time_s": 1.0e-10},
  "output": {"directory": "results"}
}
```

The example above is schematic: array lengths and edge indices must match the actual structure. The article reproduction scripts show complete programmatic constructions for BN and W2O6.

For molecular Marcus transport use `rate_model: "marcus"`, explicit `edge_transfer_integrals_eV`, and a physically justified `reorganization_energy_eV`. Carrier-specific values should be used for electron and hole runs when available.

## 9. What to refine with DFT

Refine a parameter when the scientific conclusion is sensitive to it. Typical high-priority quantities are:

- molecule--substrate level alignment;
- defect/trap site energies;
- nuclear reorganization energies;
- screened charge-transfer energetics;
- couplings near a mechanism crossover.

The transport engine does not change when EHT values are replaced by DFT/Wannier/experimental values. This is the intended hierarchy: use EHT to locate the mechanism cheaply, then spend higher-level electronic-structure effort only where it can change the conclusion.
