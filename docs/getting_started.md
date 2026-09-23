# Getting started

## 1. Create an environment

```bash
git clone https://github.com/tromer-unb/hopping3d.git
cd hopping3d
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```

Verify both command-line programs:

```bash
hopping3d-eht doctor
hopping3d --help
```

## 2. Choose the physical representation

Before running transport, decide what a hopping site represents.

- Atomic localized state: one site per atom (carbon stack, BN, W-O benchmark).
- Molecular solid: one site per molecule; EHT is projected onto frontier molecular subspaces.
- Molecule on a slab: substrate states and molecular frontier states are projected separately.

This choice is more important than any numerical KMC parameter.

## 3. Parametrize with EHT

Example for BN:

```bash
hopping3d-eht bn examples/paper/02_bn_tensor/BN_bulk.cif \
  --repeat 2 2 2 -o bn_eht.json
```

Inspect the JSON. It contains the reference coupling, edge-specific transfer integrals, and diagnostic electronic energy scales.

## 4. Run transport

The transport executable consumes a JSON configuration:

```bash
hopping3d path/to/transport.json
```

Complete programmatic examples are in `reproduction/paper/`. They deliberately show the EHT-to-rate mapping in Python rather than hiding it behind a black-box command.

## 5. Reproduce the paper cases

```bash
bash reproduction/run_all.sh
```

Use `FULL=1` for the larger statistical ensembles described in the paper.

## 6. Interpret results at the correct accuracy level

Robust uses include connectivity, principal transport directions, anisotropy/isotropy, and strong chemical-path selection. Absolute mobility, trap depths, level alignment, and reorganization energies require more accurate electronic structure when those quantities are scientifically decisive.
