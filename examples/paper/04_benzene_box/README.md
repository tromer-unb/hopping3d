# 04 - Periodic box of 100 benzene molecules

Run the quick benchmark:

```bash
bash examples/paper/04_benzene_box/run.sh
```

Production statistics:

```bash
FULL=1 python3 reproduction/paper/benzene_box.py
```

The script builds a periodic 5x5x4 molecular box (100 molecules), identifies 300 nearest-neighbor molecular pairs, solves one EHT problem per dimer, and extracts both HOMO-like (hole) and LUMO-like (electron) transfer integrals.

The paper transport demonstration uses the hole channel with Marcus kinetics and an external `lambda = 0.20 eV`. Nuclear reorganization energy is not supplied by ordinary EHT and must be provided or refined independently.
