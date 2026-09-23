# Periodic box of 100 benzene molecules

The article benchmark uses a periodic 5x5x4 center lattice containing 100 benzene molecules. The aligned and orientationally isotropic cases use identical molecular centers; only orientations differ.

Generate and run the quick benchmark:

```bash
python3 reproduction/article_cases/run_benzene_box.py
```

Production statistics:

```bash
FULL=1 python3 reproduction/article_cases/run_benzene_box.py
```

For every molecular neighbor pair, calibrated pi-EHT is projected onto the benzene HOMO and LUMO doublets, producing separate hole and electron transfer integrals. The article transport benchmark uses the hole channel with Marcus kinetics and an external `lambda = 0.20 eV`. Ordinary EHT does not determine the nuclear reorganization energy.
