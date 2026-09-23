# Reproducibility

This repository reproduces the **numerical datasets** used by the current Hopping3D study. It intentionally does not contain manuscript TeX or publication figures.

## Quick verification

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
pytest -q
bash reproduction/article_cases/run_all.sh
```

## Production statistics

```bash
FULL=1 python3 reproduction/article_cases/run_atomic_cases.py
FULL=1 python3 reproduction/article_cases/run_benzene_box.py
python3 reproduction/article_cases/run_interface.py
```

Generated files are written to `reproduction/article_cases/results/` and are ignored by Git. Random seeds and ensemble sizes are explicit in the scripts. See `docs/article_reproduction.md` for the mapping between each script and the physical case.
