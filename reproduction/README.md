# Reproduction

This directory contains **code and structural inputs only** for reproducing the numerical datasets discussed in the article. Publication TeX and final figures are intentionally not stored in this repository.

## Current article workflow

```bash
bash reproduction/article_cases/run_all.sh
```

For production statistics:

```bash
FULL=1 python reproduction/article_cases/run_atomic_cases.py
FULL=1 python reproduction/article_cases/run_benzene_box.py
python reproduction/article_cases/run_interface.py
```

Generated JSON/CSV/CIF outputs are written to `reproduction/article_cases/results/` and are ignored by Git.

`inputs/` keeps the compact structural inputs/builders required by the atomic examples.
