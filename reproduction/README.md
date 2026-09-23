# Reproduction

This directory contains **code and structural inputs only** for reproducing the numerical datasets discussed in the article. Publication TeX and final figures are intentionally not stored in this repository.

## Current article workflow

```bash
bash reproduction/run_all.sh
```

For production statistics:

```bash
FULL=1 python3 reproduction/paper/atomic_cases.py
FULL=1 python3 reproduction/paper/benzene_box.py
python3 reproduction/paper/graphene_benzene.py
```

Generated JSON/CSV/CIF outputs are written to `reproduction/paper/results/` and are ignored by Git.

All versioned structural inputs live next to their examples under `examples/paper/`.
