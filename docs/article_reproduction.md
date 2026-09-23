# Reproducing the numerical datasets used in the article

This repository intentionally contains **code and structural inputs only**. Manuscript TeX and publication figures are not stored here.

Install the package, then run:

```bash
bash reproduction/article_cases/run_all.sh
```

The scripts reproduce the numerical JSON/CSV datasets behind the article analyses:

1. layered carbon EHT interlayer scale;
2. BN direct-EHT diffusion tensor;
3. W2O6 EHT chemical-path benchmark;
4. periodic 100-benzene box with aligned/isotropic orientations;
5. single- and multi-benzene graphene interface parameterization.

For a faster smoke test use the default settings. For the production statistics used in the article:

```bash
FULL=1 python3 reproduction/article_cases/run_atomic_cases.py
FULL=1 python3 reproduction/article_cases/run_benzene_box.py
python3 reproduction/article_cases/run_interface.py
```

Outputs are written under `reproduction/article_cases/results/`. Plotting is deliberately left outside the repository: these files are the reproducible numerical source data.
