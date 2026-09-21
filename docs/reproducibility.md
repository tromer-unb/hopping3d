# Reproducibility guide

## Minimal calculations

1. `pip install -e .`
2. `hopping3d examples/first_passage/params.json`
3. `hopping3d examples/diffusion_tensor/params.json`

## Article figures

- `bash reproduction/reproduce.sh quick` runs reduced-statistics simulations suitable for verification.
- `bash reproduction/reproduce.sh full` uses the production ensemble sizes used for the reported statistics.
- `python reproduction/scripts/plot_publication_figures.py` regenerates Figures 2--4 directly from the committed reference summaries, without rerunning long ensembles.

Random seeds are explicit in every production script. Reference CSV files contain the reported ensemble summaries so numerical comparisons can be automated.

## Manuscript source

After figures are present in `reproduction/generated/`, run `make manuscript` (or `bash manuscript/build.sh`) to compile the LaTeX source. The manuscript references generated figures by relative path, so the repository contains a direct chain from CIF/parameters to data, figures, and text.
