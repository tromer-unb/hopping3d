# Article source, data, and figures

This directory contains the revised EHT-anchored manuscript source, frozen article-local datasets, and the deterministic workflow used to generate the publication figures.

The numerical code paths used to generate the article datasets are kept under `../examples/paper/` and `../reproduction/paper/`. The files in `paper/data/` are frozen inputs/processed datasets used to reconstruct the publication figures; they are intentionally versioned with the manuscript. Generated PNG/SVG files are build artifacts and are not versioned.

Rebuild the figures from the frozen data with:

```bash
python3 scripts/build_figures.py
```

Compile the article and Supplemental Material with:

```bash
bash build.sh
```

Build artifacts such as `.aux`, `.log`, `.bbl`, and compiled PDFs are not versioned. See `../docs/article_reproduction.md` for the mapping between article cases and executable reproduction scripts.
