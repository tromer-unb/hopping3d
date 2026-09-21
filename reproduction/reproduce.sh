#!/usr/bin/env bash
set -euo pipefail
MODE=${1:-quick}
export MPLBACKEND=Agg
if [[ "$MODE" == "full" ]]; then export FULL=1; else export FULL=0; fi
python reproduction/inputs/build_layered_carbon.py
python reproduction/scripts/figure1_overview.py
python reproduction/scripts/figure2_topology.py
python reproduction/scripts/figure2_kinetics.py
python reproduction/scripts/figure2_paths.py
python reproduction/scripts/figure3_bn_tensor.py
python reproduction/scripts/figure4_chemistry.py
python reproduction/scripts/plot_publication_figures.py
python reproduction/scripts/figureS_controls.py
python reproduction/scripts/figureS_convergence.py
