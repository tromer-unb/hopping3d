#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"
${PYTHON:-python3} reproduction/article_cases/run_atomic_cases.py
${PYTHON:-python3} reproduction/article_cases/run_benzene_box.py
${PYTHON:-python3} reproduction/article_cases/run_interface.py
printf '\nArticle numerical datasets reproduced in reproduction/article_cases/results/\n'
