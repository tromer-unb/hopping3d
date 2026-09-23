#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
PYTHON="${PYTHON:-python3}"
$PYTHON reproduction/paper/atomic_cases.py
$PYTHON reproduction/paper/benzene_box.py
$PYTHON reproduction/paper/graphene_benzene.py
printf '\nDone. Numerical outputs are in reproduction/paper/results/.\n'
