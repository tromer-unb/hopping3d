#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="$HERE/results"
mkdir -p "$OUT"
hopping3d-eht w2o6 "$HERE/W2O6.cif" --cutoff 3.4 -o "$OUT/w2o6_eht.json"
echo "EHT parameters written to $OUT/w2o6_eht.json"
echo "For first-passage transport run: python3 reproduction/paper/atomic_cases.py"
