#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="$HERE/results"
mkdir -p "$OUT"
hopping3d-eht bn "$HERE/BN_bulk.cif" --repeat 2 2 2 -o "$OUT/bn_eht.json"
echo "EHT parameters written to $OUT/bn_eht.json"
echo "For the diffusion calculation run: python3 reproduction/paper/atomic_cases.py"
