#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="$HERE/results"
mkdir -p "$OUT"
hopping3d-eht graphene "$HERE/stacked_graphene_4layers.cif" -o "$OUT/carbon_eht.json"
echo "Wrote $OUT/carbon_eht.json"
