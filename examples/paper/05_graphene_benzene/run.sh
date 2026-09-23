#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="$HERE/results"
mkdir -p "$OUT"
hopping3d-eht interface "$HERE/graphene_benzene.cif" -o "$OUT/interface.json" --report "$OUT/interface.md"
hopping3d-eht coverage "$HERE/graphene_4benzene.cif" -o "$OUT/coverage.json" --report "$OUT/coverage.md"
echo "Wrote single- and multi-adsorbate EHT reports in $OUT"
