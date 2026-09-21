#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
pdflatex -interaction=nonstopmode main.tex >/dev/null
bibtex main >/dev/null
pdflatex -interaction=nonstopmode main.tex >/dev/null
pdflatex -interaction=nonstopmode main.tex >/dev/null
pdflatex -interaction=nonstopmode SI.tex >/dev/null
bibtex SI >/dev/null || true
pdflatex -interaction=nonstopmode SI.tex >/dev/null
pdflatex -interaction=nonstopmode SI.tex >/dev/null
