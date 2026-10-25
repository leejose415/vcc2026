#!/usr/bin/env bash
# H1 cell-based benchmark of analysis/final/h1_pred arms (submission path: guide-stratified controls, val_final).
# Usage: bash analysis/final/bench.sh <variants json> <out folder name> [extra generator switches, e.g. --exact-totals --match-pool]
set -euo pipefail
V=$HOME/vcc; F=$V/analysis/final
cd $V/code/vcc_2026_3c; export PATH=$HOME/venv/bin:$PATH
python ${GEN_SCRIPT:-vcc3_generate.py} benchmark --cells $V/lib/h1/H1_full_2025.h5ad --genes $V/data/gene_names.csv \
  --targets $V/lib/vcc_2026_3/validation/targets.csv --predictions $F/h1_pred --eval-profile vcc2026 --n-targets 75 \
  --reference $V/lib/vcc_2026_3/val_final --variants $F/$1 --out $F/$2 --guide-col guide_id "${@:3}"
