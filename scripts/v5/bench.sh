#!/usr/bin/env bash
# One H1 cell-based benchmark process on the submission path (guide-stratified controls), reusing val_final.
# Usage: bash bench.sh <variants json> <out folder name>
set -euo pipefail
V=$HOME/vcc; RUN=$V/lib/vcc_2026_4/runs/dev_nopie
cd $V/code/vcc_2026_3c
export PATH=$HOME/venv/bin:$PATH
date -u
python vcc3_generate.py benchmark --cells $V/lib/h1/H1_full_2025.h5ad --genes $V/data/gene_names.csv \
  --targets $V/lib/vcc_2026_3/validation/targets.csv --predictions $RUN/h1_pred_v5 --eval-profile vcc2026 --n-targets 75 \
  --reference $V/lib/vcc_2026_3/val_final --variants "$1" --out $RUN/"$2" --guide-col guide_id
date -u
