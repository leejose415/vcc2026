#!/usr/bin/env bash
# v6-lin main through the submission cell path on H1 (guide-stratified controls, val_final reference): floor on (as the
# submissions) and floor off.  Usage: bash analysis/v6lin/bench_v6lin.sh
set -euo pipefail
V=$HOME/vcc; P=$V/analysis/v6lin
cat > $P/variants.json <<'EOF'
[{"method": "v6lin", "call": "z", "z": 0.0, "min_effect": 0.0, "alpha": 1.0, "detect_floor": true, "label": "v6lin_floor"},
 {"method": "v6lin", "call": "z", "z": 0.0, "min_effect": 0.0, "alpha": 1.0, "detect_floor": false, "label": "v6lin_nofloor"}]
EOF
cd $V/code/vcc_2026_3c; export PATH=$HOME/venv/bin:$PATH
python vcc3_generate.py benchmark --cells $V/lib/h1/H1_full_2025.h5ad --genes $V/data/gene_names.csv \
  --targets $V/lib/vcc_2026_3/validation/targets.csv --predictions $P/h1_pred --eval-profile vcc2026 --n-targets 75 \
  --reference $V/lib/vcc_2026_3/val_final --variants $P/variants.json --out $P/h1_bench --guide-col guide_id
