#!/usr/bin/env bash
# Resume Model 7 after the fits: H1 scoring -> cell benchmarks (2 streams) | A-C folders; then the two A-C files.
set -uo pipefail
V=$HOME/vcc; cd $V; L=$V/logs/m7; export PATH=$HOME/venv/bin:$PATH
PYTHONPATH=scripts/v5:analysis/model7 python analysis/model7/m7_abc.py model7 > $L/abc.log 2>&1 || echo "abc failed"
date -u; echo "steps 4-5"
PYTHONPATH=scripts/v5:analysis/model7 python analysis/model7/m7_score.py > $L/step45.log 2>&1 || { echo "score failed"; exit 1; }
bench() {
  cd $V/code/vcc_2026_3c
  python vcc3_generate.py benchmark --cells $V/lib/h1/H1_full_2025.h5ad --genes $V/data/gene_names.csv \
    --targets $V/lib/vcc_2026_3/validation/targets.csv --predictions $V/analysis/model7/h1_pred --eval-profile vcc2026 --n-targets 75 \
    --reference $V/lib/vcc_2026_3/val_final --variants $V/analysis/model7/$1 --out $V/analysis/model7/$2 --guide-col guide_id
}
date -u; echo "benchmarks"
bench variants_p1.json h1_bench_p1 > $L/bench_p1.log 2>&1 & p1=$!
bench variants_p2.json h1_bench_p2 > $L/bench_p2.log 2>&1 & p2=$!
wait $p1 $p2
cd $V; P=lib/vcc_2026_5/m7_pred/model7
date -u; echo "file floor"
GEN_EXTRA='--workers 6' bash scripts/make_submission.sh $P m7 $P/settings_floor.json vcc_2026_5_m7_floor > $L/file_floor.log 2>&1
date -u; echo "file no floor"
GEN_EXTRA='--workers 6' bash scripts/make_submission.sh $P m7 $P/settings_nofloor.json vcc_2026_5_m7_nofloor > $L/file_nofloor.log 2>&1
date -u; echo done
