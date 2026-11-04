#!/usr/bin/env bash
# Model 7 end to end on vcc2027 after the Step 0 tables.
#   branch A: diagnostics -> validation fit -> Model 2 reference -> H1 pseudobulk scoring -> H1 cell benchmarks (2 streams)
#   branch B: final fit (H1 included) -> A-C prediction folders; then, after branch A's benchmarks (vcc prep peaks near
#             30 GB, so never alongside cell-eval), the A-C leaderboard files with the floor on, then off.
set -uo pipefail
V=$HOME/vcc; cd $V; L=$V/logs/m7; export PATH=$HOME/venv/bin:$PATH
for n in K562 RPE1 HepG2 Jurkat H1 HCT116 HEK293T CD4_Rest CD4_Stim8hr CD4_Stim48hr; do
  until [ -f lib/vcc_2026_5/m7/${n}_report.json ]; do sleep 20; done
done
rm -f $L/benches_done
branch_b() {
  cd $V; date -u; echo "B: final fit"
  PYTHONPATH=scripts/v5:scripts python analysis/model7/m7_fit.py final > $L/step2_final.log 2>&1 || { echo "B: final failed"; return 1; }
  PYTHONPATH=scripts/v5:analysis/model7 python analysis/model7/m7_abc.py model7 > $L/abc.log 2>&1 || { echo "B: abc failed"; return 1; }
  date -u; echo "B: waiting for the benchmarks"; until [ -f $L/benches_done ]; do sleep 30; done
  P=lib/vcc_2026_5/m7_pred/model7
  date -u; echo "B: file floor"
  GEN_EXTRA='--workers 6' bash scripts/make_submission.sh $P m7 $P/settings_floor.json vcc_2026_5_m7_floor > $L/file_floor.log 2>&1
  date -u; echo "B: file no floor"
  GEN_EXTRA='--workers 6' bash scripts/make_submission.sh $P m7 $P/settings_nofloor.json vcc_2026_5_m7_nofloor > $L/file_nofloor.log 2>&1
  date -u; echo "B: done"
}
branch_b > $L/branch_b.log 2>&1 &
date -u; echo "A: step 1"
PYTHONPATH=scripts/v5:scripts python analysis/model7/m7_fit.py diag > $L/step1.log 2>&1
date -u; echo "A: step 2 validation"
PYTHONPATH=scripts/v5:scripts python analysis/model7/m7_fit.py validate > $L/step2_validation.log 2>&1 || { echo "A: validate failed"; touch $L/benches_done; wait; exit 1; }
date -u; echo "A: model 2 reference"
PYTHONPATH=analysis/v6lin:scripts/v5 python analysis/model7/m2_dump.py > $L/m2_dump.log 2>&1
date -u; echo "A: steps 4-5"
PYTHONPATH=scripts/v5:analysis/model7 python analysis/model7/m7_score.py > $L/step45.log 2>&1 || { echo "A: score failed"; touch $L/benches_done; wait; exit 1; }
bench() {
  cd $V/code/vcc_2026_3c
  python vcc3_generate.py benchmark --cells $V/lib/h1/H1_full_2025.h5ad --genes $V/data/gene_names.csv \
    --targets $V/lib/vcc_2026_3/validation/targets.csv --predictions $V/analysis/model7/h1_pred --eval-profile vcc2026 --n-targets 75 \
    --reference $V/lib/vcc_2026_3/val_final --variants $V/analysis/model7/$1 --out $V/analysis/model7/$2 --guide-col guide_id
}
date -u; echo "A: H1 cell benchmarks"
bench variants_p1.json h1_bench_p1 > $L/bench_p1.log 2>&1 & p1=$!
bench variants_p2.json h1_bench_p2 > $L/bench_p2.log 2>&1 & p2=$!
wait $p1 $p2
touch $L/benches_done; date -u; echo "A: done"
wait
