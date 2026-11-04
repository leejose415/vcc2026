#!/bin/bash
# Step 17 on vcc2027: waits until this VM's 15c chain has finished generation, scoring and boot_score
# (SCORED_AND_BOOTED in its chain.log), then predictor -> rows (H1 in 4 parts, the guardrails 1 each) -> summary.
# Usage: bash code/vcc6/step17.sh
set -uo pipefail
V=$HOME/vcc; P=$HOME/venv/bin/python; C=$V/code/vcc6; R=$V/results/stage1/2026-10-30_step17; L=$R/logs
CHAIN=$V/results/stage1/2026-10-29_step15/c/config1p_K4_1200/logs/chain.log
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 POLARS_MAX_THREADS=1
mkdir -p $L
say() { echo "[$(date -u +%H:%M:%S)] $*" | tee -a $L/step17.log; }
say "waiting for 15c (SCORED_AND_BOOTED)"
until grep -q SCORED_AND_BOOTED $CHAIN 2>/dev/null; do grep -q STOPPED $CHAIN 2>/dev/null && { say "15c chain STOPPED; not starting"; exit 1; }; sleep 30; done
say "start predictor"; $P $C/step17.py predictor > $L/predictor.log 2>&1 || { say "STOPPED at predictor"; exit 1; }; say "end predictor"
pids=()
for p in 0 1 2 3; do $P $C/step17.py rows --line H1 --part $p --nparts 4 > $L/rows_H1_p$p.log 2>&1 & pids+=($!); done
for g in K562 HCT116 HEK293T; do $P $C/step17.py rows --line $g > $L/rows_$g.log 2>&1 & pids+=($!); done
rc=0; for x in "${pids[@]}"; do wait $x || rc=1; done
[ $rc -eq 0 ] || { say "STOPPED at rows"; exit 1; }; say "end rows"
$P $C/step17.py summary > $L/summary.log 2>&1 || { say "STOPPED at summary"; exit 1; }
say STEP17_DONE
