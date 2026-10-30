#!/bin/bash
# Step 15c chain for one config on one VM, strictly in sequence (no step overlaps a table load or a scoring):
# model (table load, own process) -> generate EvalA, HCT116final -> score both -> boot vs C+G -> diag (items 2-4).
# Usage: bash code/vcc6/step15c.sh CONFIG WORKERS     (config1p_K0 on vcc2026 with 4; config1p_K4_1200 on vcc2027 with 8)
set -uo pipefail
CFG=$1; W=$2
V=$HOME/vcc; P=$HOME/venv/bin/python; C=$V/code/vcc6; R=$V/results/stage1/2026-10-29_step15/c/$CFG; L=$R/logs
FE=$V/results/stage1/final_evals/harness/lib
mkdir -p $L
say() { echo "[$(date -u +%H:%M:%S)] $*" | tee -a $L/chain.log; }
stage() { local name=$1; shift; say "start $name"; "$@" > $L/$name.log 2>&1; local rc=$?; say "end $name: exit $rc"; [ $rc -eq 0 ] || { say "STOPPED at $name"; exit 1; }; }
declare -A CTRL=([EvalA]=$V/lib/harness/evalA/context.h5ad [HCT116final]=$FE/evalHCT116final/context.h5ad)
declare -A KD=([EvalA]=$V/results/target_lists/evalA_h1_292.txt [HCT116final]=$V/results/target_lists/final_HCT116.txt)
declare -A BG=([EvalA]=H1 [HCT116final]=HCT116)
declare -A REF=([EvalA]=$V/lib/harness/evalA/real.h5ad [HCT116final]=$FE/evalHCT116final/real.h5ad)
stage model $P $C/step15c.py model --config $CFG
for l in EvalA HCT116final; do
  [ -f $R/$l/e2e/salt_none/run.json ] || stage generate_$l env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    $P $C/run_line.py predict --controls ${CTRL[$l]} --knockdowns ${KD[$l]} --background ${BG[$l]} \
    --lfc-npz $R/$l/input.npz --salt "" --workers $W --out $R/$l/e2e
done
for l in EvalA HCT116final; do
  h=$R/$l/e2e/salt_none/harness/$(basename $(dirname ${REF[$l]}))
  [ -f $h/summary.json ] || stage score_$l $P $C/run_line.py score --pred $R/$l/e2e/salt_none/prediction.h5ad --reference ${REF[$l]} --threads $W
done
stage boot $P $C/step15c.py boot --config $CFG
say "SCORED_AND_BOOTED"
stage diag $P $C/step15c.py diag --config $CFG
say CHAIN_DONE
