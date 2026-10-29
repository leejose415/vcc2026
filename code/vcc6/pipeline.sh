#!/bin/bash
# Stage 1 RPE1 dev run (results/stage1/2026-10-27_rpe1/), after Step 0 (stage1.py, the 300 list, the RPE1dev eval).
# Lane 1: Step 1 arms -> Step 1 analysis -> Step 3 inputs + profile check -> run_line zero / oracle (3 salts each)
#         -> [waits for Step 2] -> run_line B / C (3 salts each); every prediction scored as it lands.
# Lane 2: Step 2 (one load of the 117 tables, ~15 GB), started with lane 1.
# Stops at the first failing stage. Usage: bash code/vcc6/pipeline.sh
set -uo pipefail
V=$HOME/vcc; P=$HOME/venv/bin/python; R=$V/results/stage1/2026-10-27_rpe1; L=$R/logs; C=$V/code/vcc6
EVAL=$R/harness/lib/evalRPE1dev; LIST=$V/results/target_lists/rpe1_dev_300.txt
DIRS="--effects-dir $V/lib/effects_v5 --effects-dir $V/lib/effects_v5_extra"
SALTS=(--salt "" --salt s1 --salt s2)
mkdir -p $L
say() { echo "[$(date -u +%H:%M:%S)] $*" | tee -a $L/pipeline.log; }
stage() { local name=$1; shift; say "start $name"; "$@" > $L/$name.log 2>&1; local rc=$?; say "end $name: exit $rc"; [ $rc -eq 0 ] || { say "STOPPED at $name"; exit 1; }; }

( $P $C/step2_votes.py > $L/step2.log 2>&1; echo "EXIT=$?" >> $L/step2.log; say "step2 lane: $(tail -1 $L/step2.log)" ) &
stage step1_arms $P $C/step1_oracle.py arms --workers 4
stage step1_analyze $P $C/step1_analyze.py
stage step3_inputs $P $C/step3_inputs.py inputs
stage step3_profile $P $C/step3_inputs.py profile
for arm in zero oracle; do
  stage step3_${arm}_predict $P $C/run_line.py predict --controls $EVAL/context.h5ad --knockdowns $LIST --background RPE1 \
    --lfc-npz $R/step3/inputs/$arm.npz "${SALTS[@]}" --out $R/step3/$arm
  for s in none s1 s2; do stage step3_${arm}_${s}_score $P $C/run_line.py score --pred $R/step3/$arm/salt_$s/prediction.h5ad --reference $EVAL/real.h5ad; done
done
until grep -q '^EXIT=' $L/step2.log 2>/dev/null; do sleep 30; done
grep -q '^EXIT=0' $L/step2.log || { say "STOPPED: step2 failed"; exit 1; }
for arm in B C; do
  stage step3_${arm}_predict $P $C/run_line.py predict --controls $EVAL/context.h5ad --knockdowns $LIST --background RPE1 \
    $DIRS --arm $arm "${SALTS[@]}" --out $R/step3/arm$arm
  for s in none s1 s2; do stage step3_${arm}_${s}_score $P $C/run_line.py score --pred $R/step3/arm$arm/salt_$s/prediction.h5ad --reference $EVAL/real.h5ad; done
done
say PIPELINE_DONE
