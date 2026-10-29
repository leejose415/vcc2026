#!/bin/bash
# Step 4 on vcc2026: generate both arms (run_line --lfc-npz, frozen generator, salt ""), score each on Eval A (frozen
# harness), then the paired knockdown bootstrap. One job at a time; every exit code checked.
V=$HOME/vcc; P=$HOME/venv/bin/python; I=$V/results/stage1/2026-10-27_incumbent; O=$I/evalA
CTX=$V/lib/harness/evalA/context.h5ad; KD=$V/results/target_lists/evalA_h1_292.txt; REF=$V/lib/harness/evalA/real.h5ad
say() { echo "[$(date -u +%H:%M:%S)] $*" | tee -a $I/logs/step4.log; }
run() { local n=$1; shift; say "start $n"; "$@" > $I/logs/$n.log 2>&1; local rc=$?; say "end $n: exit $rc"; [ $rc -eq 0 ] || { say STOPPED; exit $rc; }; }
cd $V
for arm in M2 CG; do
  run gen_$arm $P code/vcc6/run_line.py predict --controls $CTX --knockdowns $KD --background H1 --lfc-npz $O/$arm/input.npz \
    --salt "" --workers 4 --out $O/$arm/e2e
done
for arm in M2 CG; do
  run score_$arm $P code/vcc6/run_line.py score --pred $O/$arm/e2e/salt_none/prediction.h5ad --reference $REF
done
run boot $P code/vcc6/boot_score.py $O/CG/e2e/salt_none/harness/evalA $O/M2/e2e/salt_none/harness/evalA --out $O/boot_M2_vs_CG.json
say STEP4_DONE
