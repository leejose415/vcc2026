#!/bin/bash
# 2026-10-28 08:20: the queue was stopped by mistake (to move F to vcc2027, which Jose declined). Wait for E (pid 457463),
# check its output, then run F here as before. Every exit code checked.
V=$HOME/vcc; P=$HOME/venv/bin/python; D=$V/results/stage1/2026-10-28_DEF; KD=$V/results/target_lists/final_DEF_300.txt
say() { echo "[$(date -u +%H:%M:%S)] $*" | tee -a $D/logs/run_def.log; }
while kill -0 457463 2>/dev/null; do sleep 20; done
if grep -q "generated 300/300" $D/logs/gen_E.log && [ -s $D/E/e2e/salt_none/run.json ]; then say "end gen_E: output complete (300/300, run.json written)"; else say "gen_E: output incomplete"; say STOPPED; exit 1; fi
say "start gen_F"
cd $V && $P code/vcc6/run_line.py predict --controls $V/data_final/extracted/context_F.h5ad --knockdowns $KD --background F   --lfc-npz $D/F/input.npz --salt "" --workers 4 --out $D/F/e2e > $D/logs/gen_F.log 2>&1
rc=$?; say "end gen_F: exit $rc"; [ $rc -eq 0 ] || { say STOPPED; exit $rc; }
say PREDICT_DONE
