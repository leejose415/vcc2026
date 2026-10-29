#!/bin/bash
# 2026-10-28 D/E/F final file, C+G (stage1 arm C + gaussian), Jose's choice. One context at a time on vcc2026, salt "",
# each context's own controls, the 300 final targets; every exit code checked.
V=$HOME/vcc; P=$HOME/venv/bin/python; D=$V/results/stage1/2026-10-28_DEF
say() { echo "[$(date -u +%H:%M:%S)] $*" | tee -a $D/logs/run_def.log; }
run() { local n=$1; shift; say "start $n"; "$@" > $D/logs/$n.log 2>&1; local rc=$?; say "end $n: exit $rc"; [ $rc -eq 0 ] || { say STOPPED; exit $rc; }; }
cd $V
for X in D E F; do
  run predict_$X $P code/vcc6/run_line.py predict --controls $V/data_final/extracted/context_$X.h5ad \
    --knockdowns $V/results/target_lists/final_DEF_300.txt --background $X \
    --effects-dir $V/lib/effects_v5 --effects-dir $V/lib/effects_v5_extra --arm C --magnitude gaussian --salt "" --workers 4 --out $D/$X
done
say PREDICT_DONE
