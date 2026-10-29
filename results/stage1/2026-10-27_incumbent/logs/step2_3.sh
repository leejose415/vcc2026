#!/bin/bash
# Steps 2-3 on vcc2026, one job at a time, every exit code checked.
V=$HOME/vcc; P=$HOME/venv/bin/python; I=$V/results/stage1/2026-10-27_incumbent; O=$I/evalA
CTX=$V/lib/harness/evalA/context.h5ad; KD=$V/results/target_lists/evalA_h1_292.txt
say() { echo "[$(date -u +%H:%M:%S)] $*" | tee -a $I/logs/step2_3.log; }
run() { local n=$1; shift; say "start $n"; "$@" > $I/logs/$n.log 2>&1; local rc=$?; say "end $n: exit $rc"; [ $rc -eq 0 ] || { say STOPPED; exit $rc; }; }
cd $V
run cg_model $P code/vcc6/model2_port.py cg --controls $CTX --knockdowns $KD --background H1 --out $O
run m2_port $P code/vcc6/model2_port.py port --controls $CTX --knockdowns $KD --holdout H1 --out $O
run step3 $P code/vcc6/model2_port.py compare --out $O
say STEPS_2_3_DONE
