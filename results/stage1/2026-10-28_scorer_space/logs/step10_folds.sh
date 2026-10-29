#!/bin/bash
# 2026-10-28 scorer-space part 2, vcc2027: Step 10 source filters (F0 coverage only; F1, F2 per-fold C+G + Delta pool,
# one table load each), then the K562/HEK293T DE references. One job at a time; exit codes checked.
V=$HOME/vcc; S=$V/results/stage1/2026-10-28_scorer_space; L=$S/logs; PY=$HOME/venv/bin/python
say() { echo "[$(date -u +%H:%M:%S)] $*" | tee -a $L/step10_folds.log; }
run() { local n=$1; shift; say "start $n"; "$@" > $L/$n.log 2>&1; local rc=$?; say "end $n: exit $rc"; [ $rc -eq 0 ] || { say STOPPED; exit $rc; }; }
run cov_F0 $PY $V/code/vcc6/scorer_study.py folds --filter F0 --coverage-only
run folds_F1 $PY $V/code/vcc6/scorer_study.py folds --filter F1
run folds_F2 $PY $V/code/vcc6/scorer_study.py folds --filter F2
say FOLDS_DONE
bash $L/de_refs.sh
