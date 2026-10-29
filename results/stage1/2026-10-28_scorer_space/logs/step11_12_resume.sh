#!/bin/bash
# 2026-10-28 21:30: resume after step11fit failed (HGNC table missing on vcc2027; copied from vcc2026, md5 checked).
# Waits for the Step 10 extension so no more than four proxy jobs run at once.
V=$HOME/vcc; S=$V/results/stage1/2026-10-28_scorer_space; L=$S/logs; PY=$HOME/venv/bin/python
say() { echo "[$(date -u +%H:%M:%S)] $*" | tee -a $L/step11_12_chain.log; }
until grep -qE "STEP10EXT_DONE|STOPPED" $L/step10ext.log 2>/dev/null; do sleep 30; done
say "resume: step10ext $(tail -1 $L/step10ext.log)"
$PY $V/code/vcc6/scorer_step11.py step11fit > $L/step11fit.log 2>&1; rc=$?; say "step11fit exit $rc"; [ $rc -eq 0 ] || { say STOPPED; exit $rc; }
for X in K562 HCT116 HEK293T H1; do ( $PY $V/code/vcc6/scorer_step12.py step12 --line $X > $L/step12_$X.log 2>&1; echo "EXIT $?" >> $L/step12_$X.log ) & done
wait
say "step12 lines: $(for X in K562 HCT116 HEK293T H1; do tail -1 $L/step12_$X.log; done | tr "\n" " ")"
$PY $V/code/vcc6/scorer_step12.py step12fit > $L/step12fit.log 2>&1; rc=$?; say "step12fit exit $rc"; [ $rc -eq 0 ] || { say STOPPED; exit $rc; }
say STEP12_DONE
