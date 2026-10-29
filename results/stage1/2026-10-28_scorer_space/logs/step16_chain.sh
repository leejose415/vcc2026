#!/bin/bash
# 2026-10-29 Step 16 (Jose): wider config-1 magnitude grid (A) on the four lines in parallel, its fit, then the H1
# half-split call budgets (B). Proxy only, no cells.
V=$HOME/vcc; S=$V/results/stage1/2026-10-28_scorer_space; L=$S/logs; PY=$HOME/venv/bin/python
say() { echo "[$(date -u +%H:%M:%S)] $*" | tee -a $L/step16_chain.log; }
say start
for X in K562 HCT116 HEK293T H1; do ( $PY $V/code/vcc6/scorer_step16.py step16 --line $X > $L/step16_$X.log 2>&1; echo "EXIT $?" >> $L/step16_$X.log ) & done
wait
say "lines: $(for X in K562 HCT116 HEK293T H1; do tail -1 $L/step16_$X.log; done | tr "\n" " ")"
for X in K562 HCT116 HEK293T H1; do tail -1 $L/step16_$X.log | grep -q "EXIT 0" || { say STOPPED; exit 1; }; done
$PY $V/code/vcc6/scorer_step16.py step16fit > $L/step16fit.log 2>&1; rc=$?; say "step16fit exit $rc"; [ $rc -eq 0 ] || { say STOPPED; exit $rc; }
$PY $V/code/vcc6/scorer_step16.py step16calls > $L/step16calls.log 2>&1; rc=$?; say "step16calls exit $rc"; [ $rc -eq 0 ] || { say STOPPED; exit $rc; }
say STEP16_DONE
