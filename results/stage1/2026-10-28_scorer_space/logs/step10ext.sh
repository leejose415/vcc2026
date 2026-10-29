#!/bin/bash
# 2026-10-28: arm (e) past the Step 10 grid edge (non-called m 0.3-1, F0) — side check, not part of the carry rule.
# H1/HCT116 now; K562/HEK293T once their Step 11 runs exit (at most four proxy jobs at once).
V=$HOME/vcc; S=$V/results/stage1/2026-10-28_scorer_space; L=$S/logs; PY=$HOME/venv/bin/python
say() { echo "[$(date -u +%H:%M:%S)] $*" | tee -a $L/step10ext.log; }
for X in H1 HCT116; do ( $PY $V/code/vcc6/scorer_part2.py step10ext --line $X > $L/step10ext_$X.log 2>&1; echo "EXIT $?" >> $L/step10ext_$X.log ) & done
wait
until grep -q EXIT $L/step11_K562.log 2>/dev/null && grep -q EXIT $L/step11_HEK293T.log 2>/dev/null; do sleep 30; done
for X in K562 HEK293T; do ( $PY $V/code/vcc6/scorer_part2.py step10ext --line $X > $L/step10ext_$X.log 2>&1; echo "EXIT $?" >> $L/step10ext_$X.log ) & done
wait
say "lines: $(for X in K562 HCT116 HEK293T H1; do tail -1 $L/step10ext_$X.log; done | tr "\n" " ")"
$PY $V/code/vcc6/scorer_part2.py step10extfit > $L/step10extfit.log 2>&1; rc=$?; say "step10extfit exit $rc"
[ $rc -eq 0 ] && say STEP10EXT_DONE || say STOPPED
