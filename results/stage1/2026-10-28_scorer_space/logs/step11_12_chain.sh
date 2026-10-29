#!/bin/bash
# 2026-10-28 scorer-space part 2, vcc2027: after Step 10's fit, Step 11 on K562/HEK293T (H1/HCT116 started at 18:42 on the
# provisional carried config, which the fit must confirm), Step 11's fit, Step 12 on all lines, Step 12's fit.
V=$HOME/vcc; S=$V/results/stage1/2026-10-28_scorer_space; L=$S/logs; PY=$HOME/venv/bin/python
say() { echo "[$(date -u +%H:%M:%S)] $*" | tee -a $L/step11_12_chain.log; }
until grep -q STEP10_DONE $L/step10_chain.log 2>/dev/null; do grep -qE "STOPPED|exit [1-9]" $L/step10_chain.log 2>/dev/null && { say "step10 chain failed"; exit 1; }; sleep 30; done
$PY -c "
import json; c=json.load(open('$S/tables/step10.json'))['carried']
ok = c['filter']=='F0' and c['arm']=='a' and abs(c['m']-0.5)<1e-9
print('carried', c, 'matches provisional:', ok); raise SystemExit(0 if ok else 3)" >> $L/step11_12_chain.log 2>&1 || { say "carried differs from the provisional F0/a/0.5: Step 11 H1/HCT116 must be rerun; STOPPED"; exit 3; }
for X in K562 HEK293T; do ( $PY $V/code/vcc6/scorer_step11.py step11 --line $X > $L/step11_$X.log 2>&1; echo "EXIT $?" >> $L/step11_$X.log ) & done
wait
until grep -q EXIT $L/step11_H1.log && grep -q EXIT $L/step11_HCT116.log; do sleep 30; done
say "step11 lines: $(for X in K562 HCT116 HEK293T H1; do tail -1 $L/step11_$X.log; done | tr '\n' ' ')"
$PY $V/code/vcc6/scorer_step11.py step11fit > $L/step11fit.log 2>&1; rc=$?; say "step11fit exit $rc"; [ $rc -eq 0 ] || { say STOPPED; exit $rc; }
for X in K562 HCT116 HEK293T H1; do ( $PY $V/code/vcc6/scorer_step12.py step12 --line $X > $L/step12_$X.log 2>&1; echo "EXIT $?" >> $L/step12_$X.log ) & done
wait
say "step12 lines: $(for X in K562 HCT116 HEK293T H1; do tail -1 $L/step12_$X.log; done | tr '\n' ' ')"
$PY $V/code/vcc6/scorer_step12.py step12fit > $L/step12fit.log 2>&1; rc=$?; say "step12fit exit $rc"; [ $rc -eq 0 ] || { say STOPPED; exit $rc; }
say STEP12_DONE
