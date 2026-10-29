#!/bin/bash
# 2026-10-28 scorer-space part 2, vcc2027: Step 10 proxies. H1 and HCT116 after the F1/F2 fold builds (their DE
# references exist); K562 and HEK293T after their DE references; then the fit. Proxies are pseudobulk + DE-table work
# (no cells, no table load). Exit codes checked.
V=$HOME/vcc; S=$V/results/stage1/2026-10-28_scorer_space; L=$S/logs; PY=$HOME/venv/bin/python
say() { echo "[$(date -u +%H:%M:%S)] $*" | tee -a $L/step10_chain.log; }
wait_for() { until grep -q "$2" $L/$1 2>/dev/null; do grep -q STOPPED $L/$1 2>/dev/null && { say "upstream STOPPED in $1"; exit 1; }; sleep 30; done; }
wait_for step10_folds.log FOLDS_DONE
for X in H1 HCT116; do ( $PY $V/code/vcc6/scorer_part2.py step10 --line $X > $L/step10_$X.log 2>&1; echo "EXIT $?" >> $L/step10_$X.log ) & done
wait; say "H1/HCT116 done: $(tail -1 $L/step10_H1.log) $(tail -1 $L/step10_HCT116.log)"
wait_for de_refs.log DE_REFS_DONE
for X in K562 HEK293T; do ( $PY $V/code/vcc6/scorer_part2.py step10 --line $X > $L/step10_$X.log 2>&1; echo "EXIT $?" >> $L/step10_$X.log ) & done
wait; say "K562/HEK293T done: $(tail -1 $L/step10_K562.log) $(tail -1 $L/step10_HEK293T.log)"
$PY $V/code/vcc6/scorer_part2.py step10fit > $L/step10fit.log 2>&1; say "fit exit $?"
say STEP10_DONE
