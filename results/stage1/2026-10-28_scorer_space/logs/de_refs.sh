#!/bin/bash
# 2026-10-28 scorer-space part 2, vcc2027: full DE references (harness finish: real cache, ceiling, generic baseline,
# 5 replicate-anchor splits) for K562final and HEK293Tfinal on their existing cells, so all six components can be
# proxied on every fold line. One at a time; exit codes checked.
V=$HOME/vcc; L=$V/results/stage1/2026-10-28_scorer_space/logs
say() { echo "[$(date -u +%H:%M:%S)] $*" | tee -a $L/de_refs.log; }
for X in K562final HEK293Tfinal; do
  say "start $X"
  $HOME/venv/bin/python $V/code/vcc6/final_evals.py build $X --finish-only --threads 8 > $L/de_ref_$X.log 2>&1
  rc=$?; say "end $X: exit $rc"; [ $rc -eq 0 ] || { say STOPPED; exit $rc; }
done
say DE_REFS_DONE
