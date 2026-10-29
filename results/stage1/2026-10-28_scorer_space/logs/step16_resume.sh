#!/bin/bash
# 2026-10-29: step16fit failed on a spreads() bug (whole fits dict passed); rows are complete, rerun fit + calls.
V=$HOME/vcc; S=$V/results/stage1/2026-10-28_scorer_space; L=$S/logs; PY=$HOME/venv/bin/python
say() { echo "[$(date -u +%H:%M:%S)] $*" | tee -a $L/step16_chain.log; }
say "resume after the spreads() fix"
$PY $V/code/vcc6/scorer_step16.py step16fit > $L/step16fit.log 2>&1; rc=$?; say "step16fit exit $rc"; [ $rc -eq 0 ] || { say STOPPED; exit $rc; }
$PY $V/code/vcc6/scorer_step16.py step16calls > $L/step16calls.log 2>&1; rc=$?; say "step16calls exit $rc"; [ $rc -eq 0 ] || { say STOPPED; exit $rc; }
say STEP16_DONE
