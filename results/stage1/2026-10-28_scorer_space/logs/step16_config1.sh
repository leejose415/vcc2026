#!/bin/bash
# 2026-10-29: Step 16 B on config 1 as well (config 1 out-scores the Step 16 carry held out); after the carried-arm B.
V=$HOME/vcc; S=$V/results/stage1/2026-10-28_scorer_space; L=$S/logs; PY=$HOME/venv/bin/python
say() { echo "[$(date -u +%H:%M:%S)] $*" | tee -a $L/step16_chain.log; }
until sed -n "/resume after/,\$p" $L/step16_chain.log | grep -qE "STEP16_DONE|STOPPED"; do sleep 20; done
$PY $V/code/vcc6/scorer_step16.py step16calls --magnitude config1 > $L/step16calls_config1.log 2>&1; rc=$?; say "step16calls config1 exit $rc"
say STEP16_CONFIG1_DONE
