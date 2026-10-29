#!/bin/bash
# 2026-10-29: side check past the B budget edge (K4/K5 at 1200 won), config 1 magnitudes; after the config-1 B run.
V=$HOME/vcc; S=$V/results/stage1/2026-10-28_scorer_space; L=$S/logs; PY=$HOME/venv/bin/python
say() { echo "[$(date -u +%H:%M:%S)] $*" | tee -a $L/step16_chain.log; }
until grep -q STEP16_CONFIG1_DONE $L/step16_chain.log; do sleep 20; done
$PY $V/code/vcc6/scorer_step16.py step16calls --magnitude config1 --budgets 1600,2400,3200,4800 > $L/step16calls_config1_ext.log 2>&1; rc=$?; say "step16calls config1 ext exit $rc"
say STEP16_EXT_DONE
