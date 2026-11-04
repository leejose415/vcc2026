#!/bin/bash
# Part F combination rows on the proxy, one spec per process, in order.
cd ~/vcc; export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 PYTHONUNBUFFERED=1
L=results/stage1/2026-11-01_step19/logs
for spec in '{"name": "nosource_fix", "nosource_fix": true, "genin_check": true}' \
            '{"name": "hybrid_T5", "calls": "hybrid", "T": 5}' \
            '{"name": "hybrid_T5+nosource_fix", "calls": "hybrid", "T": 5, "nosource_fix": true}' \
            '{"name": "hybrid_T2+nosource_fix", "calls": "hybrid", "T": 2, "nosource_fix": true}'; do
  echo "START partF $spec $(date -u +%T)"
  ~/venv/bin/python code/vcc6/step19.py partF --spec "$spec" >> $L/partF.log 2>&1 || { echo "FAILED partF $(date -u +%T)"; exit 1; }
done
echo "DONE partF $(date -u +%T)"
