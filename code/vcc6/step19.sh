#!/bin/bash
# Step 19 chain on vcc2026 (CPU only), parts strictly in order; each part logs to logs/<part>.log and writes
# logs/<part>.done on success. A part that fails stops the chain (later parts don't start).
# Usage: bash code/vcc6/step19.sh part0 partA ...      (Part D's grid runs 4 single-thread workers, then its summary)
set -u
cd ~/vcc
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 PYTHONUNBUFFERED=1
L=results/stage1/2026-11-01_step19/logs
PY=~/venv/bin/python
for part in "$@"; do
  [ -f $L/$part.done ] && { echo "$part already done"; continue; }
  echo "START $part $(date -u +%T)"
  if [ "$part" = partDgrid ]; then
    pids=()
    for w in 0 1 2 3; do
      OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 $PY code/vcc6/step19.py partD --worker $w > $L/partD_worker$w.log 2>&1 &
      pids+=($!)
    done
    rc=0; for p in "${pids[@]}"; do wait $p || rc=1; done
  else
    $PY code/vcc6/step19.py $part > $L/$part.log 2>&1; rc=$?
  fi
  if [ $rc -ne 0 ]; then echo "FAILED $part rc=$rc $(date -u +%T)"; exit $rc; fi
  touch $L/$part.done; echo "DONE $part $(date -u +%T)"
done
echo "CHAIN_END $(date -u +%T)"
