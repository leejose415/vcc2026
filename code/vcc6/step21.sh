#!/bin/bash
# Step 21 Track 2 chain on vcc2026 (Jose 2026-11-02): runs the named parts in order, each in its own process, stopping
# at the first failure. One chain at a time on this VM, so no table load runs next to scoring.
# usage: bash code/vcc6/step21.sh TAG part [part ...]      (TAG names this chain's log: chain_<TAG>.log)
set -u
cd ~/vcc
TAG=$1; shift
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 POLARS_MAX_THREADS=4
L=results/stage1/2026-11-02_step21/logs; mkdir -p $L
( while true; do echo "$(date -u +%T) $(free -m | awk '/Mem:/ {print "used " $3 " MB, available " $7 " MB"} /Swap:/ {print "swap used " $3 " MB"}' | tr '\n' ' ') $(ps -eo rss,args --sort=-rss | awk 'NR==2 {printf "top rss %.1f GB %s %s", $1/1e6, $3, $4}')"; sleep 60; done ) >> $L/mem_monitor_vcc2026.log 2>&1 &
MON=$!
trap 'kill $MON 2>/dev/null' EXIT
for p in "$@"; do
  [ -e $L/$p.done ] && { echo "$(date -u +%T) SKIP $p (done)" >> $L/chain_$TAG.log; continue; }
  echo "$(date -u +%T) START $p" >> $L/chain_$TAG.log
  ~/venv/bin/python code/vcc6/step21.py $p > $L/$p.log 2>&1
  rc=$?
  if [ $rc -ne 0 ]; then echo "$(date -u +%T) FAILED $p rc $rc" >> $L/chain_$TAG.log; exit $rc; fi
  touch $L/$p.done; echo "$(date -u +%T) DONE $p" >> $L/chain_$TAG.log
done
echo "$(date -u +%T) CHAIN END $TAG" >> $L/chain_$TAG.log
