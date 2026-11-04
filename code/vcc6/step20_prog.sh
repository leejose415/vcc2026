#!/bin/bash
# Step 20 Part F chain on vcc2027 (Jose 2026-11-02). F3c's disk check first (the ARCHS4 download would run in the
# background, but the check decides first); then F1, F2, F3a, F3b, F6, F3 summary, F3d, F4, F5 (one timed row, then
# two 4-thread workers), F5 summary. A failed download stops only its sub-part (F3a / F3b); anything else stops the chain.
set -u
cd ~/vcc
L=results/stage1/2026-11-02_step20/logs; mkdir -p $L results/stage1/2026-11-02_step20/tables
export POLARS_MAX_THREADS=8
( while true; do echo "$(date -u +%T) $(free -m | awk '/Mem:/ {print "used " $3 " MB, available " $7 " MB"} /Swap:/ {print "swap used " $3 " MB"}' | tr '\n' ' ') $(ps -eo rss,args --sort=-rss | awk 'NR==2 {printf "top rss %.1f GB %s %s", $1/1e6, $3, $4}')"; sleep 60; done ) >> $L/mem_monitor_vcc2027.log 2>&1 &
MON=$!
trap 'kill $MON 2>/dev/null' EXIT
run() {   # run PART THREADS [soft]
  local part=$1 th=$2 soft=${3:-}
  echo "$(date -u +%T) START $part" >> $L/chain_vcc2027.log
  OMP_NUM_THREADS=$th MKL_NUM_THREADS=$th OPENBLAS_NUM_THREADS=$th ~/venv/bin/python code/vcc6/step20.py $part > $L/$part.log 2>&1
  local rc=$?
  if [ $rc -ne 0 ]; then
    if [ -n "$soft" ]; then echo "$(date -u +%T) STOPPED_SUBPART $part rc $rc (download or check failed; the chain continues without it)" >> $L/chain_vcc2027.log; return 0; fi
    echo "$(date -u +%T) FAILED $part rc $rc" >> $L/chain_vcc2027.log; exit $rc
  fi
  touch $L/$part.done; echo "$(date -u +%T) DONE $part" >> $L/chain_vcc2027.log
}
for p in testF partF3c; do [ -e $L/$p.done ] || run $p 8; done
[ -e $L/partF1.done ] || run partF1 8
[ -e $L/partF2.done ] || run partF2 8
for p in partF3a partF3b; do [ -e $L/$p.done ] || run $p 8 soft; done
for p in partF6 partF3sum partF3d partF4 partF5time; do
  [ -e $L/$p.done ] || run $p $( [ $p = partF5time ] && echo 4 || echo 8 )
done
if [ ! -e $L/partF5.done ]; then
  echo "$(date -u +%T) START partF5 (2 workers x 4 threads)" >> $L/chain_vcc2027.log
  for w in 0 1; do
    OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 POLARS_MAX_THREADS=4 ~/venv/bin/python code/vcc6/step20.py partF5 --worker $w > $L/partF5_w$w.log 2>&1 &
    eval "P$w=$!"
  done
  wait $P0; R0=$?; wait $P1; R1=$?
  if [ $R0 -ne 0 ] || [ $R1 -ne 0 ]; then echo "$(date -u +%T) FAILED partF5 rc $R0 $R1" >> $L/chain_vcc2027.log; exit 1; fi
  touch $L/partF5.done; echo "$(date -u +%T) DONE partF5" >> $L/chain_vcc2027.log
fi
[ -e $L/partF5sum.done ] || run partF5sum 8
echo "$(date -u +%T) CHAIN END" >> $L/chain_vcc2027.log
