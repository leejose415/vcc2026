#!/bin/bash
# Step 20 chain on vcc2026 (Jose 2026-11-02): test, Part 0, A, B (pools, then scoring), C, D; then G once Part F's
# outputs are copied back (logs/partF_copied.ok) and partG is ready (logs/partG_ready.ok). Stops at the first failure.
set -u
cd ~/vcc
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 POLARS_MAX_THREADS=4
L=results/stage1/2026-11-02_step20/logs
( while true; do echo "$(date -u +%T) $(free -m | awk '/Mem:/ {print "used " $3 " MB, available " $7 " MB"}') $(ps -eo rss,args --sort=-rss | awk 'NR==2 {printf "top rss %.1f GB %s %s", $1/1e6, $3, $4}')"; sleep 60; done ) >> $L/mem_monitor_vcc2026.log 2>&1 &
MON=$!
trap 'kill $MON 2>/dev/null' EXIT
run() {
  local part=$1; shift
  echo "$(date -u +%T) START $part" >> $L/chain_vcc2026.log
  ~/venv/bin/python code/vcc6/step20.py "$@" > $L/$part.log 2>&1
  local rc=$?
  if [ $rc -ne 0 ]; then echo "$(date -u +%T) FAILED $part rc $rc" >> $L/chain_vcc2026.log; exit $rc; fi
  touch $L/$part.done; echo "$(date -u +%T) DONE $part" >> $L/chain_vcc2026.log
}
for p in test part0 partA partB partBscore partC partD; do
  [ -e $L/$p.done ] || run $p $p
done
echo "$(date -u +%T) WAITING for Part F outputs and partG" >> $L/chain_vcc2026.log
while [ ! -e $L/partF_copied.ok ] || [ ! -e $L/partG_ready.ok ]; do sleep 60; done
[ -e $L/partG.done ] || run partG partG
echo "$(date -u +%T) CHAIN END" >> $L/chain_vcc2026.log
