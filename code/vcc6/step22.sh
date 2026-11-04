#!/bin/bash
# Step 22 chain on vcc2027 (Jose 2026-11-03): runs the named parts in order, each in its own process, stopping at the
# first failure. One chain at a time on the VM, so no table load runs next to scoring. Logs a runtime estimate per part.
# usage: bash code/vcc6/step22.sh TAG part [part ...]      (TAG names this chain's log: chain_<TAG>.log)
set -u
cd ~/vcc
TAG=$1; shift
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 POLARS_MAX_THREADS=8
L=results/stage1/2026-11-03_step22/logs; mkdir -p $L
declare -A EST=([test]=1 [part0]=8 [srcdelta]=10 [P0]=40 [P1]=15 [P2]=15 [S1]=20 [S2]=10 [S3]=20 [S4]=15 [Q1]=10 [Q2]=12
                [R2]=8 [R3ranks]=8 [R3gen]=12 [R3score]=25 [R3cmp]=6 [D]=10 [P2b]=10)
( while true; do echo "$(date -u +%T) $(free -m | awk '/Mem:/ {print "used " $3 " MB, available " $7 " MB"} /Swap:/ {print "swap used " $3 " MB"}' | tr '\n' ' ') disk_free $(df -BG / | awk 'NR==2{print $4}') $(ps -eo rss,args --sort=-rss | awk 'NR==2 {printf "top rss %.1f GB %s %s", $1/1e6, $3, $4}')"; sleep 60; done ) >> $L/mem_monitor_vcc2027.log 2>&1 &
MON=$!
trap 'kill $MON 2>/dev/null' EXIT
for p in "$@"; do
  [ -e $L/$p.done ] && { echo "$(date -u +%T) SKIP $p (done)" >> $L/chain_$TAG.log; continue; }
  echo "$(date -u +%T) START $p (estimate ${EST[$p]:-?} min)" >> $L/chain_$TAG.log
  if [ $p = R3score ]; then
    SW=$(swapon --show --bytes --noheadings | awk '{s+=$3} END {print s/1e9}')
    echo "$(date -u +%T) R3score: estimated peak ~30 GB (Eval A harness scoring); swap ${SW} GB (the VM's own /swapfile); available $(free -g | awk '/Mem:/ {print $7}') GB; disk free $(df -BG / | awk 'NR==2{print $4}')" >> $L/chain_$TAG.log
    awk -v s="$SW" 'BEGIN {exit !(s >= 30)}' || { echo "$(date -u +%T) FAILED R3score swap ${SW} GB < 30 GB" >> $L/chain_$TAG.log; exit 1; }
    ~/venv/bin/python code/vcc6/run_line.py score --pred results/stage1/2026-11-03_step22/R3/EvalA_hybrid/salt_none/prediction.h5ad \
      --reference lib/harness/evalA/real.h5ad --threads 8 > $L/$p.log 2>&1
  else
    ~/venv/bin/python code/vcc6/step22.py $p > $L/$p.log 2>&1
  fi
  rc=$?
  if [ $rc -ne 0 ]; then echo "$(date -u +%T) FAILED $p rc $rc" >> $L/chain_$TAG.log; exit $rc; fi
  touch $L/$p.done; echo "$(date -u +%T) DONE $p" >> $L/chain_$TAG.log
done
echo "$(date -u +%T) CHAIN END $TAG" >> $L/chain_$TAG.log
