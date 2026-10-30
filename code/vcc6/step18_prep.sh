#!/bin/bash
# Step 18, vcc prep retry on the existing merged file (attempt 1 hit ENOSPC: vcc prep writes an uncompressed copy of the
# prediction, ~15.5 GB, plus its zstd, ~3.2 GB, to a temp dir on the same disk). Temporary 18 GB swap (attempt 1 peaked
# at 31.8 GB RAM + 13.9 GB swap), vcc prep with exit check, swap removed. No upload.
set -uo pipefail
V=$HOME/vcc; R=$V/results/stage1/2026-10-30_step18; L=$R/logs; F=$R/final; SW=/swapfile_v7
say() { echo "[$(date -u +%H:%M:%S)] $*" | tee -a $L/final_chain.log; }
{ sudo fallocate -l 18G $SW && sudo chmod 600 $SW && sudo mkswap $SW >/dev/null && sudo swapon $SW; } || { say "swap FAILED"; exit 1; }
say "retry: swap on $(free -g | awk '/Swap/{print $2}') GB; disk free $(df -h / | awk 'NR==2{print $4}')"
( while sleep 30; do echo "$(date -u +%T) $(free -m | awk '/Mem/{print "mem_used_MB",$3} /Swap/{print "swap_used_MB",$3}' | tr "\n" " ") disk_free $(df -BG / | awk 'NR==2{print $4}')"; done ) > $L/prep_mem.log & MON=$!
say "start vcc_prep (retry)"
$HOME/venv/bin/vcc prep $F/prediction_DEF.h5ad -g $V/data_final/extracted/gene_names.csv \
  --perts $V/data_final/extracted/pert_counts.csv -o $F/prediction_DEF.vcc --force > $L/vcc_prep.log 2>&1; rc=$?
kill $MON 2>/dev/null; say "end vcc_prep (retry): exit $rc"
sudo swapoff $SW && sudo rm -f $SW && say "swap removed"
[ $rc -eq 0 ] && say FINAL_DONE || say "STOPPED at vcc_prep (retry)"
