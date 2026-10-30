#!/bin/bash
# Step 18 final on vcc2026: temporary 24 GB swap (disk: ~65 GB free, ~27 GB of outputs) -> submission files + on-disk merge + checks (step18.py final) ->
# vcc prep (final gene_names.csv, pert_counts.csv; exit checked) -> swap removed. No upload.
set -uo pipefail
V=$HOME/vcc; P=$HOME/venv/bin/python; R=$V/results/stage1/2026-10-30_step18; L=$R/logs; F=$R/final; SW=/swapfile_v7
say() { echo "[$(date -u +%H:%M:%S)] $*" | tee -a $L/final_chain.log; }
off() { sudo swapoff $SW 2>/dev/null; sudo rm -f $SW; say "swap removed"; }
[ -f $SW ] || { sudo fallocate -l 24G $SW && sudo chmod 600 $SW && sudo mkswap $SW >/dev/null && sudo swapon $SW; } || { say "swap FAILED"; exit 1; }
say "swap on: $(free -g | awk '/Swap/{print $2}') GB"
say "start final"; $P $V/code/vcc6/step18.py final > $L/final.log 2>&1; rc=$?; say "end final: exit $rc"
[ $rc -eq 0 ] || { say STOPPED; off; exit 1; }
say "start vcc_prep"; ( while sleep 30; do echo "$(date -u +%T) $(free -m | awk '/Mem/{print "mem_used_MB",$3} /Swap/{print "swap_used_MB",$3}' | tr "\n" " ")"; done ) > $L/prep_mem.log & MON=$!
$HOME/venv/bin/vcc prep $F/prediction_DEF.h5ad -g $V/data_final/extracted/gene_names.csv \
  --perts $V/data_final/extracted/pert_counts.csv -o $F/prediction_DEF.vcc --force > $L/vcc_prep.log 2>&1; rc=$?
say "end vcc_prep: exit $rc"; kill $MON 2>/dev/null
off
[ $rc -eq 0 ] && say FINAL_DONE || say "STOPPED at vcc_prep"
