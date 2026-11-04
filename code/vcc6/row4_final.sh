#!/bin/bash
# Row 4 D/E/F final on vcc2026 (no upload): fetch the three generator outputs from vcc2027 (md5 both sides) ->
# merge straight from them (row4_def.py final) -> check before prep (row4_def.py check) -> delete this VM's COPIES of
# the three generator outputs (originals stay on vcc2027) -> disk gates -> temporary 18 GB swap -> vcc prep (exit
# checked; memory / disk logged) -> swap removed -> .vcc payload check (row4_def.py check_vcc). Stops at the first failure.
# usage: bash code/vcc6/row4_final.sh ATTEMPT
set -uo pipefail
A=${1:-1}
V=$HOME/vcc; P=$HOME/venv/bin/python; R=results/stage1/2026-11-03_row4_DEF; L=$V/$R/logs; F=$V/$R/final; SW=/swapfile_v7
CL=$L/final_chain_attempt$A.log; mkdir -p $L $F
SSH="gcloud compute ssh vcc2027 --zone us-central1-c --internal-ip --quiet --command"
say() { echo "[$(date -u +%H:%M:%S)] $*" >> $CL; }
fail() { say "FAILED: $*"; [ -f $SW ] && { sudo swapoff $SW 2>/dev/null; sudo rm -f $SW; say "swap removed"; }; say STOPPED; exit 1; }
cd $V
# 1. fetch
if [ ! -e $L/fetch.done ]; then
  say "fetch from vcc2027 (estimate 3 min)"
  $SSH "cd ~/vcc && find $R/D $R/E $R/F $R/tables $R/logs -type f | sort > /tmp/r4.txt && xargs md5sum < /tmp/r4.txt > /tmp/r4.md5 && tar cf - -T /tmp/r4.txt" | tar xf - -C $V || fail fetch
  $SSH "cat /tmp/r4.md5" > $L/fetch_md5_vcc2027.txt
  md5sum -c --quiet $L/fetch_md5_vcc2027.txt > $L/fetch_check.txt 2>&1 || { grep -v 'logs/' $L/fetch_check.txt | grep -q FAILED && fail "fetch md5"; }
  touch $L/fetch.done; say "fetch md5 OK ($(wc -l < $L/fetch_md5_vcc2027.txt) files); disk free $(df -BG / | awk 'NR==2{print $4}')"
fi
# 2. merge + 3. checks before prep
[ -e $L/final.done ] || { say "start final (estimate 13 min)"; env OMP_NUM_THREADS=4 $P code/vcc6/row4_def.py final > $L/final_attempt$A.log 2>&1 || fail final; touch $L/final.done; say "end final"; }
[ -e $L/check_gen.done ] || { say "start check (estimate 10 min)"; env OMP_NUM_THREADS=4 $P code/vcc6/row4_def.py check > $L/check_gen_attempt$A.log 2>&1 || fail check; touch $L/check_gen.done; say "end check: PASS"; }
# 4. delete this VM's copies of the generator outputs (contained in the merged file; originals on vcc2027)
for x in D E F; do
  f=$R/$x/e2e/salt_none/prediction.h5ad
  if [ -f $f ]; then say "delete local copy $f ($(du -h $f | cut -f1); md5 $(grep " $f\$" $L/fetch_md5_vcc2027.txt | cut -d' ' -f1); original kept on vcc2027)"; rm -f $f; fi
done
# 5. disk gates and swap
FREE=$(df -BG / | awk 'NR==2{print $4}' | tr -d G); say "disk free before swap: ${FREE} GB"
[ "$FREE" -ge 42 ] || fail "disk free ${FREE} GB < 42 GB before swap"
if [ ! -f $SW ]; then { sudo fallocate -l 18G $SW && sudo chmod 600 $SW && sudo mkswap $SW >/dev/null && sudo swapon $SW; } || fail "swap"; fi
FREE=$(df -BG / | awk 'NR==2{print $4}' | tr -d G); say "swap on: $(free -g | awk '/Swap/{print $2}') GB; disk free after swap: ${FREE} GB"
[ "$FREE" -ge 22 ] || fail "disk free ${FREE} GB < 22 GB after swap"
# 6. vcc prep
say "start vcc_prep (estimate 17 min; peak ~32 GB RAM + ~12 GB swap, ~19 GB temp disk)"
( while sleep 30; do echo "$(date -u +%T) $(free -m | awk '/Mem/{print "mem_used_MB",$3} /Swap/{print "swap_used_MB",$3}' | tr "\n" " ") disk_free $(df -BG / | awk 'NR==2{print $4}')"; done ) > $L/prep_mem_attempt$A.log & MON=$!
$HOME/venv/bin/vcc prep $F/prediction_DEF.h5ad -g $V/data_final/extracted/gene_names.csv \
  --perts $V/data_final/extracted/pert_counts.csv -o $F/prediction_DEF.vcc --force > $L/vcc_prep_attempt$A.log 2>&1; rc=$?
kill $MON 2>/dev/null; say "end vcc_prep: exit $rc"
sudo swapoff $SW && sudo rm -f $SW && say "swap removed; disk free $(df -BG / | awk 'NR==2{print $4}')"
[ $rc -eq 0 ] || { say "STOPPED at vcc_prep"; exit 1; }
# 7. payload check
say "start check_vcc (estimate 10 min)"
env OMP_NUM_THREADS=4 $P code/vcc6/row4_def.py check_vcc > $L/check_vcc_attempt$A.log 2>&1 || { say "FAILED check_vcc"; say STOPPED; exit 1; }
say "end check_vcc: PASS; .vcc md5 $(md5sum $F/prediction_DEF.vcc | cut -d' ' -f1), $(du -h $F/prediction_DEF.vcc | cut -f1)"
say FINAL_DONE
