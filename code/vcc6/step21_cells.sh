#!/bin/bash
# Step 21 Track 1 chain on vcc2027 (Jose 2026-11-02): row 4 on cells. Strictly in sequence:
# test -> f5screen -> t1inputs (checks: F5 + shift removed = the spec's input; calls = hybrid T 5 bit for bit) ->
# generate Eval A, HCT116final (frozen generator, salt "", 400 cells) -> swap check -> score both (frozen harness) ->
# boot_score (2,000 resamples) against the live entry's 15c cells. Stops at the first failure. Fixed log names.
set -uo pipefail
V=$HOME/vcc; P=$HOME/venv/bin/python; C=$V/code/vcc6
R=$V/results/stage1/2026-11-02_step21; L=$R/logs; CE=$R/cells; mkdir -p $L
FE=$V/results/stage1/final_evals/harness/lib; C15=$V/results/stage1/2026-10-29_step15/c/config1p_K0
declare -A CTRL=([EvalA]=$V/lib/harness/evalA/context.h5ad [HCT116final]=$FE/evalHCT116final/context.h5ad)
declare -A KD=([EvalA]=$V/results/target_lists/evalA_h1_292.txt [HCT116final]=$V/results/target_lists/final_HCT116.txt)
declare -A BG=([EvalA]=H1 [HCT116final]=HCT116)
declare -A REF=([EvalA]=$V/lib/harness/evalA/real.h5ad [HCT116final]=$FE/evalHCT116final/real.h5ad)
declare -A HN=([EvalA]=evalA [HCT116final]=evalHCT116final)
CL=$L/chain_vcc2027.log
say() { echo "[$(date -u +%H:%M:%S)] $*" >> $CL; }
( while true; do echo "$(date -u +%T) $(free -m | awk '/Mem:/ {print "used " $3 " MB, available " $7 " MB"} /Swap:/ {print "swap used " $3 " MB"}' | tr '\n' ' ') disk_free $(df -BG / | awk 'NR==2{print $4}') $(ps -eo rss,args --sort=-rss | awk 'NR==2 {printf "top rss %.1f GB %s %s", $1/1e6, $3, $4}')"; sleep 60; done ) >> $L/mem_monitor_vcc2027.log 2>&1 &
MON=$!; trap 'kill $MON 2>/dev/null' EXIT
stage() { local name=$1; shift; [ -e $L/$name.done ] && { say "skip $name (done)"; return 0; }
  say "START $name"; "$@" > $L/$name.log 2>&1; local rc=$?
  if [ $rc -ne 0 ]; then say "FAILED $name rc $rc"; exit 1; fi; touch $L/$name.done; say "DONE $name"; }
T8="env OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 POLARS_MAX_THREADS=8"
T1="env OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 POLARS_MAX_THREADS=1"
cd $V
stage t1_test $T8 $P $C/step21.py test
stage t1_f5screen $T8 $P $C/step21.py f5screen
stage t1_inputs $T8 $P $C/step21.py t1inputs
for l in EvalA HCT116final; do
  stage generate_$l $T1 $P $C/run_line.py predict --controls ${CTRL[$l]} --knockdowns ${KD[$l]} --background ${BG[$l]} \
    --lfc-npz $CE/$l/input.npz --salt "" --workers 8 --out $CE/$l/e2e
done
SW=$(swapon --show --bytes --noheadings | awk '{s+=$3} END {print s/1e9}')
say "swap before scoring: ${SW} GB (the VM's own /swapfile; a second 32 GB swapfile does not fit the 32 GB free disk); disk free $(df -BG / | awk 'NR==2{print $4}')"
awk -v s="$SW" 'BEGIN {exit !(s >= 30)}' || { say "FAILED swap check: ${SW} GB < 30 GB"; exit 1; }
for l in EvalA HCT116final; do
  stage score_$l $T8 $P $C/run_line.py score --pred $CE/$l/e2e/salt_none/prediction.h5ad --reference ${REF[$l]} --threads 8
done
for l in EvalA HCT116final; do
  stage boot_$l $T8 $P $C/boot_score.py $C15/$l/e2e/salt_none/harness/${HN[$l]} $CE/$l/e2e/salt_none/harness/${HN[$l]} --out $CE/$l/boot_vs_live.json
done
say CHAIN_DONE
