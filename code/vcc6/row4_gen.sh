#!/bin/bash
# Row 4 D/E/F on vcc2027: inputs (with checks) -> the frozen generator per context (run_line predict --lfc-npz, salt "",
# 400 cells, 8 workers x 1 thread), in sequence. Fixed log names. Stops at the first failure. No upload.
set -uo pipefail
V=$HOME/vcc; P=$HOME/venv/bin/python; R=$V/results/stage1/2026-11-03_row4_DEF; L=$R/logs; mkdir -p $L
CL=$L/gen_chain_vcc2027.log
say() { echo "[$(date -u +%H:%M:%S)] $*" >> $CL; }
( while true; do echo "$(date -u +%T) $(free -m | awk '/Mem:/ {print "used " $3 " MB, available " $7 " MB"} /Swap:/ {print "swap used " $3 " MB"}' | tr '\n' ' ') disk_free $(df -BG / | awk 'NR==2{print $4}')"; sleep 60; done ) >> $L/mem_monitor_vcc2027.log 2>&1 &
MON=$!; trap 'kill $MON 2>/dev/null' EXIT
cd $V
say "start inputs (estimate 3 min)"
env OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 $P code/vcc6/row4_def.py inputs > $L/inputs.log 2>&1 || { say "FAILED inputs"; exit 1; }
say "end inputs"
for x in D E F; do
  [ -f $R/$x/e2e/salt_none/run.json ] && { say "$x already generated"; continue; }
  say "start generate_$x (estimate 22 min; disk free $(df -BG / | awk 'NR==2{print $4}'))"
  env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 $P code/vcc6/run_line.py predict \
    --controls $V/data_final/extracted/context_$x.h5ad --knockdowns $V/results/target_lists/final_DEF_300.txt \
    --background $x --lfc-npz $R/$x/input.npz --salt "" --workers 8 --out $R/$x/e2e > $L/generate_$x.log 2>&1
  rc=$?; [ $rc -eq 0 ] || { say "FAILED generate_$x rc $rc"; exit 1; }
  say "end generate_$x"
done
say GEN_DONE
