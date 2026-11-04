#!/bin/bash
# Step 18 generation on one VM: for each context given, the frozen generator from the saved arrays (run_line predict
# --lfc-npz), in sequence. Usage: bash code/vcc6/step18_gen.sh WORKERS D [E F]
set -uo pipefail
W=$1; shift
V=$HOME/vcc; P=$HOME/venv/bin/python; R=$V/results/stage1/2026-10-30_step18; L=$R/logs
mkdir -p $L
say() { echo "[$(date -u +%H:%M:%S)] $*" | tee -a $L/gen_chain_$(hostname).log; }
for x in "$@"; do
  [ -f $R/config1p_K0/$x/e2e/salt_none/run.json ] && { say "$x already generated"; continue; }
  say "start generate_$x"
  env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 $P $V/code/vcc6/run_line.py predict \
    --controls $V/data_final/extracted/context_$x.h5ad --knockdowns $V/results/target_lists/final_DEF_300.txt \
    --background $x --lfc-npz $R/config1p_K0/$x/input.npz --salt "" --workers $W --out $R/config1p_K0/$x/e2e \
    > $L/generate_$x.log 2>&1
  rc=$?; say "end generate_$x: exit $rc"; [ $rc -eq 0 ] || { say STOPPED; exit 1; }
done
say GEN_DONE
