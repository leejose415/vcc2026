#!/bin/bash
# 2026-10-28 D/E/F, after the 06:24 OOM (model and generation in one process with 117 tables loaded): model step
# per context in its own process (model2_port.py cg = run_line's model step, exits), then generation from the saved
# arrays (run_line predict --lfc-npz, no tables). D's model (05:35) is reused. Every exit code checked.
V=$HOME/vcc; P=$HOME/venv/bin/python; D=$V/results/stage1/2026-10-28_DEF; KD=$V/results/target_lists/final_DEF_300.txt
say() { echo "[$(date -u +%H:%M:%S)] $*" | tee -a $D/logs/run_def.log; }
run() { local n=$1; shift; say "start $n"; "$@" > $D/logs/$n.log 2>&1; local rc=$?; say "end $n: exit $rc"; [ $rc -eq 0 ] || { say STOPPED; exit $rc; }; }
cd $V
run input_D $P -c "
import numpy as np, sys; sys.path.insert(0,'code/vcc6'); import model2_port as M
z=np.load('$D/D/model/predictions.npz', allow_pickle=True); want=[k for k in open('$KD').read().split() if k]
assert list(map(str,z['knockdowns']))==want
M.write_input('$D/D/input.npz', want, [str(g) for g in z['genes']], z['lfc'], z['de_sig'], z['n_sources']); print('ok')"
for X in E F; do
  run model_$X $P code/vcc6/model2_port.py cg --controls $V/data_final/extracted/context_$X.h5ad --knockdowns $KD --background $X --out $D/${X}_model
  mkdir -p $D/$X && cp $D/${X}_model/CG/input.npz $D/$X/input.npz
done
for X in D E F; do
  run gen_$X $P code/vcc6/run_line.py predict --controls $V/data_final/extracted/context_$X.h5ad --knockdowns $KD --background $X \
    --lfc-npz $D/$X/input.npz --salt "" --workers 4 --out $D/$X/e2e
done
say PREDICT_DONE
