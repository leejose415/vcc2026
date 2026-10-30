#!/bin/bash
# 15b.4 + addendum O3, part $1 of 8; starts after every oracle part has exited
cd ~/vcc; L=results/stage1/2026-10-29_step15/logs; p=$1
until [ "$(grep -l "^exit" $L/oracle_p*.log 2>/dev/null | wc -l)" = "6" ]; do sleep 20; done
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 POLARS_MAX_THREADS=1
~/venv/bin/python code/vcc6/step15.py predrun --pred b --part $p --nparts 8 > $L/predrun_b_p$p.log 2>&1; echo exit $? >> $L/predrun_b_p$p.log
~/venv/bin/python code/vcc6/step15.py o3 --part $p --nparts 8 > $L/o3_p$p.log 2>&1; echo exit $? >> $L/o3_p$p.log
~/venv/bin/python code/vcc6/step15.py o3 --pred b --part $p --nparts 8 > $L/o3_b_p$p.log 2>&1; echo exit $? >> $L/o3_b_p$p.log
