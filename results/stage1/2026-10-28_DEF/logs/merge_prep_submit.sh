#!/bin/bash
# 2026-10-28 09:47: the in-memory combine (def_file.py) was stopped at 30 GB RSS + 15 GB swap. The three per-context
# submission files it had already written are merged on disk (anndata concat_on_disk, streaming), then vcc prep (exit
# checked), then vcc submit (Jose's explicit OK), then the temporary swap is removed. Stops at the first failure.
V=$HOME/vcc; P=$HOME/venv/bin/python; D=$V/results/stage1/2026-10-28_DEF; F=$D/final
say() { echo "[$(date -u +%H:%M:%S)] $*" | tee -a $D/logs/run_def.log; }
say "start merge_on_disk"
$P -c "
import anndata.experimental as e, anndata as ad
F='$F'
e.concat_on_disk([f'{F}/submission_{x}.h5ad' for x in 'DEF'], f'{F}/prediction_DEF.h5ad', join='inner', merge='same', keys=list('DEF'), index_unique='-')
a=ad.read_h5ad(f'{F}/prediction_DEF.h5ad', backed='r'); print(a.shape, a.obs.context.value_counts().to_dict(), a.obs_names.is_unique, a.var_names[:3].tolist())
" > $D/logs/merge.log 2>&1; rc=$?; say "end merge_on_disk: exit $rc"; [ $rc -eq 0 ] || { say STOPPED; exit 1; }
say "start vcc_prep"
$HOME/venv/bin/vcc prep $F/prediction_DEF.h5ad -g $V/data_final/extracted/gene_names.csv --perts $V/data_final/extracted/pert_counts.csv -o $F/prediction_DEF.vcc --force > $D/logs/vcc_prep.log 2>&1; rc=$?
say "end vcc_prep: exit $rc"; [ $rc -eq 0 ] || { say STOPPED; exit 1; }
say "start vcc_submit"
$HOME/venv/bin/vcc submit $F/prediction_DEF.vcc -m LHY > $D/logs/vcc_submit.log 2>&1; rc=$?
say "end vcc_submit: exit $rc"
sudo swapoff /swapfile_eval && sudo rm /swapfile_eval && say "swap removed"
[ $rc -eq 0 ] && say SUBMITTED || say "SUBMIT FAILED"
