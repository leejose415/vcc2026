#!/bin/bash
# 2026-10-28 D/E/F: after F, combine the three contexts, vcc prep (exit checked), vcc submit (Jose: "run that on D/E/F and
# upload that"), then remove the temporary swap. Stops at the first failure; nothing is uploaded unless prep exits 0.
V=$HOME/vcc; P=$HOME/venv/bin/python; D=$V/results/stage1/2026-10-28_DEF
say() { echo "[$(date -u +%H:%M:%S)] $*" | tee -a $D/logs/run_def.log; }
until tail -1 $D/logs/run_def.log | grep -qE "PREDICT_DONE|STOPPED"; do sleep 20; done
tail -1 $D/logs/run_def.log | grep -q PREDICT_DONE || { say "prep_submit: predictions did not finish"; exit 1; }
say "swap: adding 16 GB for vcc prep"
sudo fallocate -l 16G /swapfile_eval && sudo chmod 600 /swapfile_eval && sudo mkswap /swapfile_eval > /dev/null && sudo swapon /swapfile_eval || { say "swap failed"; exit 1; }
say "start combine_prep"
cd $V && $P code/vcc6/def_file.py --runs $D --out $D/final > $D/logs/combine_prep.log 2>&1; rc=$?
say "end combine_prep: exit $rc"
if [ $rc -ne 0 ]; then say STOPPED; sudo swapoff /swapfile_eval && sudo rm /swapfile_eval; exit $rc; fi
say "start vcc_submit"
$HOME/venv/bin/vcc submit $D/final/prediction_DEF.vcc -m LHY > $D/logs/vcc_submit.log 2>&1; rc=$?
say "end vcc_submit: exit $rc"
sudo swapoff /swapfile_eval && sudo rm /swapfile_eval && say "swap removed"
[ $rc -eq 0 ] && say SUBMITTED || say "SUBMIT FAILED"
