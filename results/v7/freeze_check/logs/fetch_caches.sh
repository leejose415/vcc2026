#!/bin/bash
# Retry starting vcc2027 (zone STOCKOUT), copy the H1 proxy caches, stop vcc2027. Jose approved 2026-10-29.
Z=us-central1-c; R=vcc/results/stage1/2026-10-28_scorer_space/cache; D=$HOME/vcc/results/stage1/2026-10-28_scorer_space/cache
mkdir -p $D
for i in $(seq 1 30); do
  if gcloud compute instances start vcc2027 --zone $Z >/dev/null 2>&1; then echo "$(date -u +%T) started on attempt $i"; break; fi
  echo "$(date -u +%T) attempt $i: not started"; timeout 120 tail -f /dev/null
done
[ "$(gcloud compute instances describe vcc2027 --zone $Z --format='value(status)')" = RUNNING ] || { echo "gave up"; exit 2; }
for j in $(seq 1 20); do gcloud compute ssh vcc2027 --zone $Z --command true >/dev/null 2>&1 && break; timeout 15 tail -f /dev/null; done
gcloud compute ssh vcc2027 --zone $Z --command "cd ~/$R && ls -la fold_H1.npz s11_H1.npz && md5sum fold_H1.npz s11_H1.npz" > $D/remote_md5.txt 2>&1
gcloud compute scp --zone $Z "vcc2027:~/$R/fold_H1.npz" "vcc2027:~/$R/s11_H1.npz" $D/; rc=$?
cat $D/remote_md5.txt; (cd $D && md5sum fold_H1.npz s11_H1.npz)
gcloud compute instances stop vcc2027 --zone $Z >/dev/null 2>&1
echo "$(date -u +%T) copy rc $rc; vcc2027 $(gcloud compute instances describe vcc2027 --zone $Z --format='value(status)')"
exit $rc
