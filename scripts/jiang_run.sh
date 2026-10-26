#!/bin/bash
# vcc2027: Jiang Seurat -> per-line raw h5ad, one pathway at a time (pull rds, export with R, convert, upload, delete).
# A pathway's complete.txt (export done) becomes converted.txt once its files are uploaded, so the converter only sees
# the pathway in hand; an existing export is not pulled or exported again.
cd ~/vcc; B=gs://vcc2026/vm; W=lib/cells/_jiang_work; P=~/venv/bin/python
mkdir -p $W lib/jiang logs/new_context
for pw in IFNB IFNG INS TGFB TNFA; do
  [ -f $W/$pw/converted.txt ] && { echo "$pw: already converted"; continue; }
  r=lib/jiang/Seurat_object_${pw}_Perturb_seq.rds
  if [ ! -f $W/$pw/complete.txt ]; then
    echo "[$(date -u +%H:%M:%S)] $pw: pull"
    gsutil -q cp $B/$r $r || { echo "PULL_FAILED $pw"; exit 1; }
    echo "[$(date -u +%H:%M:%S)] $pw: export"
    /usr/bin/time -v Rscript scripts/extract_jiang_counts.R $r $W/$pw 2> $W/$pw.time || { tail -20 $W/$pw.time; echo "EXPORT_FAILED $pw"; exit 1; }
    grep "Maximum resident" $W/$pw.time
  fi
  rm -f $r
  echo "[$(date -u +%H:%M:%S)] $pw: convert"
  $P scripts/convert_new_contexts.py jiang || { echo "CONVERT_FAILED $pw"; exit 1; }
  p=$(echo $pw | tr A-Z a-z)
  gsutil -q -m cp -n lib/cells/jiang_${p}_*.h5ad $B/lib/cells/ || { echo "UPLOAD_FAILED $pw"; exit 1; }
  rm -f $W/$pw/x.bin $W/$pw/i.bin $W/$pw/p.bin; mv $W/$pw/complete.txt $W/$pw/converted.txt
  echo "[$(date -u +%H:%M:%S)] $pw: done ($(ls lib/cells/jiang_${p}_*.h5ad | wc -l) files uploaded)"; df -h ~ | tail -1
done
echo JIANG_ALL_DONE
