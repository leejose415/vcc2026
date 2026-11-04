#!/bin/bash
# Replogle et al. 2022 K562 genome-wide Perturb-seq, RAW single-cell counts (not the essential set,
# not the normalized matrix: the Replogle loader only accepts raw counts).
# Source: https://plus.figshare.com/articles/dataset/20029387  (file 35775507, 65.83 GB)
#
# 1. Streams figshare -> gs://vcc2026/replogle/ (uses no local disk), verifies the figshare MD5.
# 2. Copies the bucket object to ~/vcc/lib/replogle/ only if enough disk would remain afterwards.
# Safe to re-run: a verified bucket copy or complete local copy is skipped.
#
# Run in tmux (takes ~20-60 min):  tmux new -d -s replogle 'bash ~/vcc/scripts/download_replogle_k562_gwps.sh'
# Follow with:                      tail -f ~/vcc/logs/replogle_download.log
mkdir -p "$HOME/vcc/logs"
exec > >(tee -a "$HOME/vcc/logs/replogle_download.log") 2>&1
trap 'echo "EXIT=$?"' EXIT
set -euo pipefail

NAME=K562_gwps_raw_singlecell_01.h5ad
URL=https://ndownloader.figshare.com/files/35775507
MD5_HEX=887e3e6a8c8df6eadf7a3030a53c9546          # figshare's computed_md5 for this file
GCS=gs://vcc2026/replogle/$NAME
LOCAL_DIR=$HOME/vcc/lib/replogle
MARGIN_GB=${MARGIN_GB:-15}                          # disk that must stay free after the local copy

SIZE=$(curl -fsS https://api.figshare.com/v2/articles/20029387/files \
       | python3 -I -c 'import json,sys; print(next(f["size"] for f in json.load(sys.stdin) if f["id"]==35775507))')
MD5_B64=$(python3 -I -c 'import base64,sys; print(base64.b64encode(bytes.fromhex(sys.argv[1])).decode())' "$MD5_HEX")
echo "$(date -u +%FT%TZ) $NAME: $((SIZE / 1000000000)) GB, md5 $MD5_HEX"

gcs_md5() { gcloud storage objects describe "$GCS" --format='value(md5_hash)' 2>/dev/null || true; }

# ---- 1. figshare -> bucket ----
if [ "$(gcs_md5)" = "$MD5_B64" ]; then
  echo "bucket copy already verified: $GCS"
else
  for attempt in 1 2 3; do
    echo "$(date -u +%FT%TZ) streaming to $GCS (attempt $attempt)"
    if curl -fL --no-progress-meter "$URL" | gcloud storage cp - "$GCS" && [ "$(gcs_md5)" = "$MD5_B64" ]; then
      break
    fi
    echo "attempt $attempt failed or MD5 mismatch (bucket md5: $(gcs_md5))"
    [ "$attempt" = 3 ] && { echo "giving up"; exit 1; }
    sleep 60
  done
  echo "$(date -u +%FT%TZ) bucket copy verified: $GCS"
fi

# ---- 2. bucket -> VM (only if it fits with margin) ----
mkdir -p "$LOCAL_DIR"
if [ -f "$LOCAL_DIR/$NAME" ] && [ "$(stat -c %s "$LOCAL_DIR/$NAME")" = "$SIZE" ]; then
  echo "local copy already present: $LOCAL_DIR/$NAME"
  exit 0
fi
AVAIL=$(df -B1 --output=avail "$LOCAL_DIR" | tail -1)
NEED=$((SIZE + MARGIN_GB * 1000000000))
if [ "$AVAIL" -lt "$NEED" ]; then
  echo "NOT copying to the VM: $((AVAIL / 1000000000)) GB free, need $((NEED / 1000000000)) GB" \
       "(file + ${MARGIN_GB} GB margin). Free disk or grow the boot disk, then re-run this script."
  exit 0
fi
echo "$(date -u +%FT%TZ) copying to $LOCAL_DIR"
gcloud storage cp "$GCS" "$LOCAL_DIR/$NAME"         # gcloud verifies the hash on download
echo "$(date -u +%FT%TZ) done: $LOCAL_DIR/$NAME"
