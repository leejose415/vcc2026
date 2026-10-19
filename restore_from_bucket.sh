#!/usr/bin/env bash
# Restore archived files from gs://vcc2026 to their original local paths.
#   bash restore_from_bucket.sh <local_path or folder>   (absolute, or relative to ~/vcc)
# Looks the path up in ARCHIVE_MANIFEST.csv (a file, or every file under a folder), copies the exact recorded object
# generation back to the same local path, and verifies size, crc32c and MD5 against the manifest. Files already present
# are only verified. Stops on the first mismatch (the bad copy is renamed to <path>.restore_failed).
set -euo pipefail
V="$HOME/vcc"
MANIFEST="$V/ARCHIVE_MANIFEST.csv"
[ $# -eq 1 ] || { echo "usage: bash restore_from_bucket.sh <local_path or folder>" >&2; exit 2; }
[ -f "$MANIFEST" ] || { echo "no manifest at $MANIFEST" >&2; exit 2; }
target="$1"; [[ "$target" = /* ]] || target="$V/$target"; target="${target%/}"

rows=$(python3 -I - "$MANIFEST" "$target" <<'PY'
import csv, sys
manifest, target = sys.argv[1], sys.argv[2]
seen = {}
for r in csv.DictReader(open(manifest)):
    p = r['local_path']
    if p == target or p.startswith(target + '/'):
        seen[p] = r                      # the latest row for a path wins
for r in seen.values():
    print('\t'.join([r['local_path'], r['size_bytes'], r['crc32c'], r['md5'], r['gs_uri'], r['generation']]))
PY
)
[ -n "$rows" ] || { echo "nothing in the manifest under $target" >&2; exit 1; }

n=0
while IFS=$'\t' read -r path size crc md5 uri gen; do
  if [ ! -e "$path" ]; then
    mkdir -p "$(dirname "$path")"
    gcloud storage cp "$uri#$gen" "$path" --quiet
  fi
  got_size=$(stat -c %s "$path")
  hex=$(gcloud storage hash --hex "$path" --format='value(crc32c_hash,md5_hash)')
  got_crc=$(echo "$hex" | awk '{print $1}'); got_md5=$(echo "$hex" | awk '{print $2}')
  want_crc=$(echo "$crc" | base64 -d | od -An -tx1 | tr -d ' \n')
  if [ "$got_size" != "$size" ] || [ "$got_crc" != "$want_crc" ] || [ "$got_md5" != "$md5" ]; then
    mv "$path" "$path.restore_failed"
    echo "MISMATCH $path: size $got_size/$size crc32c $got_crc/$want_crc md5 $got_md5/$md5" >&2; exit 1
  fi
  n=$((n + 1)); echo "ok $path"
done <<< "$rows"
echo "restored and verified $n file(s)"
