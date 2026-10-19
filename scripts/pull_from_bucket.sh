#!/usr/bin/env bash
# Copy the vcc_2026_4 inputs from gs://vcc2026 into ~/vcc on a fresh VM, then verify every file's size against its
# bucket object. Bucket copies are kept. Usage: bash scripts/pull_from_bucket.sh
set -euo pipefail
V="$HOME/vcc"
# bucket prefix (or object) -> local folder
PAIRS=(
  "gs://vcc2026/vm/lib/effects_v4|lib/effects_v4"
  "gs://vcc2026/vm/lib/qc_cells|lib/qc_cells"
  "gs://vcc2026/vm/lib/pie_features|lib/pie_features"
  "gs://vcc2026/vm/data|data"
  "gs://vcc2026/vm/results/qc|results/qc"
  "gs://vcc2026/vm/results/analysis|results/analysis"
  "gs://vcc2026/vm/code|code"
  "gs://vcc2026/tcell|lib/tcell"
  "gs://vcc2026/lib/h1/H1_full_2025.h5ad|lib/h1"
)
bad=0
for pair in "${PAIRS[@]}"; do
  src="${pair%%|*}"; dst="$V/${pair##*|}"; mkdir -p "$dst"
  if [[ "$src" == *.h5ad ]]; then
    gcloud storage cp --no-clobber "$src" "$dst/" --quiet
    listing=$(gcloud storage ls -l "$src")
    base="$(dirname "$src")"
  else
    gcloud storage rsync -r "$src" "$dst" --quiet
    listing=$(gcloud storage ls -l -r "$src/**")
    base="$src"
  fi
  n=0
  while read -r size _ url; do
    [[ "$url" == gs://* ]] || continue
    rel="${url#"$base"/}"; f="$dst/$rel"
    got=$(stat -c %s "$f" 2>/dev/null || echo missing)
    if [ "$got" != "$size" ]; then echo "SIZE MISMATCH $f: local $got, bucket $size"; bad=$((bad + 1)); fi
    n=$((n + 1))
  done <<< "$listing"
  echo "ok ${pair##*|}: $n files, sizes checked"
done
[ "$bad" -eq 0 ] && echo "ALL SIZES MATCH" || { echo "$bad size mismatches"; exit 1; }
