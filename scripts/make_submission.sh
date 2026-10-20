#!/usr/bin/env bash
# Build a leaderboard file from per-line prediction folders: existing generator (vcc3_generate.py submit, 400 fresh
# control-based cells per knockdown) -> vcc prep. Never submits.
# Usage: bash scripts/make_submission.sh <pred dir with A/ B/ C/> <array, e.g. tf_tau0.5> <generator settings json> <run name>
# Settings json: a selected_generator.json from an H1 benchmark (its settings block is used; the array is given here).
set -euo pipefail
pred=$(realpath "$1"); array="$2"; settings=$(realpath "$3"); name="$4"   # absolute: the generator runs from its own folder
V="$HOME/vcc"; PY="$HOME/venv/bin/python"; out="$V/results/$name"
[ -e "$out/prediction.vcc" ] && { echo "$out already has prediction.vcc: choose a new run name"; exit 1; }
if [ -e "$out/generator_report.json" ]; then
  echo "generator output already complete in $out: packaging only"
else
  [ -e "$out/prediction.h5ad" ] && { echo "$out has an unfinished prediction.h5ad: choose a new run name"; exit 1; }
  for line in A B C; do [ -f "$pred/$line/$array.npy" ] || { echo "missing $pred/$line/$array.npy"; exit 1; }; done
  free=$(df --output=avail -BG "$V" | tail -1 | tr -dc 0-9); [ "$free" -ge 40 ] || { echo "only ${free} GB free; need ~40"; exit 1; }
  cd "$V/code/vcc_2026_3c"
  # GEN_EXTRA: extra generator switches, e.g. "--exact-totals --workers 7"
  "$PY" vcc3_generate.py submit --data "$V/data" --lines A B C --delta-array "$array" --settings "$settings" \
    --delta "A=$pred/A" --delta "B=$pred/B" --delta "C=$pred/C" --out "$out" ${GEN_EXTRA:-}
fi
export TMPDIR="$V/tmp_prep"; mkdir -p "$TMPDIR"     # vcc prep writes a full temporary copy; /tmp is a small tmpfs here
mem=$(awk '/MemAvailable|SwapFree/ {s += $2} END {print int(s / 1e6)}' /proc/meminfo)
[ "$mem" -ge 34 ] || { echo "only ${mem} GB RAM + swap free; vcc prep peaks near 30 GB"; exit 1; }
"$HOME/venv/bin/vcc" prep "$out/prediction.h5ad" -g "$V/data/gene_names.csv" --perts "$V/data/pert_counts.csv" -o "$out/prediction.vcc"
rm -rf "$TMPDIR"
ls -la "$out"
echo "READY (not submitted): $out/prediction.vcc"
