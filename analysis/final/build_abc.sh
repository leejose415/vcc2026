#!/usr/bin/env bash
# Build an A-C leaderboard file (never submits): generator (vcc3_generate.py or gen_f1.py) -> vcc prep -> checks.
# Usage: bash analysis/final/build_abc.sh <pred dir with A/ B/ C/> <array> <settings json> <run name> [generator switches...]
#   GEN_SCRIPT=analysis/final/gen_f1.py for the F1 emission (default: the standard generator)
set -euo pipefail
pred=$(realpath "$1"); array="$2"; settings=$(realpath "$3"); name="$4"; shift 4
V="$HOME/vcc"; PY="$HOME/venv/bin/python"; out="$V/results/$name"; S=$V/scripts/v5
GEN=${GEN_SCRIPT:-$V/code/vcc_2026_3c/vcc3_generate.py}; GEN=$(realpath "$GEN")
[ -e "$out/prediction.vcc" ] && { echo "$out already has prediction.vcc"; exit 1; }
free=$(df --output=avail -BG "$V" | tail -1 | tr -dc 0-9); [ "$free" -ge 38 ] || { echo "only ${free} GB free"; exit 1; }
cd "$V/code/vcc_2026_3c"
date -u; "$PY" "$GEN" submit --data "$V/data" --lines A B C --delta-array "$array" --settings "$settings" \
  --delta "A=$pred/A" --delta "B=$pred/B" --delta "C=$pred/C" --out "$out" "$@"
export TMPDIR="$V/tmp_prep_$name"; mkdir -p "$TMPDIR"
date -u; "$HOME/venv/bin/vcc" prep "$out/prediction.h5ad" -g "$V/data/gene_names.csv" --perts "$V/data/pert_counts.csv" -o "$out/prediction.vcc"
rm -rf "$TMPDIR"
date -u; cd "$V"; PYTHONPATH=$S "$PY" $S/filecheck.py check "$out" lib/vcc_2026_5/old_0638_pbulk.npz results/vcc_2026_4_final_nopie_top5000/prediction.vcc "$pred" "$array" > "$out/checks.log" 2>&1 || echo "checks failed"
tail -40 "$out/checks.log"; ls -la "$out"; echo "READY (not submitted): $out/prediction.vcc"
