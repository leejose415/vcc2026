#!/usr/bin/env bash
# pie prep (controls only) + pie infer for one line, timed. Inputs from scripts/pie_contexts.py.
# Usage: bash scripts/pie_run.sh <line>   (K562, RPE1, HepG2, Jurkat, A, B, C)
# Output: lib/pie/predictions/<line>.parquet and lib/pie/predictions/<line>_timing.json
set -euo pipefail
line="$1"; low=$(echo "$line" | tr '[:upper:]' '[:lower:]')
case "$line" in HCT116|HEK293T|H1) echo "PIE is leaked for $line; refusing"; exit 1;; esac
V="$HOME/vcc"; PIE="$HOME/pie"
export PATH="$HOME/.local/bin:$PATH"; source "$HOME/.config/pie/common.sh"
prep="$V/lib/pie/prep/$line"; out="$V/lib/pie/predictions/$line.parquet"; mkdir -p "$V/lib/pie/predictions"
cd "$PIE"
t0=$(date +%s)
uv run pie prep name="vcc_$low" h5ad="$V/lib/pie/contexts/$line.h5ad" output_dir="$prep" controls_only=true \
  obs.context_col=context obs.pert_col=target_gene obs.control_label=non-targeting overwrite=true
t1=$(date +%s)
uv run pie infer experiment_name=replogle_xdataset rows_kind=query rows_path="$V/lib/pie/queries/$line.json" \
  "preprocessed_dirs=[$prep]" output_path="$out" overwrite=true device=cuda batch_size=16
t2=$(date +%s)
rows=$(python3 -c "import json; print(sum(len(v) for v in json.load(open('$V/lib/pie/queries/$line.json')).values()))")
echo "{\"line\": \"$line\", \"rows\": $rows, \"prep_seconds\": $((t1 - t0)), \"infer_seconds\": $((t2 - t1))}" | tee "$V/lib/pie/predictions/${line}_timing.json"
