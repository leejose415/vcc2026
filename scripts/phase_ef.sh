#!/usr/bin/env bash
# Phases E-F on vcc2027, unattended: wait for the PIE chain, then
#   dev_pie (validation, PIE on) and final_nopie (final mode, PIE off) in parallel -> final_pie (final, PIE on),
#   effect-space + H1 pipeline for dev_pie, A-C prediction folders from both final runs, PIE-only fallback folders.
# Choosing the model / tau and building the submission are separate, manual steps. Logs: ~/vcc/logs/vcc_2026_4/.
set -eo pipefail
cd "$HOME/vcc"; PY="$HOME/venv/bin/python"; L=logs/vcc_2026_4; R=lib/vcc_2026_4/runs
export PYTHONPATH=code
PIELOG=$L/${PIE_CHAIN:-pie_c3}.log
until grep -q '^EXIT=' $PIELOG 2>/dev/null; do sleep 120; done
grep -q '^EXIT=0' $PIELOG || { echo "PIE chain failed; see $PIELOG"; exit 1; }
for l in K562 RPE1 HepG2 Jurkat A B C; do [ -f lib/pie_tables/pie_$l.h5ad ] || { echo "missing pie_$l"; exit 1; }; done
date -u; echo "== training dev_pie and final_nopie in parallel"
[ -f $R/dev_pie/done.json ] || "$PY" -m transformer.train train --config code/transformer/configs/dev.yaml --name dev_pie --pie on > $L/train_dev_pie.log 2>&1 &
p1=$!
[ -f $R/final_nopie/done.json ] || "$PY" -m transformer.train train --config code/transformer/configs/final.yaml --name final_nopie --pie off > $L/train_final_nopie.log 2>&1 &
p2=$!
wait $p1; wait $p2
date -u; echo "== final_pie training, with dev_pie evaluation alongside"
[ -f $R/final_pie/done.json ] || "$PY" -m transformer.train train --config code/transformer/configs/final.yaml --name final_pie --pie on > $L/train_final_pie.log 2>&1 &
p3=$!
"$PY" -m transformer.evaluate effect --run $R/dev_pie > $L/eval_dev_pie.log 2>&1
"$PY" -m transformer.evaluate predict-lines --run $R/final_nopie --out $R/final_nopie/pred_lines > $L/pred_final_nopie.log 2>&1
"$PY" -m transformer.evaluate predict-lines --pie-only --config code/transformer/configs/final.yaml --out lib/vcc_2026_4/pie_only/pred_lines > $L/pred_pie_only.log 2>&1
"$PY" -m transformer.evaluate h1-pipeline --run $R/dev_pie --taus 0.5 0.7 0.9 --topn 1000 2000 --floor --tag _floor --no-baseline > $L/h1_dev_pie.log 2>&1
wait $p3
"$PY" -m transformer.evaluate predict-lines --run $R/final_pie --out $R/final_pie/pred_lines > $L/pred_final_pie.log 2>&1
date -u; echo "== phase E-F runs done"
