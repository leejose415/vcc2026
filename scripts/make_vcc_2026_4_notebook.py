import nbformat as nbf
from pathlib import Path

nb = nbf.v4.new_notebook()
md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
C = []

C.append(md("""# vcc_2026_4: transformer with PIE evidence

Goal: a transformer submission to the 2026 leaderboard before it closes (rules: "until October 22, 2026", two entries per
UTC day; planned cutoff 21 Oct 23:59 UTC). Runs on **vcc2027** (L4 24 GB). The logic lives in `scripts/` and
`code/transformer/`; cells call them through `run` (short steps) and `bg` / `status` / `wait` (long steps in tmux, logs in
`~/vcc/logs/vcc_2026_4/`).

| Phase | What | Where |
|---|---|---|
| 0 | GPU VM, driver, leaderboard rules | vcc2027 |
| A | setup, data copy, provenance checks, data inventory | vcc2027 |
| B | baselines v2 (shared gene set S) | CPU |
| C | PIE: install, control-only contexts, K562 benchmark, PIE tables | GPU |
| D | transformer: loader, model, losses, tests, smoke test | GPU |
| E | dev run, effect-space and H1 pipeline validation | GPU |
| F | final-mode retrain, A-C predictions, generator, `vcc prep` | GPU |

**Fixed decisions (effects_v4):** on-target lfc <= -0.5, MIN_CELLS 30, SE delta method floored at 0.05, unmeasured = NaN;
H1's QC-filtered table is training/validation data, the scoring harness keeps the unfiltered H1 cells. PIE is leaked for
HCT116, HEK293T and H1 (never used there); clean for K562, RPE1, HepG2, Jurkat."""))

C.append(code("""import os, sys, shutil, subprocess, shlex, time, json
from pathlib import Path
import pandas as pd
from IPython.display import display, Image
V = Path.home() / 'vcc'
ENV = dict(CODE=str(V), DATA=str(V / 'data'), LIB=str(V / 'lib'), LOGS=str(V / 'logs' / 'vcc_2026_4'),
           RES=str(V / 'results' / 'vcc_2026_4'), B='gs://vcc2026', PY=sys.executable, PIE=str(Path.home() / 'pie'))
os.environ.update(ENV)
os.environ['PATH'] = str(Path.home() / '.local/bin') + os.pathsep + os.environ['PATH']
for d in ('LOGS', 'RES'):
    os.makedirs(ENV[d], exist_ok=True)
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40)

def run(cmd):
    \"\"\"Short step: run now (from ~/vcc), stream the output here, stop the cell on error.\"\"\"
    p = subprocess.Popen(['bash', '-c', 'set -eo pipefail; cd "$CODE"\\n' + cmd], env=os.environ,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    for line in p.stdout:
        print(line, end='')
    if p.wait():
        raise RuntimeError(f'step failed (exit {p.returncode}) - see the output above')

def bg(name, cmd):
    \"\"\"Long step: run in tmux session `name`; output goes to ~/vcc/logs/vcc_2026_4/<name>.log.\"\"\"
    log, sh = Path(ENV['LOGS']) / f'{name}.log', Path(ENV['LOGS']) / f'{name}.sh'
    exports = ''.join(f'export {k}={shlex.quote(v)}\\n' for k, v in ENV.items())
    sh.write_text(f'#!/bin/bash\\n{exports}export PATH={shlex.quote(os.environ["PATH"])}\\nexec > {shlex.quote(str(log))} 2>&1\\n'
                  f'trap \\'echo "EXIT=$?"\\' EXIT\\nset -eo pipefail\\ncd "$CODE"\\n{cmd}\\n')
    subprocess.run(['tmux', 'kill-session', '-t', name], stderr=subprocess.DEVNULL)
    subprocess.run(['tmux', 'new-session', '-d', '-s', name, 'bash', str(sh)], check=True)
    print(f'started "{name}" in the background. Check with: status("{name}")')

def status(name, n=20):
    alive = subprocess.run(['tmux', 'has-session', '-t', name], stderr=subprocess.DEVNULL).returncode == 0
    log = Path(ENV['LOGS']) / f'{name}.log'
    lines = log.read_text().splitlines() if log.exists() else []
    state = 'RUNNING' if alive else ('FINISHED OK' if lines and lines[-1] == 'EXIT=0' else 'STOPPED - read the log below')
    print(f'== {name}: {state}\\n' + '\\n'.join(lines[-n:]))
    return state

def wait(name, every=60):
    \"\"\"Block until tmux session `name` ends (so the notebook runs top to bottom); raise if it failed.\"\"\"
    while subprocess.run(['tmux', 'has-session', '-t', name], stderr=subprocess.DEVNULL).returncode == 0:
        time.sleep(every)
    if status(name, n=5) != 'FINISHED OK':
        raise RuntimeError(f'{name} failed - see {ENV["LOGS"]}/{name}.log')

print('python:', sys.executable)
print('tmux:', shutil.which('tmux') or 'MISSING')
run('nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader || echo "no GPU here"')"""))

C.append(md("""## Phase 0: vcc2027

- vcc2027: us-central1-c, g2-standard-8 (8 vCPU, 31 GB RAM), 1x NVIDIA L4 (23 GB), 300 GB pd-balanced, Ubuntu 26.04 minimal.
- Driver: Google's `cuda_installer.pyz` (repo mode) installs NVIDIA 615, whose DKMS module does not build on kernel
  7.0.0-1014-gcp (`dma_buf_move_notify` undefined). Replaced by Ubuntu's prebuilt, signed modules:
  `linux-modules-nvidia-595-server-open-gcp` + `nvidia-headless-595-server-open` (driver 595.91.07, CUDA 13.2).
- Leaderboard: entries "until October 22, 2026", at most twice per day, days start at midnight UTC; no clock time given,
  so plan for 21 Oct 23:59 UTC. Final entries: test lines released 22 Oct, due 5 Nov 23:59 UTC."""))
C.append(code("""run('nvidia-smi; swapon --show; free -g; df -h /')"""))

C.append(md("""## Phase A: setup, data, provenance, inventory

- Python 3.12 venv (`~/venv`) with the vcc2026 packages and torch 2.14.1+cu130; flash and memory-efficient SDPA kernels run
  in bf16 on a 5,000-token batch (flash 1.4 ms, math 29 ms).
- Data come from the bucket (`scripts/pull_from_bucket.sh`, every file's size checked; bucket copies kept).
- Provenance (`scripts/check_provenance.py`, refuses to continue on failure): every effects_v4 table has se_floor 0.05 and
  fdr 0.05 in `uns`; it was built from `lib/qc_cells/` (cell sources) or its source table (Orion pseudobulk, Marson CD4 T);
  on-target -0.5 and MIN_CELLS 30 come from the `perturbation_qc.json` next to the knockdown_qc that `uns` names (the
  tables themselves do not record them); rows = kept knockdowns; the old `lib/effects/` is not used by any v4 code."""))
C.append(code("""run('bash scripts/pull_from_bucket.sh 2>&1 | grep -E "^ok|MISMATCH|ALL SIZES"')"""))
C.append(code("""run('"$PY" scripts/check_provenance.py')"""))
C.append(md("""### Data inventory
Roles: **target** = trained and validated on (validation mode); **validation** = H1, a target only in final mode;
**evidence-only** = CD4 T. PIE: **clean** lines get PIE tables; **leaked** lines (in PIE's training/validation) never do."""))
C.append(code("""run('"$PY" scripts/data_inventory.py > /dev/null')
display(pd.read_csv(V / 'results/vcc_2026_4/data_inventory.csv'))
display(pd.read_csv(V / 'results/vcc_2026_4/challenge_neighbours.csv'))"""))

C.append(md("""## Phase B: baselines v2

S = the 6,013 panel genes measured by every effects_v4 source. Per control cell, library_S = counts over S and
CPM_S = counts / library_S x 1e6. Per gene: log1p(mean CPM_S), detection rate, percentile rank among the S genes, in_S,
measured (`lib/effects_v4/*_baseline_v2.csv`; older baselines untouched; lfc / se not recomputed). CD4 T has no cells:
its pseudobulk control CPM is rescaled over S; detection rate NaN with `detection_missing = 1`. Orion uses its per-cell
controls. Fallback if B is late for D: detection rate + percentile rank from the panel-only CPMs, no log CPM."""))
C.append(code("""if not all((V / f'lib/effects_v4/{s}_baseline_v2.csv').exists() for s in
           ['replogle_k562_gw', 'arc_h1', 'orion_hct116', 'marson_cd4_stim48hr', 'challenge_A', 'challenge_C']):
    bg('baseline_v2', '"$PY" scripts/build_baseline_v2.py'); wait('baseline_v2')
b2 = json.loads((V / 'lib/effects_v4/baseline_v2.json').read_text())
print('S genes:', b2.pop('_S_genes'))
display(pd.DataFrame(b2).T[['control_cells', 'library_S_genes', 'median_library_S', 'median_ratio_v2_over_panel',
                            'detection_rate_present', 'measured_genes']])"""))

C.append(md("""## Phase C: PIE

PIE (Arc, `replogle_xdataset` checkpoint, best_auprc) in its own env (`~/pie`, uv, Python 3.12, torch 2.10 cu128), with the
pinned sources downloaded on first use. Control-only contexts (`scripts/pie_contexts.py`): the line's control cells,
log1p CP10k, genes = the panel genes the line measured that PIE's gene_text has under the same symbol (18,262 for A-C).
K562, RPE1, HepG2, Jurkat use PIE's own context names (k562, rpe1, hepg2, jurkat) so its context_text applies; A-C have
no context_text (masked). `scripts/pie_run.sh` = `pie prep` + `pie infer` (timed); `scripts/pie_tables.py` writes
`lib/pie_tables/pie_<line>.h5ad` (layers lfc_pred, p_de; NaN where PIE gave no prediction). PIE is leaked for HCT116,
HEK293T and H1: no PIE context or table is ever made for them (the scripts refuse).

**Result (19-20 Oct):** all seven tables built (K562 / RPE1 / HepG2 / Jurkat: their measured genes with gene_text, e.g.
7,681 for K562; A-C: 18,262 genes). Timing on the L4: the first run (line A, 300 knockdowns) took 49 min, almost all
one-off asset download (~68 GB) and building PIE's evidence cache (22 GB); after that K562 (9,245 knockdowns) 644 s,
RPE1 377 s, Jurkat 383 s, each challenge line ~4-6 min, so all PIE-clean lines + A-C cost ~45 GPU-minutes. Two PIE
processes at once need > 50 GB RAM (an OOM kill ended the first chain): run one at a time. **K562 benchmark: PIE is
weaker than the cross-source baseline** on the same 8,868 knockdowns: Pearson (all) -0.006 vs 0.015, Pearson (DE) 0.296
vs 0.460, sign accuracy (DE) 0.613 vs 0.750, AUPRC 0.023 vs 0.053; on strong knockdowns (n_de_bh >= 50) AUPRC is level
(0.128 vs 0.121) but correlations trail (0.044 vs 0.160). PIE therefore enters only as one evidence source (Phase E
compares on vs off); a PIE-only submission would be a weak fallback."""))
C.append(code("""PIE_LINES = ['A', 'K562', 'RPE1', 'HepG2', 'Jurkat', 'B', 'C']
if not all((V / f'lib/pie/contexts/{l}.h5ad').exists() for l in PIE_LINES):
    bg('pie_contexts', '"$PY" scripts/pie_contexts.py --lines ' + ' '.join(PIE_LINES)); wait('pie_contexts')
for line in PIE_LINES:                                    # line A first: it times the run
    if not (V / f'lib/pie/predictions/{line}.parquet').exists():
        bg(f'pie_{line}', f'bash scripts/pie_run.sh {line}'); wait(f'pie_{line}')
    if not (V / f'lib/pie_tables/pie_{line}.h5ad').exists():
        run(f'"$PY" scripts/pie_tables.py --lines {line}')
display(pd.DataFrame([json.loads((V / f'lib/pie/predictions/{l}_timing.json').read_text()) for l in PIE_LINES]))"""))
C.append(md("""### PIE on K562 (its clean test line) next to the cross-source baseline
Effect space, per knockdown then the median; truly DE = BH 5% in K562. Baseline = precision-weighted lfc of every other
source (T cells as Stim48hr, H1 included), score |z|; PIE = lfc_pred, score p_de."""))
C.append(code("""if not (V / 'results/analysis/k562_benchmark_pie.csv').exists():
    run('"$PY" scripts/k562_benchmark.py --pie lib/pie_tables/pie_K562.npz')
cols = ['r_all', 'r_de', 'sign_de', 'auprc', 'auprc_random']
rows = {}
for m in ('baseline', 'pie'):
    d = pd.read_csv(V / f'results/analysis/k562_benchmark_{m}.csv')
    both = set(pd.read_csv(V / 'results/analysis/k562_benchmark_pie.csv').target) & set(pd.read_csv(V / 'results/analysis/k562_benchmark_baseline.csv').target)
    d = d[d.target.isin(both)]
    rows[(m, 'all')] = d[cols].median(); rows[(m, 'n_de_bh>=50')] = d.loc[d.n_de_bh >= 50, cols].median()
display(pd.DataFrame(rows).T.round(3))"""))
C.append(md("""## Phase D: transformer build and smoke test

Code: `code/transformer/` (config, data, model, losses, train, evaluate), configs `code/transformer/configs/*.yaml`,
tests `tests/test_transformer.py`. Rows = (knockdown k, target line t); targets K562, RPE1, HepG2, Jurkat, HCT116,
HEK293T (validation mode) + H1 (final mode); CD4 T evidence only. Holdout 20 % of target genes (seeded); dev gene set =
our 300 targets + top genes to 5,000 by number of target lines measuring them, then mean |lfc/se|. Evidence per (gene,
source != t) + a PIE slot (off for HCT116 / HEK293T / H1); basal features from baseline v2; efficiency = median own-gene
lfc over evidence sources (+ N(0, 0.2) in training). Model: PIE gene features -> d = 128, evidence MLP + source
attention, 4-layer pre-norm TransformerEncoder (4 heads, FF 512, dropout 0.1, bf16); delta_hat = Linear + gamma x
precision-weighted evidence mean; p_de logit. Loss: weighted Huber + 0.1 (1 - cosine) + BCE (pos_weight <= 50).

**Result (19 Oct):** 8/8 tests pass on the real tables (no target-line data in its own evidence, no H1 in validation
mode, no PIE for HCT116/HEK293T/H1, held-out genes never trained on, efficiency never from t, DE labels = the tables'
n_de_bh, NaN-free tensors and correct masks, finite forward/backward). Tiny model overfits 64 examples: loss
7.39 -> 0.012-0.026 in 8,000 steps (600 steps were not enough: most knockdown effects are noise to memorise). 200 real
steps: 4.3 -> ~1.5, finite. Full-size model (1.63 M parameters) on the L4: 500 genes batch >= 256, 0.6 min/epoch;
**5,000 genes max batch 64 (128 OOM), 8.3 min/epoch**; 18,533 genes max batch 16 (32 OOM), 58 min/epoch. Dev runs
therefore use 5,000 genes, batch 32 (< 3 h for 20 epochs)."""))
C.append(code("""run('PYTHONPATH=code "$PY" -m pytest tests/test_transformer.py -q')"""))
C.append(code("""if not (V / 'results/vcc_2026_4/smoke.json').exists():
    bg('smoke', 'PYTHONPATH=code "$PY" -m transformer.train smoke --config code/transformer/configs/tiny.yaml --out results/vcc_2026_4/smoke.json')
    wait('smoke')
sm = json.loads((V / 'results/vcc_2026_4/smoke.json').read_text())
display(pd.DataFrame({'overfit first': sm['overfit']['first'], 'overfit last': sm['overfit']['last']}).T)
print('200 real steps:', sm['real'])"""))
C.append(code("""if not (V / 'results/vcc_2026_4/bench.json').exists():
    bg('bench', 'PYTHONPATH=code "$PY" -m transformer.train bench --config code/transformer/configs/dev.yaml --out results/vcc_2026_4/bench.json')
    wait('bench')
display(pd.DataFrame(json.loads((V / 'results/vcc_2026_4/bench.json').read_text())).round(3))"""))
C.append(md("""## Phase E: dev run and validation

Validation mode (H1 unseen), 5,000 genes (D's timing: 8.3 min/epoch, so 20 epochs < 3 h), batch 32. Two runs: PIE
evidence off (`dev_nopie`) and on (`dev_pie`, once the PIE tables for K562/RPE1/HepG2/Jurkat exist). Effect space on
held-out genes (all six target lines) and on all 300 H1 knockdowns, next to the baseline (precision-weighted evidence
mean). H1 pipeline: call p_de > tau, size by delta_hat, uncalled = no change -> existing generator (400 fresh control-based
cells per knockdown) -> official scorer (cell-eval2, vcc2026 profile) on the unfiltered H1 cells, on the same 75-target
reference as the earlier H1 results (best 0.118 = submission 3's recipe; 0.1315 = submission 4's, which did worse on the
leaderboard). PIE on vs off is judged on held-out K562/RPE1/HepG2/Jurkat genes (H1 cannot judge PIE).

**Result so far (20 Oct):** `dev_nopie` early-stopped at epoch 15 (best 12, held-out loss 0.792). Effect space, model vs
baseline: held-out genes Pearson (all) 0.083 vs 0.024, Pearson (DE) 0.681 vs 0.618, sign (DE) 0.888 vs 0.857, AUPRC
0.198 vs 0.111 (better in every line); H1 zero-shot Pearson (all) 0.156 vs 0.095, Pearson (DE) 0.417 vs 0.352, sign and
AUPRC level. **H1 pipeline, p_de > tau with delta_hat sizes: ~0 at every tau** (0.015, 0.011, 0.008, -0.009) vs the
evidence-mean arm 0.099: the model's delta_hat is shrunk toward 0 (median |delta| of calls 0.028 log2), so the generated
cells show almost no DE.

**Sizing fixed with the detection floor (20 Oct).** Same calls, each called gene raised to the smallest change a 400-cell
test detects (the submission-3 recipe's sizing): tau 0.5 / 0.7 / 0.9 + floor 0.067 / 0.037 / -0.006; top-N genes by p_de
+ floor: N = 1,000 0.088, 2,000 0.118, 3,000 0.136, **5,000 (every gene of the dev set) 0.1486**; alpha 4 before the floor
0.129 / 0.140 / 0.1475; alpha 16 < 0 (nmae collapses). Best previous: 0.1177 (submission 3 recipe), 0.1315 (submission 4's
call-everything recipe, which then did worse on the leaderboard: H1 has misranked call-everything recipes before).

**PIE on vs off** (`dev_pie` vs `dev_nopie`): best held-out loss 0.7918 vs 0.7917; on held-out genes of PIE's clean lines
(K562, RPE1, HepG2, Jurkat; 3,047 knockdowns) Pearson (all) 0.096 vs 0.097, Pearson (DE) 0.693 vs 0.684, sign (DE) 0.862
vs 0.847, AUPRC 0.221 vs 0.215 (sign and AUPRC better in every one of the four lines). Small but consistent: **PIE on**.

**Chosen for Phase F:** final mode with PIE evidence (`final_pie`, 12 epochs; `final_nopie` as backup), every gene of the
5,000-gene set called, sized by max(|delta_hat|, detection floor), alpha 1; genes outside the set = no change."""))
C.append(code("""RUNS = V / 'lib/vcc_2026_4/runs'
if not (RUNS / 'dev_nopie/done.json').exists():
    bg('train_dev_nopie', 'PYTHONPATH=code "$PY" -m transformer.train train --config code/transformer/configs/dev.yaml --name dev_nopie --pie off')
    wait('train_dev_nopie', every=300)
pie_ready = all((V / f'lib/pie_tables/pie_{l}.h5ad').exists() for l in ['K562', 'RPE1', 'HepG2', 'Jurkat'])
if pie_ready and not (RUNS / 'dev_pie/done.json').exists():
    bg('train_dev_pie', 'PYTHONPATH=code "$PY" -m transformer.train train --config code/transformer/configs/dev.yaml --name dev_pie --pie on')
    wait('train_dev_pie', every=300)
runs = [r for r in ('dev_nopie', 'dev_pie') if (RUNS / r / 'done.json').exists()]
print('runs:', runs)"""))
C.append(code("""import matplotlib.pyplot as plt
TRAIN, HELD = '#2a78d6', '#eb6834'          # categorical slots 1-2 (dataviz default palette)
fig, axes = plt.subplots(1, len(runs), figsize=(5.2 * len(runs), 3.6), sharey=True, squeeze=False)
for ax, r in zip(axes[0], runs):
    log = pd.read_csv(RUNS / r / 'log.csv')
    for col, color, label in (('train_loss', TRAIN, 'train'), ('val_loss', HELD, 'held-out genes')):
        ax.plot(log.epoch, log[col], color=color, lw=2, marker='o', ms=4, label=label)
        ax.annotate(f'{log[col].iloc[-1]:.3f}', (log.epoch.iloc[-1], log[col].iloc[-1]), xytext=(4, 0),
                    textcoords='offset points', va='center', fontsize=8, color='#444')
    best = log.val_loss.idxmin(); ax.axvline(log.epoch[best], color='#bbb', lw=1, ls=':')
    ax.set_title(f'{r}: best epoch {int(log.epoch[best])}', fontsize=10); ax.set_xlabel('epoch')
    ax.grid(axis='y', color='#e6e6e6', lw=0.8); ax.spines[['top', 'right']].set_visible(False)
axes[0][0].set_ylabel('loss (Huber + 0.1 cos + BCE)'); axes[0][0].legend(frameon=False, fontsize=9)
plt.tight_layout(); plt.show()"""))
C.append(code("""from IPython.display import display
eff = {}
for r in runs:
    if not (RUNS / r / 'effect_summary.json').exists():
        bg(f'eval_{r}', f'PYTHONPATH=code "$PY" -m transformer.evaluate effect --run lib/vcc_2026_4/runs/{r}'); wait(f'eval_{r}')
    for part in ('heldout', 'h1'):
        d = json.loads((RUNS / r / 'effect_summary.json').read_text())[part]
        for k, v in d.items():
            meth, sub = k.split('|'); eff[(r, part, sub, meth)] = v
effect = pd.DataFrame(eff).T[['r_all', 'r_de', 'sign_de', 'auprc', 'auprc_random', 'knockdowns']]
effect.index.names = ['run', 'rows', 'knockdowns', 'method']; display(effect.round(3))
per_line = []
for r in runs:
    h = pd.read_csv(RUNS / r / 'effect_heldout.csv')
    per_line.append(h.groupby(['line', 'method']).r_all.median().unstack().assign(run=r))
display(pd.concat(per_line).round(3))"""))
C.append(code("""# H1 grids actually run: dev_nopie h1_bench (tau, no floor), h1_bench_floor (tau / top-N + floor), h1_bench_g3 (top-N x alpha
# + floor); dev_pie h1_bench_floor. A run with no grid yet gets the floor grid (about 1.5 h on CPU).
h1 = []
for r in runs:
    grids = sorted((RUNS / r).glob('h1_bench*/scaled_scores.csv'))
    if not grids:
        bg(f'h1_{r}', f'PYTHONPATH=code "$PY" -m transformer.evaluate h1-pipeline --run lib/vcc_2026_4/runs/{r} --taus 0.5 0.7 0.9 '
                      f'--topn 1000 2000 5000 --floor --tag _floor --no-baseline'); wait(f'h1_{r}', every=300)
        grids = sorted((RUNS / r).glob('h1_bench*/scaled_scores.csv'))
    for g in grids:
        h1.append(pd.read_csv(g).assign(run=r, grid=g.parent.name))
h1 = pd.concat(h1).drop_duplicates(['run', 'arm']).set_index(['run', 'grid', 'arm']).sort_values('overall', ascending=False)
print('reference: model2_top2000_floor (submission 3 recipe) 0.1177; model2_all_alpha0.75 (submission 4) 0.1315')
display(h1.round(4))"""))
C.append(md("""## Phase F: submission

Final mode (H1 becomes a target and an evidence source), the chosen config and epoch count; predictions for A, B, C (our
300 targets); calling with the tau chosen on H1 -> generator -> `vcc prep` (`scripts/make_submission.sh`, needs ~34 GB
RAM + swap). The summary is shown before anything is submitted; submitting is a separate, manual step. Fallback: a
PIE-only submission (PIE lfc_pred / p_de for A-C, same calling + generator).

*Result: filled in when the phase runs.*"""))
C.append(code("""SETTINGS = RUNS / 'dev_nopie/h1_bench_g3/selected_generator.json'      # tf_top5000 + floor, alpha 1 (H1 0.1486)
chosen = json.loads(SETTINGS.read_text()); print(chosen['arm'], chosen['delta_array'], round(chosen['scaled_overall'], 4))
for r in ('final_pie', 'final_nopie'):
    if not (RUNS / r / 'done.json').exists():
        bg(f'train_{r}', f'PYTHONPATH=code "$PY" -m transformer.train train --config code/transformer/configs/final.yaml --name {r} --pie {"on" if r == "final_pie" else "off"}')
        wait(f'train_{r}', every=300)
    if not (RUNS / r / 'pred_top5000/C/tf_top5000.npy').exists():
        run(f'PYTHONPATH=code "$PY" -m transformer.evaluate predict-lines --run lib/vcc_2026_4/runs/{r} --taus --topn 5000 --out lib/vcc_2026_4/runs/{r}/pred_top5000')"""))
C.append(code("""SUBS = {r: V / f'results/vcc_2026_4_{r}_top5000' for r in ('final_pie', 'final_nopie')}
for r, out in SUBS.items():
    if not (out / 'prediction.vcc').exists():
        bg(f'sub_{r}', f'bash scripts/make_submission.sh lib/vcc_2026_4/runs/{r}/pred_top5000 tf_top5000 {SETTINGS} vcc_2026_4_{r}_top5000')
        wait(f'sub_{r}', every=300)
summary = []
for r, out in SUBS.items():
    rep = json.loads((out / 'generator_report.json').read_text())
    summary.append(dict(run=r, file=str(out / 'prediction.vcc'), size_gb=round((out / 'prediction.vcc').stat().st_size / 1e9, 2),
                        lines=list(rep['lines']), cells=sum(v['n_cells'] for v in rep['lines'].values()),
                        ontarget_measured={k: v['measured_ontarget'] for k, v in rep['lines'].items()},
                        array=rep['delta_array'], alpha=rep['settings']['alpha'], floor=rep['settings']['detect_floor']))
display(pd.DataFrame(summary))
print('Not submitted. Submitting is a separate manual step: vcc submit <file> -m <name> --wait')"""))

nb['cells'] = C
nb['metadata'] = {'kernelspec': {'display_name': 'Python 3 (venv)', 'language': 'python', 'name': 'python3'},
                  'language_info': {'name': 'python'}}
nbf.write(nb, str(Path.home() / 'vcc/vcc_2026_4.ipynb'))
print('written', len(C), 'cells')
