"""Write vcc_2026_6.ipynb: the end-to-end runner of the stage-1 RPE1 dev run (Steps 0-3). All model code is in code/vcc6/.

  ~/venv/bin/python code/vcc6/make_notebook.py
"""
from pathlib import Path

import nbformat

V = Path.home() / 'vcc'
nb = nbformat.v4.new_notebook()
md = lambda s: nb.cells.append(nbformat.v4.new_markdown_cell(s.strip('\n')))
code = lambda s: nb.cells.append(nbformat.v4.new_code_cell(s.strip('\n')))

md(r'''
# vcc_2026_6 — Stage 1, end to end (RPE1 dev run)

This notebook runs and shows the stage-1 dev run: the precision-weighted transfer model, the frozen cell generator and the frozen cell-eval2 harness. **No model code lives here.** One copy of everything is in `code/vcc6/`:

| file | what |
|---|---|
| `stage1.py` | the model: sources, pooling (arms A / B / C), prior, shrinkage, calls, generator handoff, self-test |
| `gen_cells.py` | a whole arm through the frozen generator (`code/cellgen`), in parallel, seed = `knockdown_seed(name, salt)` |
| `rpe1_eval.py` | the RPE1dev eval (300 knockdowns) built and scored by the frozen harness (`scripts/harness`, not edited) |
| `step0_identity.py`, `step1_oracle.py`, `step1_analyze.py`, `step2_votes.py`, `step3_inputs.py` | Steps 0–3 |
| `run_line.py` | the end-to-end path an unseen line runs: controls + knockdown list → model → generator → file → harness |
| `pipeline.sh` | Steps 1–3 in order, logged (the way to run them for real) |
| `report.py` | the tables below |

Outputs: `results/stage1/2026-10-27_rpe1/`. The 300 knockdowns are fixed (`results/target_lists/rpe1_dev_300.txt`, drawn once with numpy seed 2026). Held out: the whole RPE1 line (every table whose background is RPE1) as a source; the 300 are what is scored. H1 is a source here; it is held out only when H1 is the target.

**Running.** Each step cell reads existing results. Set its `RUN[...]` flag to recompute. Steps 1 and 3 take hours on vcc2026 (the generator is about 30 s per knockdown); a kernel disconnect kills a notebook run, so for a full rerun prefer `bash code/vcc6/pipeline.sh` in the background and use this notebook to read the results.
''')
code(r'''
import json, subprocess, sys
from pathlib import Path
import pandas as pd
from IPython.display import Markdown, display

V = Path.home() / "vcc"
C6 = V / "code/vcc6"
R = V / "results/stage1/2026-10-27_rpe1"
PY = str(Path.home() / "venv/bin/python")
sys.path.insert(0, str(C6))
RUN = dict(step0=False, step1=False, step2=False, step3=False)
pd.set_option("display.width", 200, "display.max_columns", 40)

def run(name, *args):
    """One logged step (results/.../logs/<name>.log); stops on failure."""
    (R / "logs").mkdir(parents=True, exist_ok=True)
    with open(R / "logs" / f"{name}.log", "w") as fh:
        rc = subprocess.run([*args], stdout=fh, stderr=subprocess.STDOUT).returncode
    print(f"{name}: exit {rc}")
    if rc: raise RuntimeError(f"{name} failed, see {R / 'logs' / (name + '.log')}")

def show_json(path, keys=None):
    d = json.loads(Path(path).read_text())
    print(json.dumps({k: d[k] for k in keys} if keys else d, indent=1, default=str)[:6000])
''')
md(r'''
## Model

`code/vcc6/stage1.py`, moved verbatim from this notebook on 2026-10-27. Changes since: `knockdown_seed(name, salt="")` (the default reproduces the old seeds), `calibration_stats` → `summary["fit"]`, the Step 2 arms (`BackgroundVotes`, `combine_within`, `arm_votes`, `predict_line(arm=...)`), and `load_line_for_controls` for `run_line.py`. Arm A is unchanged: its predictions on 20 RPE1 knockdowns are byte-identical to the notebook's code (Step 0).
''')
code(r'''
import stage1 as S
S.self_test()
''')
md(r'''
## Step 0 — one copy of the model, a fixed list, the eval

- (a) `predictions.npz` on 20 RPE1 knockdowns from the notebook's model cells (verbatim copy kept in `step0/stage1_notebook_cells_verbatim.py`) and from `stage1.py`, byte-compared.
- (b) `stage1_model.py` → `archive/2026-10-26/` (the working copy had already been deleted; the surviving snapshot is there).
- (c) the 300 list, and the RPE1dev reference built by the frozen harness (`build_cells` / `finish`).
''')
code(r'''
if RUN["step0"]:
    run("step0_sample", PY, str(C6 / "rpe1_eval.py"), "sample")
    run("step0_build_eval", PY, str(C6 / "rpe1_eval.py"), "build", "--threads", "4")
    run("step0_before", PY, str(C6 / "step0_identity.py"), "before", str(R / "step0/stage1_notebook_cells_verbatim.py"))
    run("step0_after", PY, str(C6 / "step0_identity.py"), "after")
    run("step0_compare", PY, str(C6 / "step0_identity.py"), "compare")
show_json(R / "step0/identity.json", ["file_bytes_identical", "member_bytes_identical", "summary_keys_added"])
print((V / "archive/2026-10-26/README.md").read_text())
kd = pd.read_csv(R / "step1/knockdowns_300.csv", index_col=0) if (R / "step1/knockdowns_300.csv").exists() else None
display(kd.describe() if kd is not None else "knockdown table is written by Step 1's analysis")
''')
md(r'''
## Step 1 — oracle diagnosis (report only)

Arms on the 300, scored by the frozen harness: **A** the reference's own cells, **B** zero, **C** the table's lfc and de_sig through `prepare_generator_inputs` and the frozen generator (seed `knockdown_seed(name)`), **D** no cells: the scorer's pseudobulk of control mean × the table's fold change (D2: the table's own `log1p_bulk_50k_matched`).
''')
code(r'''
if RUN["step1"]:
    run("step1_arms", PY, str(C6 / "step1_oracle.py"), "arms", "--workers", "4")
    run("step1_analyze", PY, str(C6 / "step1_analyze.py"))
import report as RP
t = RP.summary_table()
display(t[t.index.str.startswith("step1")].T)
s1 = json.loads((R / "step1/step1.json").read_text())
display(pd.DataFrame({a: s1[a]["sums"] | {"norm": s1[a]["norm_ratio_of_sums"], "median r(gen, real)": s1[a]["median_pearson_gen_real"]}
                      for a in ("A", "B", "C")}))
display(pd.DataFrame({k: s1[k] for k in ("D", "D2_table_bulk_matched")}).T)
show_json(R / "step1/step1.json", ["C_slopes", "checks"])
''')
md('### Step 1 by effect-size tercile, gene groups, shared error')
code(r'''
for a in ("A", "B", "C"):
    display(Markdown(f"**arm {a}**")); display(pd.DataFrame(s1["by_tercile"][a]).T)
cs = s1["C_splits"]
display(pd.Series(cs["panel"], name="arm C share of raw squared error")); display(pd.Series(cs["expression_decile"], name="by control-expression decile"))
print("called vs uncalled:", cs["called_vs_uncalled"], " shared error share:", round(cs["shared_error_share"], 3))
display(pd.DataFrame(cs["shared_top_genes"]))
''')
md(r'''
## Step 2 — how sources vote (effect space)

**A** vcc_2026_6 as it is (every table row votes). **B** one vote per background (DerSimonian–Laird within a background, variance floored at its best single row), no drug tables, `fit_prior` on cross-background pairs. **C** B plus Tahoe, one vote per background, prior as in B.
''')
code(r'''
if RUN["step2"]:
    run("step2", PY, str(C6 / "step2_votes.py"))
display(pd.read_csv(R / "step2/backgrounds.csv").drop(columns="tables_list"))
show_json(R / "step2/structure.json", ["n_backgrounds", "fit_prior_A", "tahoe_knockdowns_in_300", "tahoe_share_median_over_touched", "tahoe_share_bottom_decile_over_touched", "orion"])
a2 = json.loads((R / "step2/arms.json").read_text())
rows = []
for arm, v in a2.items():
    for sub in ("all_300", "tahoe_touched"):
        es = v[sub]["effect_space"] or {}; f = v[sub]["fit"]
        rows.append(dict(arm=arm, subset=sub, knockdowns=v[sub]["knockdowns"], nmae=es.get("nmae"), nmae_zero=es.get("nmae_zero"), mse=es.get("mse"),
                         sign_acc_true_de=es.get("sign_acc_true_de"), auroc=es.get("auroc"), n_called=v[sub]["median_n_called"], n_true=v[sub]["median_n_true"],
                         call_count_pearson=es.get("call_count_pearson"), r_shrunk=f["shrunk"]["pearson_r"], calib_slope_shrunk=f["shrunk"]["slope_true_on_pred"],
                         transmission_shrunk=f["shrunk"]["slope_pred_on_true"], r_unshrunk=f["unshrunk"]["pearson_r"], tau2=v["prior"]["tau2"], v0=v["prior"]["v0"], pi=v["prior"]["pi"]))
display(pd.DataFrame(rows).set_index(["arm", "subset"]))
''')
md(r'''
### Step 2 result and decision (2026-10-27)

| | A every row votes | B one vote per background, no drugs | C B + Tahoe |
|---|---|---|---|
| NMAE (zero = 1.0) | 0.9992 | 0.9988 | 0.9988 |
| r shrunk / unshrunk | 0.053 / 0.100 | 0.065 / 0.097 | 0.065 / 0.097 |
| calls (true median 407) | 12.5 | 21 | 21 |
| prior π / τ² | 0.005 / 1e-4 | 0.020 / 1e-4 | 0.020 / 1e-4 |
| end to end, scaled overall | — | −0.20 | −0.20 |

The three arms are within noise of each other: Tahoe covers 94 genes, 1 of the 300 here (6 of H1's 300, 1 of the competition list), so C = B. In all three the slab weight shrinks effects about 300-fold. **Magnitude test** (no refit, arm B arrays): G = Gaussian shrink only, `lfc_unshrunk × v0 / (v0 + v_model)`, gives r 0.099, calls 403.5, and **+0.122** scaled overall end to end (S −0.197, zero −0.405, salt noise range 0.0035).

**Decision: the stage-1 baseline is arm C with the Gaussian shrink (no slab).**
''')
md(r'''
## Step 3 — end to end, the path an unseen line runs

`run_line.py predict` gets only RPE1's control cells (the eval's context half) and the 300 list; it computes the quantiles, baseline, gate and library size from them, runs the model (or the zero / oracle arrays), the generator with `knockdown_seed(name, salt)`, writes the file and runs `vcc prep`. `run_line.py score` runs the frozen harness. Zero, oracle, B and C, salts "", s1, s2.
''')
code(r'''
if RUN["step3"]:
    ev, lst = R / "harness/lib/evalRPE1dev", V / "results/target_lists/rpe1_dev_300.txt"
    run("step3_inputs", PY, str(C6 / "step3_inputs.py"), "inputs")
    run("step3_profile", PY, str(C6 / "step3_inputs.py"), "profile")
    dirs = ["--effects-dir", str(V / "lib/effects_v5"), "--effects-dir", str(V / "lib/effects_v5_extra")]
    salts = ["--salt", "", "--salt", "s1", "--salt", "s2"]
    for arm in ("zero", "oracle", "B", "C"):
        out = R / "step3" / (arm if arm in ("zero", "oracle") else f"arm{arm}")
        src = ["--lfc-npz", str(R / f"step3/inputs/{arm}.npz")] if arm in ("zero", "oracle") else dirs + ["--arm", arm]
        run(f"step3_{arm}_predict", PY, str(C6 / "run_line.py"), "predict", "--controls", str(ev / "context.h5ad"), "--knockdowns", str(lst),
            "--background", "RPE1", *src, *salts, "--out", str(out))
        for s in ("none", "s1", "s2"):
            run(f"step3_{arm}_{s}_score", PY, str(C6 / "run_line.py"), "score", "--pred", str(out / f"salt_{s}/prediction.h5ad"), "--reference", str(ev / "real.h5ad"))
show_json(R / "step3/profile_check.json")
t = RP.summary_table()
display(t[t.index.str.startswith("step3")])
display(RP.salt_spread(t))
display(RP.tercile_table())
''')
md('## Report')
code(r'''
# standalone: runs without the setup cell above
from pathlib import Path
from IPython.display import Markdown, display
rep = Path.home() / "vcc/results/stage1/2026-10-27_rpe1/REPORT.md"
display(Markdown(rep.read_text()) if rep.exists() else Markdown(f"`{rep}` is not written yet."))
''')
nbformat.validate(nb)
nbformat.write(nb, V / 'vcc_2026_6.ipynb')
print(f'{V / "vcc_2026_6.ipynb"}: {len(nb.cells)} cells')
