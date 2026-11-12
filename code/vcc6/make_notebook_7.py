"""Write vcc_2026_7.ipynb: config 1 end to end on one page -- 1 data loader, 2 model, 3 cell generator, 4 internal eval
(Eval A, HCT116final), 5 vcc upload (D/E/F). The stage-1 model code, the frozen generator (code/cellgen) and the frozen
harness (scripts/harness) are imported, never copied.

  ~/venv/bin/python code/vcc6/make_notebook_7.py
"""
from pathlib import Path

import nbformat

V = Path.home() / 'vcc'
nb = nbformat.v4.new_notebook()
md = lambda s: nb.cells.append(nbformat.v4.new_markdown_cell(s.strip('\n')))
code = lambda s: nb.cells.append(nbformat.v4.new_code_cell(s.strip('\n')))

md(r'''
# vcc_2026_7 — config 1 (frozen 2026-10-29)

**Config 1** is C+G (stage-1 arm C + Gaussian shrink, all 117 tables, the line's own tables excluded) with new magnitudes handed to the frozen generator.

**Values** (`CONFIG` in the setup cell; the same as `CONFIG1` in `code/vcc6/scorer_step16.py`):

| Setting | Value | Meaning |
|---|---|---|
| s_c | 0.5 | called genes (C+G's calls, K0) get s_c × C+G's LFC |
| m | 0.5 | non-called genes get control profile + m × c × Δ_pool |
| κ | 1 | c = z²/(z² + κ) for genes with ≥ 2 source lines; z = \|mean\| / (SD / √n) between source lines |
| c₁ | 1 | c for a gene with one source line (0 with none) |
| calls | K0 | C+G's own calls, unchanged |
| tables | F0 | all 117 effect tables |

Δ_pool is the sources' knockdown-minus-control pseudobulk, pooled by the same stage-1 machinery.

**Results, Eval A (H1 held out), scaled.** Config 1 has not yet been run end to end; its column is the noise-free proxy, which matched C+G's end-to-end overall to 0.0013. Section 4 gives the real score.

| Metric (scaled) | C+G, end-to-end run | C+G, proxy | Config 1, proxy (held out) |
|---|---|---|---|
| MSE | 0 (raw 1.496) | 0 | 0.113 |
| PDS | 0.723 | 0.724 | 0.808 |
| NMAE | 0.040 | 0.040 | 0.001 |
| Direction fidelity | 0.380 | 0.371 | 0.371 |
| Direction reach | 0.070 | 0.074 | 0.073 |
| Jaccard | 0.097 | 0.095 | 0.095 |
| **Overall** | **0.219** | 0.217 | **0.244** |

- **Sources:** `results/stage1/2026-10-28_scorer_space/tables/step11_heldout.csv` (arm i_lofo, H1) and `step9.json` (EvalA_CG); full precision in `results/7/config1_record.json`.
- **Frozen:** `code/vcc6/CONFIG1_FROZEN.md`, checked by `tests/test_config1.py` (CONFIG, this notebook's data-loader and model cells, `stage1.py`, `run_line.py`, the 117-table set). Toggling the run flags is fine; anything else breaks the freeze, and new work goes in `vcc_2026_8.ipynb`.
- **Freeze check:** this notebook's own cells, run for Eval A, reproduce the proxy's config 1 exactly: identical calls, means within 1.1e-7, proxy overall 0.24361700 both ways (`results/7/freeze_check/`).
- **Risk:** NMAE +0.001 is fragile (+0.008 and −0.006 on the two H1 knockdown halves).
- **D/E/F:** section 5 builds the upload with this config. D/E/F has no truth here; the table above is Eval A only.

**Sections:**
1. Data loader
2. Model
3. Cell generator
4. Internal eval (Eval A, HCT116final, against C+G)
5. VCC upload (D/E/F)

**Imported, not copied:** `code/vcc6/stage1.py` (model), `code/vcc6/run_line.py` and `gen_cells.py` (generation, scoring, file), the **frozen** `code/cellgen/cell_generator.py` and `scripts/harness/`, and `code/vcc6/boot_score.py`.

**Running:** about 5–6 h on vcc2026. Run it in tmux, since a dropped connection kills a notebook run:
`jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=-1 vcc_2026_7.ipynb --output results/7/vcc_2026_7.run.ipynb`

Each step skips work whose output already exists, so a rerun resumes.
''')
code(r'''
import gc, json, subprocess, sys, time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem
from IPython.display import display

V = Path.home() / "vcc"
for p in (V / "code/vcc6", V / "scripts", V / "scripts/harness", V / "code/cellgen"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
import stage1 as S
import run_line as RL

CONFIG = dict(s_c=0.5, m=0.5, kappa=1.0, c1=1.0)   # called scale, non-called scale, agreement kappa, one-source weight
RUN_EVAL = True      # section 4
RUN_UPLOAD = True    # section 5: build and vcc prep the D/E/F file
SUBMIT = False       # vcc submit -- only with Jose's explicit OK (the last entry counts)
MAKE_SWAP = False    # temporary 32 GB swap (sudo) for scoring and vcc prep; removed at the end
OUT = V / "results/7"; OUT.mkdir(parents=True, exist_ok=True)
PY, VCC = sys.executable, str(Path(sys.executable).parent / "vcc")

FE, DEF = V / "results/stage1/final_evals/harness/lib", V / "data_final/extracted"
LINES = {   # background = the line's own source tables, which never vote
    "EvalA": dict(controls=V / "lib/harness/evalA/context.h5ad", knockdowns=V / "results/target_lists/evalA_h1_292.txt",
                  background="H1", reference=V / "lib/harness/evalA/real.h5ad",
                  incumbent=V / "results/stage1/2026-10-27_incumbent/evalA/CG/e2e/salt_none/harness/evalA"),
    "HCT116final": dict(controls=FE / "evalHCT116final/context.h5ad", knockdowns=V / "results/target_lists/final_HCT116.txt",
                        background="HCT116", reference=FE / "evalHCT116final/real.h5ad",
                        incumbent=V / "results/stage1/e2e_final/CG/HCT116/salt_none/harness/evalHCT116final"),
    **{x: dict(controls=DEF / f"context_{x}.h5ad", knockdowns=V / "results/target_lists/final_DEF_300.txt",
               background=x, reference=None, incumbent=None) for x in "DEF"},
}
RUN = [l for l, s in LINES.items() if (s["reference"] and RUN_EVAL) or (not s["reference"] and RUN_UPLOAD)]
GENES = RL.panel_genes(); G_ = np.asarray(GENES, dtype=object)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S', time.gmtime())}] {msg}", flush=True)


def sh(cmd, logfile):
    """One step in its own process; stop on a non-zero exit."""
    Path(logfile).parent.mkdir(parents=True, exist_ok=True)
    log(f"{Path(str(cmd[1])).name if len(cmd) > 1 else cmd[0]} ... (log {logfile})")
    with open(logfile, "w") as f:
        rc = subprocess.run(list(map(str, cmd)), stdout=f, stderr=subprocess.STDOUT, cwd=V).returncode
    if rc:
        raise RuntimeError(f"exit {rc}: {Path(logfile).read_text()[-2000:]}")


def swap(on):
    sw = Path("/swapfile_v7")
    cmds = ([f"sudo fallocate -l 32G {sw}", f"sudo chmod 600 {sw}", f"sudo mkswap {sw}", f"sudo swapon {sw}"]
            if on and not sw.exists() else [f"sudo swapoff {sw}", f"sudo rm {sw}"] if not on and sw.exists() else [])
    for c in cmds:
        subprocess.run(c, shell=True, check=True, capture_output=True)


kds = lambda l: [k for k in Path(LINES[l]["knockdowns"]).read_text().split() if k]
todo = [l for l in RUN if not (OUT / l / "input.npz").exists()]
log(f"config {CONFIG}; lines {RUN}; model still to build for {todo}")
''')

md('''
## 1. Data loader
- **Effect tables:** all 117 effect tables (about 15 GB of RAM), loaded once and shared.
- **Δ rows:** each table's Δ rows (`log1p_bulk_50k_full` minus its control profile) for the knockdowns needed.
- **Controls:** each line's control profile, the same one generation recomputes: m0 (mean CPM), the 5-CPM gate and the control pseudobulk.
''')
code(r'''
if todo:
    TABLES = S.discover_tables([str(V / "lib/effects_v5"), str(V / "lib/effects_v5_extra")])
    BY_STEM = {p.stem: p for p in TABLES}
    SOURCES = S.load_tables(TABLES)

    need = set().union(*(kds(l) for l in todo))
    DSRC = []                                   # Delta copies of the sources (as stage1 sources: lfc := Delta)
    for src in SOURCES:
        p = BY_STEM[src.name]
        ctrl = pd.read_csv(p.with_name(src.name + "_baseline.csv"), index_col=0)["log1p_bulk_50k_full"].reindex(GENES).to_numpy(float)
        with h5py.File(p, "r") as f:
            kd = read_elem(f["obs"])["target_gene"].astype(str).to_numpy()
            rows = [i for i, k in enumerate(kd) if k in need]
            if not rows or not np.isfinite(ctrl).any():
                continue
            full, se = np.asarray(f["layers/log1p_bulk_50k_full"][rows], float), np.asarray(f["layers/se"][rows], np.float32)
        if not np.isfinite(full).any():
            continue
        D, k_, n = (full - ctrl).astype(np.float32), kd[rows], np.full(len(rows), np.nan, np.float32)
        if src.modality != "drug":
            k_, D, se, n = S._collapse_duplicate_knockdowns(k_, D, se, n)
        DSRC.append(S.Source(src.name, src.background, src.modality, src.reference, src.genes, k_, D, se, n, src.mult))

    CONTROLS = {}
    for l in todo:
        X = RL.read_controls(LINES[l]["controls"], GENES)
        m0, gate, _, prof = RL.control_profile(X, OUT / l / "controls")
        raw = np.asarray(X.sum(0)).ravel().astype(float)
        CONTROLS[l] = dict(m0=m0, gate=gate, cprof=np.log1p(5e4 * raw / raw.sum()))
    log(f"{len(SOURCES)} tables; {len(DSRC)} with Delta rows; controls for {todo}")
''')

md('''
## 2. Model
Per line, from the sources other than the line's own:
- **C+G:** calls and LFC, as `run_line predict --arm C --magnitude gaussian`.
- **Δ pool:** the same model on the Δ rows.
- **Agreement z per gene:** computed across the source lines.
- **Config 1:** turned into generator LFCs, so the frozen generator's mean inversion gives exactly config 1's means.

Writes `results/7/<line>/input.npz`.
''')
code(r'''
def config1(cg_lfc, cg_sig, dpool, n_bg, z, m0, cprof):
    z2 = np.nan_to_num(z.astype(float)) ** 2
    c = np.where((n_bg >= 2) & np.isfinite(z), z2 / (z2 + CONFIG["kappa"]), np.where(n_bg == 1, CONFIG["c1"], 0.0))
    mu = np.maximum(20.0 * np.expm1(cprof + CONFIG["m"] * c * dpool), 0.0)          # non-called: scorer-space mean (CPM)
    return np.where(cg_sig, CONFIG["s_c"] * cg_lfc, np.log2((mu + 1e-9) / (m0 + 1e-9))).astype(np.float32)


for l in todo:
    targets, ctl = kds(l), CONTROLS[l]
    key = S._line_key(LINES[l]["background"])
    others = [s for s in SOURCES if S._line_key(s.background) != key]
    basals = {}
    for s in others:
        bp = BY_STEM[s.name].with_name(s.name + "_baseline.csv")
        got = S.load_baseline(bp, s.genes) if bp.exists() else None
        if got is not None: basals[s.name] = (s.genes, got[0])
    votes, prior_votes = S.arm_votes(others, "C")
    calib, prior = S.fit_se_calibration(others, basals), S.fit_prior(others, votes=prior_votes)
    dvotes = S.BackgroundVotes([s for s in DSRC if s.name in {o.name for o in others}])

    K, G = len(targets), len(GENES)
    lfc_gen, sig, n_src = np.zeros((K, G), np.float32), np.zeros((K, G), bool), np.zeros((K, G), np.int16)
    for i, k in enumerate(targets):
        cg = S.predict_knockdown(others, k, G_, ctl["m0"], ctl["gate"], calib, prior, drop_gene=k, votes=votes, magnitude="gaussian")
        dp = S.predict_knockdown(None, k, G_, ctl["m0"], ctl["gate"], calib, prior, drop_gene=k, votes=dvotes, magnitude="gaussian")
        dpool = np.where(dp["n_sources"] > 0, dp["lfc"], 0.0)
        Y, _, _ = dvotes.gather(k, G_)
        n_bg, z = np.zeros(G, int), np.full(G, np.nan)
        if Y:
            Y = np.stack(Y); fin = np.isfinite(Y); n_bg = fin.sum(0)
            mean = np.where(fin, Y, 0.0).sum(0) / np.maximum(n_bg, 1)
            sd = np.sqrt(np.where(fin, (Y - mean) ** 2, 0.0).sum(0) / np.maximum(n_bg - 1, 1))
            z = np.where(n_bg >= 2, np.abs(mean) / np.maximum(sd / np.sqrt(np.maximum(n_bg, 1)), 1e-12), np.nan)
        lfc_gen[i] = config1(cg["lfc"], cg["de_sig"], dpool, n_bg, z, ctl["m0"], ctl["cprof"])
        sig[i], n_src[i] = cg["de_sig"], cg["n_sources"]
    np.savez_compressed(OUT / l / "input.npz", knockdowns=np.array(targets, dtype=object), genes=G_,
                        lfc=lfc_gen, de_sig=sig, n_sources=n_src)
    log(f"{l}: {K} knockdowns, {len(others)} tables vote; calls per knockdown, median {np.median(sig.sum(1)):.0f}")

for v in ("SOURCES", "DSRC", "TABLES", "BY_STEM", "CONTROLS"):    # free the tables before generation
    globals().pop(v, None)
gc.collect()
''')

md('''
## 3. Cell generator
- **Generator:** the frozen generator, run through `run_line predict --lfc-npz`.
- **Settings:** 400 cells per knockdown, at the controls' median library, salt "".
- **Throughput:** about 50 min per line on 4 CPUs.
''')
code(r'''
for l in RUN:
    if not (OUT / l / "e2e/salt_none/run.json").exists():
        s = LINES[l]
        sh([PY, V / "code/vcc6/run_line.py", "predict", "--controls", s["controls"], "--knockdowns", s["knockdowns"],
            "--background", s["background"], "--lfc-npz", OUT / l / "input.npz", "--salt", "", "--workers", 4,
            "--out", OUT / l / "e2e"], OUT / l / "generate.log")
    run = json.loads((OUT / l / "e2e/salt_none/run.json").read_text())
    log(f"{l}: generated; calls per knockdown requested {run['sig_requested_median']:.0f}, realised {run['sig_realized_median']:.0f}")
''')

md('''
## 4. Internal eval
- **Scoring:** the frozen harness scores Eval A and HCT116final.
- **Comparison:** paired knockdown bootstrap (2,000 resamples) against the incumbent C+G run on the same eval.
- **Adoption rule:** on Eval A, MSE > 0, every other component > 0, and overall above C+G with the CI's lower bound > −0.01. On HCT116, overall must be within its CI.
''')
code(r'''
if RUN_EVAL:
    swap(MAKE_SWAP)
    rows = []
    for l in [l for l in RUN if LINES[l]["reference"]]:
        s = LINES[l]; hd = OUT / l / "e2e/salt_none/harness" / s["reference"].parent.name
        if not (hd / "summary.json").exists():
            sh([PY, V / "code/vcc6/run_line.py", "score", "--pred", OUT / l / "e2e/salt_none/prediction.h5ad",
                "--reference", s["reference"]], OUT / l / "score.log")
        boot = OUT / l / "boot_vs_CG.json"
        if not boot.exists():
            sh([PY, V / "code/vcc6/boot_score.py", s["incumbent"], hd, "--out", boot], OUT / l / "boot.log")
        cmp_ = next(iter(json.loads(boot.read_text())["compare"].values()))
        new, old = (json.loads((d / "summary.json").read_text())["scaled"] for d in (hd, s["incumbent"]))
        for m in new:
            rows.append(dict(line=l, metric=m, CG=old[m], config1=new[m], diff=cmp_[m]["diff"],
                             ci_lo=cmp_[m]["ci"][0], ci_hi=cmp_[m]["ci"][1]))
    T = pd.DataFrame(rows); display(T.round(4))
    a = T[T.line == "EvalA"].set_index("metric")
    ok = (a.loc[a.index != "overall", "config1"] > 0).all() and a.loc["overall", "diff"] > 0 and a.loc["overall", "ci_lo"] > -0.01
    print("Eval A adoption rule:", "PASS" if ok else "FAIL -- C+G stays")
''')

md('''
## 5. VCC upload
- **Build:** D, E and F without controls and with `context` set, merged on disk (360,000 cells × 18,533 genes).
- **Validate:** `vcc prep` against the final gene list and `pert_counts.csv`. It peaks around 34 GiB, so use `MAKE_SWAP`.
- **Upload:** `vcc submit` only with `SUBMIT = True`, after Jose's explicit OK. The deadline is Nov 5 23:59 UTC, and the last entry counts.
''')
code(r'''
if RUN_UPLOAD:
    import anndata as ad
    import anndata.experimental as AE
    swap(MAKE_SWAP)
    fin = OUT / "final"; fin.mkdir(exist_ok=True)
    for x in "DEF":
        if not (fin / f"submission_{x}.h5ad").exists():
            RL.submission_file(OUT / x / "e2e/salt_none/prediction.h5ad", fin / f"submission_{x}.h5ad", context=x)
    merged, vccf = fin / "prediction_DEF.h5ad", fin / "prediction_DEF.vcc"
    if not merged.exists():
        AE.concat_on_disk([str(fin / f"submission_{x}.h5ad") for x in "DEF"], str(merged), join="inner", merge="same",
                          keys=list("DEF"), index_unique="-")
    a = ad.read_h5ad(merged, backed="r"); assert a.shape == (360_000, len(GENES)), a.shape
    log(f"merged {a.shape}, contexts {a.obs['context'].value_counts().to_dict()}"); a.file.close()
    if not vccf.exists():
        sh([VCC, "prep", merged, "-g", DEF / "gene_names.csv", "--perts", DEF / "pert_counts.csv", "-o", vccf, "--force"],
           fin / "vcc_prep.log")
    log(f"vcc prep OK: {vccf} ({vccf.stat().st_size / 1e9:.2f} GB)")
    if SUBMIT:
        sh([VCC, "submit", vccf, "-m", "LHY"], fin / "vcc_submit.log"); print((fin / "vcc_submit.log").read_text()[-1000:])
    else:
        print("Not submitted (SUBMIT = False).")
swap(False)
''')

out = V / 'vcc_2026_7.ipynb'
nbformat.write(nb, out)
print(out)
