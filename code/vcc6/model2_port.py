"""Model 2's pooling (the v3 leaderboard recipe, +0.091, rank 749; commit c208334, vcc2_meta.py) on the current
pipeline, against C+G, differing only in how sources are pooled into a predicted LFC (Jose 2026-10-27, revised task).

Model 2 pooling, exactly as the recipe: one source per cell line, 8 lines (v3's SOURCES), the target line held out;
per-gene DerSimonian-Laird tau2 pooled over the panel targets seen in >= 2 lines (own column excluded, floored at 0,
median where not estimable); re = sum w lfc / sum w with w = 1 / (SE^2 + tau2); no shrinkage. Changes, only these:
  tables   each line's panel source (lib/sources/<line>/effects_panel.h5ad, deleted 10-21) -> that line's v5 table(s)
           (Orion: the single-cell table and its disjoint pseudobulk complement together are the line)
  calls    C+G's rule: stage1.predictive_calls (BH at 0.05 on |LFC| / the target's 400-cell measurement SE, over the
           controls' 5-CPM gate, own gene excluded), with C+G's own measurement SE for the same controls
  others   non-called genes as C+G: kept at their predicted LFC; a gene no source measured is 0 and never called
           (stage1.prepare_generator_inputs does this from n_sources)
  panel    v3 pooled tau2 over its 575 challenge + validation targets; here the targets being predicted
           (D/E/F 300 + Eval A 292 = 564)
Both arms are then generated through run_line predict --lfc-npz (the frozen generator; p_sign only orders a saved
reach list the generator never reads, so this is C+G's own generation path) and scored by the frozen harness.

  cg       ~/venv/bin/python code/vcc6/model2_port.py cg --controls CTX --knockdowns LIST --background H1 --out DIR
           C+G (stage1 arm C + gaussian) exactly as run_line predict builds it, without generating: DIR/CG/
  port     ~/venv/bin/python code/vcc6/model2_port.py port --controls CTX --knockdowns LIST --holdout H1 --out DIR
           Model 2 pooling + C+G's calls (needs DIR/CG): DIR/M2/
  compare  ~/venv/bin/python code/vcc6/model2_port.py compare --out DIR        Step 3 (no cells)
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

V = Path.home() / 'vcc'
for p in (V / 'code/vcc6', V / 'scripts'):
    sys.path.insert(0, str(p))
import stage1 as S                                                                         # noqa: E402

V5 = {'HCT116': ['orion_hct116', 'orion_hct116_pseudobulk'], 'HEK293T': ['orion_hek293t', 'orion_hek293t_pseudobulk'],
      'K562_gwps': ['replogle_k562_gw'], 'CD4_Stim48hr': ['marson_cd4_stim48hr'], 'H1': ['arc_h1'],
      'RPE1': ['replogle_rpe1'], 'HepG2': ['nadig_hepg2'], 'Jurkat': ['nadig_jurkat']}
PANEL = [V / 'data_final/extracted/pert_counts.csv', V / 'results/target_lists/evalA_h1_292.txt']


def log(msg):
    import time
    print(f'[{time.strftime("%H:%M:%S")}] {msg}', flush=True)


def panel_targets():
    t = list(pd.read_csv(PANEL[0])['target_gene'].astype(str))
    t += [k for k in PANEL[1].read_text().split() if k]
    return list(dict.fromkeys(t))


def load_line(line, targets, genes):
    """One line's v5 table(s) as targets x genes (lfc, se); NaN where absent. Several tables of one line hold disjoint
    knockdowns (Orion single-cell + pseudobulk complement); the first table holding a target supplies it."""
    by = {p.stem: p for p in S.discover_tables([str(V / 'lib/effects_v5'), str(V / 'lib/effects_v5_extra')])}
    lfc = np.full((len(targets), len(genes)), np.nan, np.float32); se = lfc.copy()
    filled = np.zeros(len(targets), bool); used = {}
    for name in V5[line]:
        s = S.load_source(by[name])
        idx = np.array([s.col_of.get(g, -1) for g in genes]); ok = idx >= 0
        n = 0
        for i, t in enumerate(targets):
            rows = s.rows_of.get(t)
            if filled[i] or not rows: continue
            assert len(rows) == 1, f'{name}: {t} has {len(rows)} rows'
            lfc[i, ok] = s.lfc[rows[0], idx[ok]]; se[i, ok] = s.se[rows[0], idx[ok]]; filled[i] = True; n += 1
        used[name] = n
        del s
    return lfc, se, used


def tau2_dl(lfc, se, own_col):
    """vcc2_meta.tau2_dl, verbatim."""
    _, T, G = lfc.shape
    Q, df, C = np.zeros(G), np.zeros(G), np.zeros(G)
    for t in range(T):
        y, v = lfc[:, t].astype(np.float64), se[:, t].astype(np.float64) ** 2
        ok = np.isfinite(y) & np.isfinite(v) & (v > 0)
        if own_col[t] >= 0: ok[:, own_col[t]] = False
        k = ok.sum(0); use = k >= 2
        if not use.any(): continue
        w = np.where(ok, 1 / np.where(ok, v, 1), 0); y = np.where(ok, y, 0)
        sw = w.sum(0); mu = np.divide((w * y).sum(0), sw, out=np.zeros(G), where=sw > 0)
        Q[use] += (w * (y - mu) ** 2).sum(0)[use]; df[use] += k[use] - 1
        C[use] += (sw - np.divide((w ** 2).sum(0), sw, out=np.zeros(G), where=sw > 0))[use]
    tau2 = np.full(G, np.nan); est = C > 0
    tau2[est] = np.maximum(0, (Q[est] - df[est]) / C[est])
    return tau2, df


def combine(lfc, se, tau2):
    """vcc2_meta.combine, verbatim: sources x rows x genes -> (re, re_se, k)."""
    y = lfc.astype(np.float64); total = se.astype(np.float64) ** 2 + tau2
    ok = np.isfinite(y) & np.isfinite(total)
    k = ok.sum(0); y = np.where(ok, y, 0)
    w = np.where(ok, 1 / np.where(ok, total, 1), 0); sw = w.sum(0)
    with np.errstate(divide='ignore', invalid='ignore'):
        out = ((w * y).sum(0) / sw, np.sqrt(1 / sw))
    return tuple(np.where(k > 0, x, np.nan).astype(np.float32) for x in out) + (k.astype(np.uint8),)


def control_gate(controls, genes, out):
    import run_line as RL
    (out / 'controls').mkdir(parents=True, exist_ok=True)
    return RL.control_profile(RL.read_controls(controls, genes), out / 'controls')


def write_input(path, want, genes, lfc, sig, n_src):
    np.savez_compressed(path, knockdowns=np.array(want, dtype=object), genes=np.array(genes, dtype=object),
                        lfc=lfc.astype(np.float32), de_sig=sig, n_sources=n_src.astype(np.int16))


def cg(a):
    """run_line predict's model step, verbatim, then its arrays as an --lfc-npz input."""
    import run_line as RL
    out = Path(a.out) / 'CG'; out.mkdir(parents=True, exist_ok=True)
    genes = RL.panel_genes()
    m0, gate, _, prof = control_gate(a.controls, genes, out)
    dirs = [str(V / 'lib/effects_v5'), str(V / 'lib/effects_v5_extra')]
    line = S.load_line_for_controls(dirs, a.background, genes, m0, gate)
    S.predict_line(line, str(out / 'model'), knockdowns_file=a.knockdowns, arm='C', magnitude='gaussian')
    z = np.load(out / 'model' / 'predictions.npz', allow_pickle=True)
    want = [k for k in Path(a.knockdowns).read_text().split() if k]
    assert list(map(str, z['knockdowns'])) == want
    write_input(out / 'input.npz', want, genes, z['lfc'], z['de_sig'], z['n_sources'])
    log(f'C+G: {len(want)} knockdowns, calls median {np.median(z["de_sig"].sum(1)):.0f}')


def port(a):
    import run_line as RL
    root = Path(a.out); out = root / 'M2'; out.mkdir(parents=True, exist_ok=True)
    genes = RL.panel_genes()
    want = [k for k in Path(a.knockdowns).read_text().split() if k]
    targets = panel_targets()
    assert not [k for k in want if k not in targets]
    lines = [l for l in V5 if l != a.holdout]
    log(f'Model 2 pooling: {len(lines)} lines vote {lines} (held out: {a.holdout}); panel {len(targets)} targets')
    L, E, used = [], [], {}
    for line in lines:
        x, e, u = load_line(line, targets, genes); L.append(x); E.append(e); used[line] = u
        log(f'  {line}: {u}')
    lfc, se = np.stack(L), np.stack(E); del L, E
    gp = {g: i for i, g in enumerate(genes)}
    tau2, df = tau2_dl(lfc, se, np.array([gp.get(t, -1) for t in targets]))
    est = np.isfinite(tau2); fallback = float(np.median(tau2[est]))
    tau2_used = np.where(est, tau2, fallback)
    log(f'tau2: estimated for {est.sum():,} genes ({(tau2[est] == 0).sum():,} at 0), median {fallback:.4f}')
    rows = [targets.index(k) for k in want]
    sub_l = lfc[:, rows].copy(); sub_s = se[:, rows]
    for j, k in enumerate(want):
        if k in gp: sub_l[:, j, gp[k]] = np.nan                            # own gene is not averaged
    re, re_se, k_bg = combine(sub_l, sub_s, tau2_used)
    # C+G's call rule with C+G's measurement SE for these controls (stage1.predict_knockdown's call)
    c = np.load(root / 'CG' / 'model' / 'predictions.npz', allow_pickle=True)
    assert list(map(str, c['knockdowns'])) == want and list(map(str, c['genes'])) == genes
    gate = pd.read_csv(root / 'CG' / 'controls' / 'controls_baseline.csv', index_col=0)['de_tested'].to_numpy(bool)
    D = np.where(np.isfinite(re), re, 0.0)
    sig = np.zeros(D.shape, bool)
    for i, k in enumerate(want):
        t = gate.copy()
        if k in gp: t[gp[k]] = False
        sig[i], _ = S.predictive_calls(D[i], c['se_measurement'][i], t)
    write_input(out / 'input.npz', want, genes, D, sig, k_bg)
    np.savez_compressed(out / 'pool.npz', re=re, re_se=re_se, n_lines=k_bg, tau2=tau2_used)
    rep = dict(recipe='v3 Model 2 pooling (commit c208334 vcc2_meta.py) on effects_v5, C+G calls, non-called kept',
               holdout=a.holdout, lines=lines, tables_used=used, panel_targets=len(targets), knockdowns=len(want),
               tau2_estimated=int(est.sum()), tau2_zero=int((tau2[est] == 0).sum()), tau2_median=fallback,
               calls=dict(mean=float(sig.sum(1).mean()), median=float(np.median(sig.sum(1)))))
    (out / 'model2_v5.json').write_text(json.dumps(rep, indent=1, default=str) + '\n')
    log(json.dumps(rep, default=str))


def compare(a):
    """Step 3: calls, |LFC| of called genes, per-knockdown Pearson r of the two LFCs over genes any source measured
    in either arm (and the target's gate), and the Jaccard of the call sets."""
    root = Path(a.out)
    m, c = (np.load(root / x / 'input.npz', allow_pickle=True) for x in ('M2', 'CG'))
    gate = pd.read_csv(root / 'CG' / 'controls' / 'controls_baseline.csv', index_col=0)['de_tested'].to_numpy(bool)
    res = {}
    for name, z in (('Model2_port', m), ('CG', c)):
        n = z['de_sig'].sum(1); l = np.abs(z['lfc'][z['de_sig']])
        res[name] = dict(calls_mean=float(n.mean()), calls_median=float(np.median(n)),
                         called_abs_lfc_median=float(np.median(l)) if l.size else None,
                         abs_lfc_median_measured=float(np.median(np.abs(z['lfc'][z['n_sources'] > 0]))))
    r, jac, r_gate = [], [], []
    for i in range(len(m['knockdowns'])):
        meas = (m['n_sources'][i] > 0) | (c['n_sources'][i] > 0)
        x, y = m['lfc'][i][meas].astype(float), c['lfc'][i][meas].astype(float)
        r.append(np.corrcoef(x, y)[0, 1] if x.std() > 0 and y.std() > 0 else np.nan)
        g = meas & gate; x, y = m['lfc'][i][g].astype(float), c['lfc'][i][g].astype(float)
        r_gate.append(np.corrcoef(x, y)[0, 1] if x.std() > 0 and y.std() > 0 else np.nan)
        a_, b_ = m['de_sig'][i], c['de_sig'][i]; u = (a_ | b_).sum()
        jac.append((a_ & b_).sum() / u if u else np.nan)
    q = lambda v: dict(median=float(np.nanmedian(v)), q25=float(np.nanquantile(v, 0.25)), q75=float(np.nanquantile(v, 0.75)))
    res['pearson_r_measured'] = q(r); res['pearson_r_measured_and_gated'] = q(r_gate); res['call_jaccard'] = q(jac)
    (root / 'step3_compare.json').write_text(json.dumps(res, indent=1) + '\n')
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)
    x = sub.add_parser('cg'); x.add_argument('--controls', required=True); x.add_argument('--knockdowns', required=True)
    x.add_argument('--background', required=True); x.add_argument('--out', required=True)
    x = sub.add_parser('port'); x.add_argument('--controls', required=True); x.add_argument('--knockdowns', required=True)
    x.add_argument('--holdout', required=True, choices=list(V5)); x.add_argument('--out', required=True)
    x = sub.add_parser('compare'); x.add_argument('--out', required=True)
    a = p.parse_args()
    {'cg': cg, 'port': port, 'compare': compare}[a.cmd](a)
