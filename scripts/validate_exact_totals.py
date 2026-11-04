"""Validate an exact-totals submission against the previous one (same model, same generator otherwise).

  python scripts/validate_exact_totals.py --old <old run dir> --new <new run dir> --pred <pred dir with A/ B/ C/> --array tf_top5000
Run dirs hold prediction.h5ad; the new one also expected_totals_<line>.npz (written by vcc3_generate.py --exact-totals).
Writes <new>/validation.json and prints a PASS/FAIL table:
  0  called genes per knockdown (the generator's call rule: nonzero calibrated delta, > 5 CPM in the line's controls,
     not the knocked-down gene), median and IQR per line
  1  structure: n_obs, n_vars, gene order, obs columns, targets per line; integer counts, min >= 0, max counts/cell <= 1e6
  2  totals: per knockdown, column sums == round(expected totals); max |error|
  3  pseudobulk: 10 random knockdowns per line, MSE of log1p CPM (x 5e4) against the expected profile, old vs new
  4  per-cell shape on the same knockdowns: detection rate, library sizes, per-gene variance (median ratio new / old)
  5  DE sanity: 3 strong + 3 weak knockdowns per line, Wilcoxon of the 400 cells vs the line's real controls (per-cell
     CPM, genes with mean control CPM > 5, BH 0.05 per knockdown): significant genes old vs new; FAIL if new > 1.5 x old
     for most of them
"""
import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import mannwhitneyu

V = Path.home() / 'vcc'
sys.path.insert(0, str(V / 'code/vcc_2026_3c'))
from vcc2_common import expression_gate                                   # noqa: E402

LINES = ['A', 'B', 'C']


def bh(p):
    n = len(p); o = np.argsort(p); q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]; out = np.empty(n); out[o] = np.minimum(q, 1); return out


def rows_of(path):
    a = ad.read_h5ad(path, backed='r')
    try:
        obs = a.obs[['target_gene', 'context']].astype(str).copy(); var = list(a.var_names.astype(str)); cols = list(a.obs.columns)
    finally:
        a.file.close()
    idx = {}
    for (c, t), g in obs.groupby(['context', 'target_gene'], sort=False).indices.items():
        g = np.sort(g); assert g[-1] - g[0] + 1 == len(g), 'knockdown rows must be contiguous'
        idx[(c, t)] = (int(g[0]), int(g[-1]) + 1)
    return obs, var, cols, idx


def read_rows(path, a, b):
    with h5py.File(path, 'r') as f:
        X = f['X']; ip = X['indptr'][a:b + 1]; lo, hi = int(ip[0]), int(ip[-1])
        return sp.csr_matrix((X['data'][lo:hi], X['indices'][lo:hi], ip - lo), shape=(b - a, int(X.attrs['shape'][1])))


def controls(line, genes):
    a = ad.read_h5ad(V / f'data/context_{line}.h5ad')
    a = a[a.obs['target_gene'].astype(str) == 'non-targeting']
    assert list(a.var_names.astype(str)) == list(genes)
    X = a.X.tocsr().astype(np.float64); lib = np.asarray(X.sum(1)).ravel()
    cpm = sp.diags(1e6 / np.maximum(lib, 1)) @ X
    mean = np.asarray(cpm.mean(0)).ravel()
    return cpm.tocsc(), mean


def de_count(args):
    cells, ctrl_cols, gene_idx = args
    lib = np.asarray(cells.sum(1)).ravel(); cpm = (sp.diags(1e6 / np.maximum(lib, 1)) @ cells.astype(np.float64)).tocsc()
    x = cpm[:, gene_idx].toarray(); y = ctrl_cols
    p = mannwhitneyu(x, y, axis=0, alternative='two-sided', method='asymptotic').pvalue
    p = np.where(np.isfinite(p), p, 1.0)
    return int((bh(p) < 0.05).sum())


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--old', required=True); ap.add_argument('--new', required=True)
    ap.add_argument('--pred', required=True); ap.add_argument('--array', default='tf_top5000')
    ap.add_argument('--workers', type=int, default=6); a = ap.parse_args()
    old, new = Path(a.old) / 'prediction.h5ad', Path(a.new) / 'prediction.h5ad'
    genes = pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str).tolist(); gp = {g: i for i, g in enumerate(genes)}
    res = {}
    # 0 called genes
    called = {}
    for L in LINES:
        ctrl_cpm, ctrl_mean = controls(L, genes); gate = expression_gate(ctrl_mean)
        D = np.load(Path(a.pred) / L / f'{a.array}.npy'); t = np.load(Path(a.pred) / L / 'targets.npy').astype(str)
        n = []
        for i, tg in enumerate(t):
            c = (D[i] != 0) & gate
            if tg in gp: c[gp[tg]] = False
            n.append(int(c.sum()))
        n = np.array(n); called[L] = dict(median=float(np.median(n)), q25=float(np.percentile(n, 25)), q75=float(np.percentile(n, 75)),
                                         min=int(n.min()), max=int(n.max()), gate_genes=int(gate.sum()))
    res['0_called_genes'] = dict(rule='nonzero calibrated delta (every gene of the 5,000-gene set with delta != 0, sized by '
                                      'max(|delta|, detection floor)), > 5 CPM in the line controls, not the knocked-down gene', **called)
    # 1 structure
    oo, ov, oc, oidx = rows_of(old); no, nv, nc, nidx = rows_of(new)
    s = dict(n_obs=(len(oo), len(no)), n_vars=(len(ov), len(nv)), gene_order_same=ov == nv, obs_columns_same=oc == nc,
             targets_same=all(set(oo[oo.context == L].target_gene) == set(no[no.context == L].target_gene) for L in LINES),
             rows_same=list(oo.index) == list(no.index) if hasattr(oo, 'index') else None)
    mx, mn, nonint = 0.0, 0.0, 0
    with h5py.File(new, 'r') as f:
        X = f['X']; ip = X['indptr'][:]
        for s0 in range(0, len(ip) - 1, 20000):
            e0 = min(s0 + 20000, len(ip) - 1); lo, hi = int(ip[s0]), int(ip[e0])
            d = X['data'][lo:hi]; nonint += int((d != np.round(d)).sum()); mn = min(mn, float(d.min(initial=0)))
            M = sp.csr_matrix((d, X['indices'][lo:hi], ip[s0:e0 + 1] - lo), shape=(e0 - s0, len(nv)))
            mx = max(mx, float(np.asarray(M.sum(1)).max()))
    s.update(non_integer_values=nonint, min_value=mn, max_counts_per_cell=mx)
    s['PASS'] = bool(s['n_obs'][0] == s['n_obs'][1] and s['n_vars'][0] == s['n_vars'][1] and s['gene_order_same'] and s['obs_columns_same']
                     and s['targets_same'] and nonint == 0 and mn >= 0 and mx <= 1e6)
    res['1_structure'] = s
    # 2 totals, 3 pseudobulk, 4 shape
    rng = np.random.default_rng(2026); tot_err = 0.0; pb = []; shape = []; strength = {}
    for L in LINES:
        z = np.load(Path(a.new) / f'expected_totals_{L}.npz'); E = z['expected'].astype(float); tl = z['targets'].astype(str)
        for i, tg in enumerate(tl):
            r0, r1 = nidx[(L, tg)]; Y = read_rows(new, r0, r1)
            tot_err = max(tot_err, float(np.abs(np.asarray(Y.sum(0)).ravel() - np.maximum(0, np.rint(E[i]))).max()))
        Dl = np.load(Path(a.pred) / L / f'{a.array}.npy'); tp = np.load(Path(a.pred) / L / 'targets.npy').astype(str)
        strength[L] = {tg: float(np.abs(Dl[k][Dl[k] != 0]).mean()) if (Dl[k] != 0).any() else 0.0 for k, tg in enumerate(tp)}
        for i in rng.choice(len(tl), 10, replace=False):
            tg = tl[i]; Xo = read_rows(old, *oidx[(L, tg)]); Xn = read_rows(new, *nidx[(L, tg)])
            lc = lambda v: np.log1p(v / max(v.sum(), 1) * 5e4)
            e = lc(E[i]); so = lc(np.asarray(Xo.sum(0)).ravel().astype(float)); sn = lc(np.asarray(Xn.sum(0)).ravel().astype(float))
            pb.append(dict(line=L, target=tg, mse_old=float(np.mean((so - e) ** 2)), mse_new=float(np.mean((sn - e) ** 2))))
            Do, Dn = Xo.toarray().astype(float), Xn.toarray().astype(float)
            mid = (Do.mean(0) >= 0.5)
            lo_, ln_ = Do.sum(1), Dn.sum(1)
            shape.append(dict(line=L, target=tg, detection_old=float((Do > 0).mean()), detection_new=float((Dn > 0).mean()),
                              lib_median_old=float(np.median(lo_)), lib_median_new=float(np.median(ln_)),
                              lib_iqr_old=float(np.subtract(*np.percentile(lo_, [75, 25]))),
                              lib_iqr_new=float(np.subtract(*np.percentile(ln_, [75, 25]))),
                              var_ratio_median=float(np.median(Dn[:, mid].var(0) / np.maximum(Do[:, mid].var(0), 1e-12)))))
    res['2_totals'] = dict(max_abs_error=tot_err, PASS=tot_err == 0)
    pbd = pd.DataFrame(pb)
    res['3_pseudobulk'] = dict(per_line={L: dict(mse_old_median=float(g.mse_old.median()), mse_new_median=float(g.mse_new.median()))
                                         for L, g in pbd.groupby('line')}, PASS=bool((pbd.mse_new < pbd.mse_old).all()),
                               rows=pb)
    sd = pd.DataFrame(shape)
    res['4_shape'] = dict(detection_old=float(sd.detection_old.mean()), detection_new=float(sd.detection_new.mean()),
                          lib_median_old=float(sd.lib_median_old.median()), lib_median_new=float(sd.lib_median_new.median()),
                          lib_iqr_old=float(sd.lib_iqr_old.median()), lib_iqr_new=float(sd.lib_iqr_new.median()),
                          var_ratio_median=float(sd.var_ratio_median.median()),
                          PASS=bool(abs(sd.detection_new.mean() - sd.detection_old.mean()) < 0.01 and 0.9 < sd.var_ratio_median.median() < 1.1),
                          rows=shape)
    # 5 DE sanity
    de = []
    for L in LINES:
        ctrl_cpm, ctrl_mean = controls(L, genes); gi = np.flatnonzero(ctrl_mean > 5); Yc = ctrl_cpm[:, gi].toarray()
        order = sorted(strength[L], key=strength[L].get); picks = [('weak', t) for t in order[:3]] + [('strong', t) for t in order[-3:]]
        jobs = []
        for kind, tg in picks:
            jobs.append((read_rows(old, *oidx[(L, tg)]), Yc, gi)); jobs.append((read_rows(new, *nidx[(L, tg)]), Yc, gi))
        with ProcessPoolExecutor(a.workers) as ex:
            counts = list(ex.map(de_count, jobs))
        for k, (kind, tg) in enumerate(picks):
            de.append(dict(line=L, target=tg, kind=kind, strength=strength[L][tg], sig_old=counts[2 * k], sig_new=counts[2 * k + 1]))
        print(f'DE sanity {L} done', flush=True)
    ded = pd.DataFrame(de); ratio = (ded.sig_new + 1) / (ded.sig_old + 1)
    res['5_de_sanity'] = dict(rows=de, median_ratio_new_over_old=float(ratio.median()),
                              share_over_1p5x=float((ded.sig_new > 1.5 * ded.sig_old.clip(lower=1)).mean()),
                              PASS=bool((ded.sig_new > 1.5 * ded.sig_old.clip(lower=1)).mean() <= 0.5))
    (Path(a.new) / 'validation.json').write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != 'rows'} for k, v in res.items()}, indent=1, default=float))
    print(pd.DataFrame(de).to_string(index=False))


if __name__ == '__main__':
    main()
