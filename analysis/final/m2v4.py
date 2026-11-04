"""Step 1: Model 2 v4 on the LFC tables (lib/effects_v4, scripts/build_effects_v4.py).

  cd ~/vcc && PYTHONPATH=analysis/v6lin:scripts/v5 ~/venv/bin/python analysis/final/m2v4.py

Sources: K562 genome-wide, RPE1, HepG2, Jurkat, HCT116, HEK293T, CD4 Stim48hr (one T-cell vote); H1 is added as a
source for A-C only. Row QC weight q: knockdown < 30% dropped, 30-50% and unmeasurable 0.5 (knockdown % from the
target's own lfc). Per target x gene: weights q / (SE^2 + tau2_g); prediction = weighted mean; SE_pred =
sqrt(sum w^2 (SE^2 + tau2)) / sum w; re_z = |prediction| / SE_pred; n_src; sign_agree = share of sources whose lfc has
the prediction's sign. Own gene excluded; genes with no source: 0.
tau2_g (method of moments, pooled over targets with >= 2 sources): mean over targets of
[mean over sources of (y - weighted mean)^2 - mean over sources of SE^2], floored at 0.
"""
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from v6lin import LINES, Line, V

OUT = V / 'analysis/final'
CHUNK = 128


def log(*a): print(time.strftime('%H:%M:%S'), *a, flush=True)


def stack(lines, names, ts, gp):
    Y = np.stack([lines[n].raw(np.array([lines[n].row.get(t, -1) for t in ts]), 'lfc') for n in names]).astype(np.float64)
    S = np.stack([lines[n].raw(np.array([lines[n].row.get(t, -1) for t in ts]), 'se') for n in names]).astype(np.float64)
    Q = np.stack([np.array([lines[n].w[lines[n].row[t]] if t in lines[n].row else 0.0 for t in ts]) for n in names])
    for i, t in enumerate(ts):
        if t in gp: Y[:, i, gp[t]] = np.nan
    ok = np.isfinite(Y) & np.isfinite(S) & (Q[:, :, None] > 0)
    return Y, S, Q, ok


def tau2_mom(lines, names, genes):
    gp = {g: i for i, g in enumerate(genes)}; G = len(genes)
    multi = sorted({t for n in names for t in lines[n].row if sum(t in lines[m].row for m in names) >= 2})
    acc = np.zeros(G); cnt = np.zeros(G)
    for a in range(0, len(multi), CHUNK):
        ts = multi[a:a + CHUNK]; Y, S, Q, ok = stack(lines, names, ts, gp)
        k = ok.sum(0); w = np.where(ok, Q[:, :, None] / np.maximum(S, 1e-3) ** 2, 0); sw = w.sum(0)
        mu = np.divide((w * np.nan_to_num(Y)).sum(0), sw, out=np.zeros_like(sw), where=sw > 0)
        dev = np.where(ok, (np.nan_to_num(Y) - mu) ** 2, 0).sum(0); se2 = np.where(ok, S ** 2, 0).sum(0)
        use = k >= 2
        val = np.where(use, dev / np.maximum(k, 1) - se2 / np.maximum(k, 1), 0)
        acc += np.where(use, val, 0).sum(0); cnt += use.sum(0)
    t2 = np.maximum(np.where(cnt > 0, acc / np.maximum(cnt, 1), 0), 0)
    log(f'tau2 MoM over {len(multi)} targets with >= 2 sources: median {np.median(t2[cnt > 0]):.4f}, '
        f'mean {t2[cnt > 0].mean():.4f}, share > 0 {(t2[cnt > 0] > 0).mean():.3f}')
    return t2, cnt


def predict(lines, names, targets, genes, t2):
    gp = {g: i for i, g in enumerate(genes)}; G = len(genes); n = len(targets)
    out = {k: np.zeros((n, G), np.float32) for k in ('re', 're_se', 're_z', 'n_src', 'sign_agree')}
    for a in range(0, n, CHUNK):
        ts = targets[a:a + CHUNK]; Y, S, Q, ok = stack(lines, names, ts, gp)
        tot = S ** 2 + t2[None, None]; w = np.where(ok, Q[:, :, None] / np.where(ok, tot, 1), 0); sw = w.sum(0)
        mu = np.divide((w * np.nan_to_num(Y)).sum(0), sw, out=np.zeros_like(sw), where=sw > 0)
        se = np.divide(np.sqrt((w ** 2 * np.where(ok, tot, 0)).sum(0)), sw, out=np.full_like(sw, np.inf), where=sw > 0)
        k = ok.sum(0); agree = np.divide((ok & (np.sign(np.nan_to_num(Y)) == np.sign(mu)[None])).sum(0), k, out=np.zeros_like(sw), where=k > 0)
        sl = slice(a, a + len(ts))
        out['re'][sl] = mu; out['re_se'][sl] = np.where(np.isfinite(se), se, np.inf); out['re_z'][sl] = np.where(sw > 0, np.abs(mu) / se, 0)
        out['n_src'][sl] = k; out['sign_agree'][sl] = agree
    return out


def main():
    genes = pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str).to_numpy()
    lines = {k: Line(k, s, genes) for k, s in LINES.items()}; log('loaded sources')
    h1 = Line('H1', 'arc_h1', genes); lines_f = dict(lines, H1=h1)
    rep = {}
    t2v, cv = tau2_mom(lines, list(lines), genes)
    rep['tau2_validation'] = dict(median=float(np.median(t2v[cv > 0])), mean=float(t2v[cv > 0].mean()), share_positive=float((t2v[cv > 0] > 0).mean()))
    h1t = [str(t) for t in h1.tg]
    pv = predict(lines, list(lines), h1t, genes, t2v)
    np.savez(OUT / 'm2v4_h1.npz', targets=np.array(h1t), tau2=t2v, **pv); log('H1 predictions saved')
    t2f, cf = tau2_mom(lines_f, list(lines_f), genes)
    rep['tau2_final'] = dict(median=float(np.median(t2f[cf > 0])), mean=float(t2f[cf > 0].mean()))
    tg = pd.read_csv(V / 'data/pert_counts.csv')['target_gene'].astype(str).tolist()
    pf = predict(lines_f, list(lines_f), tg, genes, t2f)
    np.savez(OUT / 'm2v4_abc.npz', targets=np.array(tg), tau2=t2f, **pf); log('A-C predictions saved')
    rep['coverage'] = dict(h1_with_source=int((pv['n_src'].max(1) > 0).sum()), abc_with_source=int((pf['n_src'].max(1) > 0).sum()),
                           qc=dict({k: dict(rows=l.n, dropped=int((l.w == 0).sum()), half=int((l.w == 0.5).sum())) for k, l in lines_f.items()}))
    (OUT / 'm2v4_report.json').write_text(json.dumps(rep, indent=1, default=float)); print(json.dumps(rep, indent=1, default=float))


if __name__ == '__main__':
    main()
