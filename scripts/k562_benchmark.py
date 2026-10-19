"""K562 effect-space benchmark: a cross-source baseline (and, with --pie, PIE) against effects_v4 replogle_k562_gw. Read only;
writes results/analysis/k562_benchmark_<method>.csv (per knockdown) and prints the medians.

Truth: K562 lfc over genes measured in K562; truly DE = BH 5% over those genes (two-sided normal p from lfc / se), the
knocked-down gene excluded, as for n_de_bh.
Baseline: per knockdown and gene, the 1/se^2-weighted mean of lfc over every other source that has the knockdown (RPE1,
HepG2, Jurkat, HCT116, HEK293T, H1, CD4 T Stim48hr; the three T-cell states count once, as Stim48hr); z = mean / sqrt(1 / sum w).
PIE (--pie npz with targets, genes, lfc, p_de): lfc_pred and p_de as given.
Per knockdown, over K562-measured genes (a gene the method does not predict counts as lfc 0 and score 0):
  r_all  Pearson of lfc;  r_de  Pearson over truly DE genes (>= 3);  sign_de  sign accuracy on truly DE genes;
  auprc  average precision for truly DE genes (score |z| for the baseline, p_de for PIE).
Knockdowns: every K562 knockdown the method covers (the baseline needs >= 1 other source); medians for all, and for
those strong in K562 (n_de_bh >= 50).
"""
import argparse
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem
from scipy.stats import norm
from sklearn.metrics import average_precision_score

V = Path.home() / 'vcc'
E = V / 'lib/effects_v4'
OUT = V / 'results/analysis'
TRUTH = 'replogle_k562_gw'
OTHERS = ['replogle_rpe1', 'nadig_hepg2', 'nadig_jurkat', 'orion_hct116', 'orion_hek293t', 'arc_h1', 'marson_cd4_stim48hr']
STRONG = 50


def key(obs):
    ids = obs['target_gene_id'].to_numpy(); sym = obs['target_gene'].astype(str).to_numpy()
    return np.array([str(t) if isinstance(t, str) and t.startswith('ENSG') else f'symbol:{s}' for t, s in zip(ids, sym)])


def load(source, rows=None):
    with h5py.File(E / f'{source}.h5ad', 'r') as f:
        obs = read_elem(f['obs']); k = key(obs)
        if rows is None:
            return obs, k
        idx = np.sort(np.flatnonzero(np.isin(k, rows)))
        return obs.iloc[idx], k[idx], f['layers/lfc'][idx].astype(np.float32), f['layers/se'][idx].astype(np.float32)


def bh(p):
    n = len(p); o = np.argsort(p); q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]; out = np.empty(n); out[o] = np.minimum(q, 1); return out


def metrics(true, de, pred, score):
    r = dict(n_genes=len(true), n_de=int(de.sum()))
    r['r_all'] = float(np.corrcoef(true, pred)[0, 1]) if np.std(pred) > 0 else np.nan
    if de.sum() >= 3:
        r['r_de'] = float(np.corrcoef(true[de], pred[de])[0, 1]) if np.std(pred[de]) > 0 else np.nan
    else:
        r['r_de'] = np.nan
    r['sign_de'] = float((np.sign(pred[de]) == np.sign(true[de])).mean()) if de.any() else np.nan
    r['auprc'] = float(average_precision_score(de, score)) if de.any() else np.nan
    r['auprc_random'] = float(de.mean())
    return r


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--pie', type=Path, help='npz: targets, genes, lfc, p_de')
    ap.add_argument('--sample', type=int, default=0, help='random knockdowns (0 = all)'); a = ap.parse_args()
    genes = pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str).to_numpy(); gpos = {g: i for i, g in enumerate(genes)}
    tobs, tkey = load(TRUTH)
    pick = tkey
    if a.sample:
        pick = np.random.default_rng(2026).choice(tkey, min(a.sample, len(tkey)), replace=False)
    tobs, tkey, tl, ts = load(TRUTH, pick)
    G = len(genes); T = len(tkey); pos = {k: i for i, k in enumerate(tkey)}
    W = np.zeros((T, G), np.float32); WL = np.zeros((T, G), np.float32); nsrc = np.zeros(T, int)
    for s in OTHERS:
        obs, k, l, se = load(s, tkey)
        rows = np.array([pos[x] for x in k]); w = np.where(np.isfinite(l), 1.0 / np.square(se), 0).astype(np.float32)
        W[rows] += w; WL[rows] += np.where(np.isfinite(l), l * w, 0); nsrc[rows] += 1
        print(f'{s}: {len(k)} of the {T} K562 knockdowns', flush=True)
    with np.errstate(invalid='ignore', divide='ignore'):
        base_l = np.where(W > 0, WL / W, np.nan); base_z = np.where(W > 0, base_l * np.sqrt(W), np.nan)
    pie = None
    if a.pie:
        z = np.load(a.pie, allow_pickle=False); pg = {g: i for i, g in enumerate(z['genes'].astype(str))}
        pt = {t: i for i, t in enumerate(z['targets'].astype(str))}; cols = np.array([pg.get(g, -1) for g in genes])
        pie = (pt, cols, z['lfc'], z['p_de'])
    out = {'baseline': [], 'pie': []}
    for i in range(T):
        l, se = tl[i], ts[i]; meas = np.isfinite(l)
        tgt = gpos.get(str(tobs['target_gene'].iloc[i]))
        if tgt is not None: meas[tgt] = False
        p = 2 * norm.sf(np.abs(l[meas] / se[meas])); de = bh(p) < 0.05; true = l[meas].astype(float)
        info = dict(target=tkey[i], target_gene=str(tobs['target_gene'].iloc[i]), n_de_bh=int(tobs['n_de_bh'].iloc[i]))
        if nsrc[i]:
            bl = np.nan_to_num(base_l[i][meas]); bz = np.nan_to_num(np.abs(base_z[i][meas]))
            out['baseline'].append(dict(info, n_sources=int(nsrc[i]), covered=float(np.isfinite(base_l[i][meas]).mean()),
                                        **metrics(true, de, bl, bz)))
        if pie is not None and str(tobs['target_gene'].iloc[i]) in pie[0]:
            j = pie[0][str(tobs['target_gene'].iloc[i])]; cols = pie[1][meas]; ok = cols >= 0
            pl = np.zeros(meas.sum()); ps = np.zeros(meas.sum())
            pl[ok] = pie[2][j][cols[ok]]; ps[ok] = pie[3][j][cols[ok]]
            out['pie'].append(dict(info, covered=float(ok.mean()), **metrics(true, de, pl, ps)))
    OUT.mkdir(parents=True, exist_ok=True)
    for m, rows in out.items():
        if not rows: continue
        df = pd.DataFrame(rows); df.to_csv(OUT / f'k562_benchmark_{m}.csv', index=False)
        cols = ['r_all', 'r_de', 'sign_de', 'auprc', 'auprc_random']
        print(f'\n{m}: {len(df)} knockdowns (median per-knockdown values)')
        print(pd.DataFrame({'all': df[cols].median(), f'n_de_bh >= {STRONG} (n={int((df.n_de_bh >= STRONG).sum())})':
                            df.loc[df.n_de_bh >= STRONG, cols].median()}).round(3).to_string())


if __name__ == '__main__':
    main()
