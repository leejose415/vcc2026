"""S5: target overlap, same-lab vs cross-lab transfer within the training data, and where H1's off-panel mass lives.

  cd ~/vcc && PYTHONPATH=scripts/v5:scripts ~/venv/bin/python analysis/mse_metric/s5_transfer.py

Pseudobulk space = the scorer's log1p(5e4 x share), approximated per line from the cached per-knockdown sums of per-cell
CPM (share ~ mean CPM / 1e6, so P = log1p(0.05 x mean CPM)); its sampling variance by the delta method from sum CPM^2.
Cell lines (QC'd cells, kept constructs of kept knockdowns, as scripts/build_effects_v4.py pools them): K562, RPE1
(Replogle), HepG2, Jurkat (Nadig), H1 (Arc). Table lines: HCT116, HEK293T (X-Atlas/Orion pseudobulk sums; kept rows =
the effects_v4 table's). CD4 (Marson) has no per-knockdown means and is left out.
Pair score (source S -> destination D), shared kept knockdowns, genes measured in both, each row's own target left out:
prediction = D's control + S's Delta; per knockdown num = ||pred - P_D,k||^2 - C_D,k; den = ||P_D,k - P_D,ctrl||^2 -
C_D,k - C_D,ctrl; raw = sum num / sum den (control = 1). Also: the same with S's own sampling noise removed from num
(what a noise-free S measurement would score), and D's own mean response as the prediction (D's generic baseline).
"""
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem

V = Path.home() / 'vcc'
OUT = V / 'analysis/mse_metric'
sys.path.insert(0, str(V / 'scripts'))
from build_effects_v4 import panel  # noqa: E402

CELL = {'K562': 'replogle_k562_gw', 'RPE1': 'replogle_rpe1', 'HepG2': 'nadig_hepg2', 'Jurkat': 'nadig_jurkat', 'H1': 'arc_h1'}
TABLE = {'HCT116': 'orion_hct116', 'HEK293T': 'orion_hek293t'}
STUDY = {'K562': 'Replogle', 'RPE1': 'Replogle', 'HepG2': 'Nadig', 'Jurkat': 'Nadig', 'HCT116': 'Orion', 'HEK293T': 'Orion', 'H1': 'Arc'}
LAB = {'Replogle': 'Weissman', 'Nadig': 'Weissman', 'Orion': 'Xaira', 'Arc': 'Arc'}
F = 0.05
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 100)


def log(*a): print(*a, flush=True)


def finish(name, tg, s1, s2, n, c1, c2, cn, measured, genes):
    m = s1 / np.maximum(n, 1)[:, None]; m0 = c1 / max(cn, 1)
    vm = np.maximum(s2 / np.maximum(n, 1)[:, None] - m ** 2, 0) / np.maximum(n - 1, 1)[:, None]
    vm0 = np.maximum(c2 / max(cn, 1) - m0 ** 2, 0) / max(cn - 1, 1)
    P = np.log1p(F * m); P0 = np.log1p(F * m0)
    vP = (F / (1 + F * m)) ** 2 * vm; vP0 = (F / (1 + F * m0)) ** 2 * vm0
    P[:, ~measured] = np.nan; P0[~measured] = np.nan
    gp = {g: i for i, g in enumerate(genes)}
    d = dict(name=name, tg=np.asarray(tg), row={t: i for i, t in enumerate(tg)}, P=P.astype(np.float32), vP=vP.astype(np.float32),
             P0=P0, vP0=vP0, n=n, measured=measured, own=np.array([gp.get(t, -1) for t in tg]))
    d['D'] = d['P'] - P0[None, :].astype(np.float32)
    d['Dbar'] = np.nanmean(d['D'], 0)
    log(f'{name}: {len(tg)} kept knockdowns, {int(measured.sum())} measured genes, median cells {np.median(n):.0f}, controls {cn:.0f}')
    return d


def cell_line(name, src, var_panel, genes):
    qc = V / 'results/qc' / src
    guides = pd.read_csv(qc / 'guide_qc.csv'); kq = pd.read_csv(qc / 'knockdown_qc.csv')
    guides['target_gene_id'] = guides['target_gene_id'].astype(str); kq['target_gene_id'] = kq['target_gene_id'].astype(str)
    tot = np.load(qc / 'construct_totals.npz', allow_pickle=False)
    with h5py.File(V / f'lib/qc_cells/{src}_qc.h5ad', 'r') as f: var = read_elem(f['var'])
    src_ids = var.index.astype(str).to_numpy() if 'gene_name' in var else var['gene_id'].astype(str).to_numpy()
    pid = {g: i for i, g in enumerate(var_panel.gene_id) if g}; to_panel = np.array([pid.get(g, -1) for g in src_ids]); ok = to_panel >= 0
    kept = kq[kq.knockdown_keep.astype(bool)].reset_index(drop=True)
    kg = guides[guides.keep.astype(bool) & guides.target_gene_id.isin(set(kept.target_gene_id))]
    row_of = {t: i for i, t in enumerate(kept.target_gene_id)}
    Pm = sp.csr_matrix((np.ones(len(kg)), ([row_of[t] for t in kg.target_gene_id], kg.index.to_numpy())), shape=(len(kept), len(guides)))
    G = len(genes)
    def expand(X):
        out = np.zeros((X.shape[0], G)); out[:, to_panel[ok]] = X[:, ok]; return out
    s1, s2 = expand(Pm @ tot['s1'][:-1]), expand(Pm @ tot['s2'][:-1]); n = Pm @ tot['n'][:-1]
    c1, c2 = expand(tot['s1'][-1:])[0], expand(tot['s2'][-1:])[0]
    measured = np.zeros(G, bool); measured[to_panel[ok]] = True
    return finish(name, kept.target_gene.astype(str).to_numpy(), s1, s2, n, c1, c2, tot['n'][-1], measured, genes)


def table_line(name, src, genes):
    with h5py.File(V / f"lib/orion/{name}_pseudobulk.h5ad", 'r') as f:
        obs = read_elem(f['obs']); var = read_elem(f['var']); S1 = f['layers/sum_cpm']; S2 = f['layers/sum_cpm2']
        names = obs.index.astype(str).to_numpy(); pos = {t: i for i, t in enumerate(names)}
        with h5py.File(V / f'lib/effects_v4/{src}.h5ad', 'r') as g: kept = read_elem(g['obs'])['target_gene'].astype(str).to_numpy()
        kept = np.array([t for t in kept if t in pos]); ci = pos['non-targeting']
        rows = np.array(sorted(pos[t] for t in kept)); tg = names[rows]
        dense = lambda ds, r: ds[r]                                       # dense float32 layers; rows sorted
        s1 = dense(S1, rows).astype(np.float64); s2 = dense(S2, rows).astype(np.float64)
        c1 = dense(S1, [ci]).astype(np.float64)[0]; c2 = dense(S2, [ci]).astype(np.float64)[0]
    n = obs['n_cells'].to_numpy(float); measured = var['measured'].to_numpy(bool)
    assert list(var.index.astype(str)) == list(genes)
    return finish(name, tg, s1, s2, n[rows], c1, c2, n[ci], measured, genes)


def rows_r(A, B):
    A = A - A.mean(1, keepdims=True); B = B - B.mean(1, keepdims=True)
    return (A * B).sum(1) / np.sqrt((A * A).sum(1) * (B * B).sum(1))


def pair(S, Dl):
    shared = [t for t in S['tg'] if t in Dl['row']]
    if len(shared) < 10: return None
    iS = np.array([S['row'][t] for t in shared]); iD = np.array([Dl['row'][t] for t in shared])
    g = S['measured'] & Dl['measured']
    DS, DD = S['D'][iS][:, g].astype(np.float64), Dl['D'][iD][:, g].astype(np.float64)
    own = Dl['own'][iD]; gi = np.flatnonzero(g); colpos = {c: j for j, c in enumerate(gi)}
    mask = np.ones_like(DD, bool)
    for r, o in enumerate(own):
        if o in colpos: mask[r, colpos[o]] = False
    DS, DD = np.where(mask, DS, 0), np.where(mask, DD, 0)
    vD = np.where(mask, Dl['vP'][iD][:, g], 0).sum(1); vD0 = (Dl['vP0'][g][None, :] * mask).sum(1)
    vS = np.where(mask, S['vP'][iS][:, g], 0).sum(1) + (S['vP0'][g][None, :] * mask).sum(1)
    err = ((DS - DD) ** 2).sum(1)
    num = err - vD; num_src_clean = err - vD - vS; den = (DD ** 2).sum(1) - vD - vD0
    base = (np.where(mask, Dl['Dbar'][g][None, :], 0) - DD) ** 2
    num_base = base.sum(1) - vD
    dS = DS - np.where(mask, S['Dbar'][g][None, :], 0); dD = DD - np.where(mask, Dl['Dbar'][g][None, :], 0)
    return dict(source=S['name'], dest=Dl['name'], n_shared=len(shared), n_genes=int(g.sum()),
                median_r_Delta=float(np.nanmedian(rows_r(DS, DD))), median_r_d=float(np.nanmedian(rows_r(dS, dD))),
                raw_mse=float(num.sum() / den.sum()), raw_mse_source_noise_removed=float(num_src_clean.sum() / den.sum()),
                raw_mse_dest_mean_response=float(num_base.sum() / den.sum()),
                dest_noise_share_of_den=float(vD.sum() / ((DD ** 2).sum(1)).sum()),
                tag=('same study (' + STUDY[S['name']] + ')') if STUDY[S['name']] == STUDY[Dl['name']] else
                    ('same lab, different study' if LAB[STUDY[S['name']]] == LAB[STUDY[Dl['name']]] else 'cross-lab'))


def main():
    res = {}
    genes = pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str).to_numpy(); G = len(genes)
    var_panel, _ = panel()
    L = {k: cell_line(k, s, var_panel, genes) for k, s in CELL.items()}
    for k, s in TABLE.items(): L[k] = table_line(k, s, genes)
    # ---------------- 5a
    tg_ch = pd.read_csv(V / 'data/pert_counts.csv')['target_gene'].astype(str).tolist()
    all_tested = {}
    for k, s in CELL.items():
        all_tested[k] = set(pd.read_csv(V / 'results/qc' / s / 'guide_qc.csv')['target_gene'].astype(str))
    for k in TABLE:
        with h5py.File(V / f'lib/orion/{k}_pseudobulk.h5ad', 'r') as f: all_tested[k] = set(read_elem(f['obs']).index.astype(str))
    with h5py.File(V / 'lib/effects_v4/marson_cd4_stim48hr.h5ad', 'r') as f: cd4 = set(read_elem(f['obs'])['target_gene'].astype(str))
    ov = {k: dict(tested=sum(t in all_tested[k] for t in tg_ch), kept_qc=sum(t in L[k]['row'] for t in tg_ch)) for k in L}
    ov['CD4_Stim48hr (effects table)'] = dict(tested=None, kept_qc=sum(t in cd4 for t in tg_ch))
    ov['Replogle or Nadig (kept)'] = dict(kept_qc=sum(any(t in L[k]['row'] for k in ('K562', 'RPE1', 'HepG2', 'Jurkat')) for t in tg_ch))
    ov['none of the lines (kept)'] = dict(kept_qc=sum(not any(t in L[k]['row'] for k in L) and t not in cd4 for t in tg_ch))
    sys.path.insert(0, str(V / 'code/vcc_2026_3c'))
    from vcc3_generate import Line
    gates = {}
    for c in 'ABC':
        gates[c] = Line.from_h5ad(V / f'data/context_{c}.h5ad', list(genes)).gate
    RN = L['K562']['measured'] | L['RPE1']['measured'] | L['HepG2']['measured'] | L['Jurkat']['measured']
    meas = {'H1': L['H1']['measured'], 'K562': L['K562']['measured'], 'RPE1': L['RPE1']['measured'], 'HepG2': L['HepG2']['measured'],
            'Jurkat': L['Jurkat']['measured'], 'Replogle or Nadig': RN, 'HCT116/HEK293T (Orion)': L['HCT116']['measured'] | L['HEK293T']['measured']}
    gene_tab = {k: dict(all_18533=int(v.sum()), **{f'expressed_in_{c} (>5 CPM)': int((v & gates[c]).sum()) for c in 'ABC'}) for k, v in meas.items()}
    gene_tab['A-C panel / expressed'] = dict(all_18533=G, **{f'expressed_in_{c} (>5 CPM)': int(gates[c].sum()) for c in 'ABC'})
    gene_tab['in H1 but not Replogle/Nadig'] = dict(all_18533=int((meas['H1'] & ~RN).sum()), **{f'expressed_in_{c} (>5 CPM)': int((meas['H1'] & ~RN & gates[c]).sum()) for c in 'ABC'})
    res['5a'] = dict(targets=ov, genes=gene_tab, D_F='not available on disk (data/ holds A-C only; all three share one 300-target pert_counts.csv)')
    log(pd.DataFrame(ov).T.to_string()); log(pd.DataFrame(gene_tab).T.to_string())
    # ---------------- 5b
    rows = []
    for a in L:
        for b in L:
            if a == b: continue
            r = pair(L[a], L[b])
            if r: rows.append(r); log(f"{a:8s} -> {b:8s} n={r['n_shared']:5d} r_Delta {r['median_r_Delta']:.3f} r_d {r['median_r_d']:.3f} raw {r['raw_mse']:.3f} ({r['raw_mse_source_noise_removed']:.3f}) base {r['raw_mse_dest_mean_response']:.3f}")
    df = pd.DataFrame(rows); df.to_csv(OUT / 's5b_pairs.csv', index=False)
    summ = df.groupby('tag')[['median_r_Delta', 'median_r_d', 'raw_mse', 'raw_mse_source_noise_removed', 'raw_mse_dest_mean_response']].median()
    summ['pairs'] = df.groupby('tag').size()
    res['5b'] = dict(pairs=df.to_dict('records'), by_tag=summ.to_dict('index'))
    log(summ.round(3).to_string())
    # ---------------- 5c (H1's 75-target reference, scorer space)
    from pb import H1
    h = H1(); T = h.R - h.ctrl[None]
    E2 = T ** 2
    for i in h.rows:
        if h.own[i] >= 0: E2[i, h.own[i]] = 0
    groups = {'transformer 5,000 panel': h.run_set, 'Replogle/Nadig-measured, off panel': RN & ~h.run_set,
              'not measured by Replogle/Nadig, measured by Orion': ~RN & ~h.run_set & meas['HCT116/HEK293T (Orion)'],
              'measured only in H1 (of these lines)': ~RN & ~h.run_set & ~meas['HCT116/HEK293T (Orion)']}
    rk = pd.Series(h.ctrl).rank(method='first').to_numpy(); dec = np.minimum((rk - 1) * 10 // G, 9).astype(int)
    out5c = {}
    for q, lab in [(2, 'strong tercile'), (None, 'all 75')]:
        sel = np.arange(h.T) if q is None else np.flatnonzero(h.terc == q); tot = E2[sel].sum()
        out5c[lab] = dict(genes={k: dict(n_genes=int(m.sum()), share=float(E2[sel][:, m].sum() / tot)) for k, m in groups.items()},
                          deciles={f'decile {k + 1}': float(E2[sel][:, dec == k].sum() / tot) for k in range(10)},
                          noise_share_of_total=float(h.Creal[sel].sum() / tot))
    res['5c'] = out5c; log(json.dumps(out5c, indent=1))
    (OUT / 's5_results.json').write_text(json.dumps(res, indent=1, default=float))


if __name__ == '__main__':
    main()
