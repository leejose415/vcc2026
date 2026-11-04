"""Scorer-space study part 2, Step 11 inputs (Jose 2026-10-28, results/stage1/2026-10-28_scorer_space/), on vcc2027.

No table load: the Delta votes are rebuilt from the per-source Delta files (delta/, plus delta_nb/ for complex
partners) and the tables' se rows read by h5py, each source a stand-in carrying only its metadata (stage1.table_meta);
the prior is the fold's own (cache/fold_<line><F>.npz), and the Delta pool is stage1's pool_effects + gaussian_posterior,
as predict_knockdown computes it. Recomputing the plain Delta pool must reproduce the fold cache's (checked).

  partners --filter F     S3: response cosine between knockdowns from K562 genome-wide Delta (all its rows; gated
                          genes), RPE1 (essential-wide) for the K562 fold, where K562 is not a source; the needed
                          knockdowns' 10 nearest. Delta rows for partner knockdowns of every table -> delta_nb/.
  pool --filter F         per fold line: per-(k, g) source agreement (number of backgrounds, sign agreement with the
                          pooled sign, agreement z = |mean| / (between-source SD / sqrt n)); the Delta pool with the
                          n nearest partners' votes added at weight w_nb (vote variance / w_nb), n in {3, 10},
                          w_nb in {0.25, 0.5, 1}; the fold sources' top-50 Delta PCs (uncentred, the line's gated
                          genes). Writes cache/s11_<line><F>.npz and tables/step11_partners*.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

V = Path.home() / 'vcc'
sys.path.insert(0, str(V / 'code/vcc6'))
import scorer_study as SS                                                                  # noqa: E402

OUT = SS.OUT
NB_N, NB_W = [3, 10], [0.25, 0.5, 1.0]
N_PC = 50


def tables_by_stem():
    import stage1 as S
    import screen as SC
    return {p.stem: p for p in S.discover_tables(SC.DIRS)}


def full_delta(path):
    """Delta for every row of a table (log1p_bulk_50k_full - its control profile), duplicate knockdowns averaged; the
    table's gated genes (de_tested, > 5 CPM)."""
    import h5py
    from anndata.io import read_elem
    genes = list(pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str))
    b = pd.read_csv(path.with_name(path.stem + '_baseline.csv'), index_col=0)
    ctrl = b['log1p_bulk_50k_full'].reindex(genes).to_numpy(np.float64)
    gate = (b['de_tested'].reindex(genes).fillna(False).astype(bool) & (b['control_mean_cpm'].reindex(genes) > 5)).to_numpy()
    with h5py.File(path, 'r') as f:
        kd = read_elem(f['obs'])['target_gene'].astype(str).to_numpy()
        full = np.concatenate([f['layers/log1p_bulk_50k_full'][s:s + 1000] for s in range(0, len(kd), 1000)])
    D = (full - ctrl[None, :]).astype(np.float32)
    u, inv = np.unique(kd, return_inverse=True)
    S_ = np.zeros((len(u), D.shape[1]), np.float64); C_ = np.zeros_like(S_)
    np.add.at(S_, inv, np.nan_to_num(D)); np.add.at(C_, inv, np.isfinite(D))
    return u, np.where(C_ > 0, S_ / np.maximum(C_, 1), np.nan).astype(np.float32), gate


def nearest(kd, D, gate, want, n=10):
    """Each wanted knockdown's n nearest other knockdowns by cosine of Delta over gated genes."""
    X = np.nan_to_num(D[:, gate]).astype(np.float64)
    X /= np.maximum(np.linalg.norm(X, axis=1, keepdims=True), 1e-12)
    pos = {k: i for i, k in enumerate(kd)}
    out = {}
    for k in want:
        if k not in pos: continue
        c = X @ X[pos[k]]; c[pos[k]] = -np.inf
        top = np.argsort(-c)[:n]
        out[k] = [(str(kd[j]), float(c[j])) for j in top]
    return out


def partners(a):
    lists = SS.knockdown_lists()
    need = set().union(*map(set, lists.values()))
    by = tables_by_stem()
    res, nb_all = {}, set()
    for src_table in ('replogle_k562_gw', 'replogle_rpe1'):
        kd, D, gate = full_delta(by[src_table])
        nb = nearest(kd, D, gate, sorted(need))
        res[src_table] = nb
        for v in nb.values(): nb_all |= {p for p, _ in v}
        print(src_table, len(kd), 'knockdowns,', len(nb), 'needed knockdowns held', flush=True)
    json.dump(res, open(OUT / 'tables/step11_partners.json', 'w'))
    defs = lists['DEF']; k5 = res['replogle_k562_gw']
    share = dict(def_targets=len(defs), def_in_k562_gw=sum(k in k5 for k in defs),
                 def_with_partner_cos_gt_0_5=sum(k in k5 and k5[k][0][1] > 0.5 for k in defs))
    (OUT / 'tables/step11_partner_share.json').write_text(json.dumps(share, indent=1) + '\n')
    print(share, flush=True)
    # Delta rows of the partner knockdowns that delta/ does not already hold, every table
    import h5py
    from anndata.io import read_elem
    extra = nb_all - need
    genes = list(pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str))
    (OUT / 'delta_nb').mkdir(exist_ok=True)
    for name, p in by.items():
        if not (OUT / 'delta' / f'{name}.npz').exists():
            continue                                                   # no usable Delta (Marson) or not needed
        base = pd.read_csv(p.with_name(name + '_baseline.csv'), index_col=0)
        ctrl = base['log1p_bulk_50k_full'].reindex(genes).to_numpy(np.float64)
        with h5py.File(p, 'r') as f:
            kd = read_elem(f['obs'])['target_gene'].astype(str).to_numpy()
            pick = sorted(i for i, k in enumerate(kd) if k in extra)
            if not pick: continue
            full = f['layers/log1p_bulk_50k_full'][pick]
        np.savez_compressed(OUT / 'delta_nb' / f'{name}.npz', knockdowns=kd[pick], genes=np.array(genes),
                            delta=(np.asarray(full, np.float64) - ctrl[None, :]).astype(np.float32))
        print(name, len(pick), 'partner rows', flush=True)


def stand_ins(filt, drop):
    """Metadata-only stand-ins of the sources a filter keeps, minus the fold's lines (no layers loaded)."""
    import h5py
    import stage1 as S
    from anndata.io import read_elem
    by = tables_by_stem(); genes = np.array(pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str), dtype=object)
    out, masks = [], {}
    for name, p in by.items():
        with h5py.File(p, 'r') as fh:
            obs = read_elem(fh['obs']); uns = {'effects_v5': dict(read_elem(fh['uns']['effects_v5']))} if 'effects_v5' in fh['uns'] else {}
        mod, ref, bg = S.table_meta(obs, uns, name)
        if not SS.table_kept(name, mod, ref, filt) or S._line_key(bg) in drop:
            continue
        masks[name] = SS.row_mask(p, filt)
        out.append(S.Source(name, bg, mod, ref, genes, np.zeros(0, dtype=object), np.zeros((0, 0)), np.zeros((0, 0)),
                            np.zeros(0), S.multiplier_for(mod, ref)))
    return out, masks


def delta_sources_nb(others, masks):
    """scorer_study.delta_sources over delta/ and delta_nb/ together (disjoint knockdowns), rows in table order."""
    import h5py
    import stage1 as S
    from anndata.io import read_elem
    by = tables_by_stem(); out = []
    for src in others:
        parts = [np.load(d / f'{src.name}.npz', allow_pickle=True) for d in (OUT / 'delta', OUT / 'delta_nb')
                 if (d / f'{src.name}.npz').exists()]
        if not parts: continue
        have = {}
        for z in parts:
            for k, row in zip(z['knockdowns'].astype(str), z['delta']):
                have.setdefault(k, []).append(row)
        with h5py.File(by[src.name], 'r') as fh:
            allk = read_elem(fh['obs'])['target_gene'].astype(str).to_numpy()
            rows = sorted(i for i, k in enumerate(allk) if k in have)
            se = fh['layers/se'][rows]
        kd = allk[rows]; cnt = {}
        D = []
        for k in kd:
            j = cnt.get(k, 0); D.append(have[k][j]); cnt[k] = j + 1
        D = np.stack(D)
        m = masks.get(src.name)
        if m is not None:
            j = np.array([r for r, i in enumerate(rows) if m[i]], dtype=int)
            if len(j) == 0: continue
            kd, D, se = kd[j], D[j], se[j]
        n = np.full(len(kd), np.nan, np.float32)
        if src.modality != 'drug':
            kd, D, se, n = S._collapse_duplicate_knockdowns(kd, D, se, n)
        out.append(S.Source(src.name, src.background, src.modality, src.reference, src.genes, kd,
                            np.asarray(D, np.float32), np.asarray(se, np.float32), n, src.mult))
    return out


class NbVotes:
    """The Delta BackgroundVotes with each partner's background votes appended at vote variance / w_nb."""
    def __init__(self, votes, partners, n, w):
        self.v, self.p, self.n, self.w = votes, partners, n, float(w)

    def gather(self, k, genes):
        Y, SE, M = self.v.gather(k, genes)
        Y, SE, M = list(Y), list(SE), list(M)
        for q, _ in self.p.get(str(k), [])[:self.n]:
            y2, s2, m2 = self.v.gather(q, genes)
            Y += list(y2); SE += [s / np.sqrt(self.w) for s in s2]; M += list(m2)
        return Y, SE, M


def pooled(votes, k, genes, prior):
    import stage1 as S
    Y, SE, M = votes.gather(k, genes)
    if not Y:
        return np.zeros(len(genes)), np.zeros(len(genes), np.int32), None
    mu, var_fe, tau2, _, n = S.pool_effects(np.stack(Y), np.stack(SE), np.asarray(M), prior['tau2'])
    t2 = np.where(np.isfinite(tau2), tau2, prior['tau2'])
    lfc, _, _, _ = S.gaussian_posterior(mu, var_fe, t2, prior['v0'])
    return lfc, n, np.stack(Y)


def pool(a):
    import stage1 as S
    from sklearn.utils.extmath import randomized_svd
    filt = a.filter; suf = '' if filt == 'F0' else f'_{filt}'
    P = json.load(open(OUT / 'tables/step11_partners.json'))
    genes = np.array(pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str), dtype=object)
    for f_, fl in SS.FOLD_LINES.items():
        drop = {SS.LINE_KEY[x] for x in fl}
        others, masks = stand_ins(filt, drop)
        dsrc = delta_sources_nb(others, masks); votes = S.BackgroundVotes(dsrc)
        part = P['replogle_rpe1'] if f_ == 'K562' else P['replogle_k562_gw']
        for line in fl:
            z = np.load(OUT / 'cache' / f'fold_{line}{suf}.npz', allow_pickle=True)
            prior = json.loads(str(z['prior'])); targets = list(z['targets']); K, G = len(targets), len(genes)
            gate = z['gate']; own = {i: list(genes).index(t) for i, t in enumerate(targets) if t in set(genes)}
            out = dict(sign_share=np.full((K, G), np.nan, np.float32), z_agree=np.full((K, G), np.nan, np.float32),
                       n_bg=np.zeros((K, G), np.int16), dpool_check=np.zeros((K, G), np.float32))
            for n in NB_N:
                for w in NB_W: out[f'dpool_nb{n}_w{w}'] = np.zeros((K, G), np.float32)
            for i, k in enumerate(targets):
                lfc, n_, Y = pooled(votes, k, genes, prior)
                out['dpool_check'][i] = lfc
                if Y is not None:
                    fin = np.isfinite(Y); nn = fin.sum(0); out['n_bg'][i] = nn
                    Yz = np.where(fin, Y, 0.0); mean = Yz.sum(0) / np.maximum(nn, 1)
                    sd = np.sqrt(np.where(fin, (Y - mean) ** 2, 0.0).sum(0) / np.maximum(nn - 1, 1))
                    sgn = np.sign(lfc)
                    out['sign_share'][i] = np.where((nn > 0) & (sgn != 0), (np.where(fin, np.sign(Y), 0) == sgn).sum(0) / np.maximum(nn, 1), np.nan)
                    out['z_agree'][i] = np.where(nn >= 2, np.abs(mean) / np.maximum(sd / np.sqrt(np.maximum(nn, 1)), 1e-12), np.nan)
                for n in NB_N:
                    for w in NB_W:
                        out[f'dpool_nb{n}_w{w}'][i] = pooled(NbVotes(votes, part, n, w), k, genes, prior)[0]
            ref_dp = np.where(z['dpool_n'] > 0, z['dpool'], 0.0); chk = np.where(z['dpool_n'] > 0, out['dpool_check'], 0.0)
            diff = float(np.abs(ref_dp - chk).max())
            # PCA: the fold sources' Delta rows (every knockdown they hold here), uncentred, on this line's gated genes
            M = np.concatenate([np.nan_to_num(s.lfc[:, gate]) for s in dsrc]).astype(np.float32)
            _, s, Vt = randomized_svd(M, N_PC, random_state=2026)
            np.savez_compressed(OUT / 'cache' / f's11_{line}{suf}.npz', **out, pc_Vt=Vt.astype(np.float32), pc_s=s,
                                pc_rows=M.shape[0], gate=gate, dpool_check_maxdiff=diff, partners=json.dumps(part))
            print(f'{filt} {f_} {line}: {K} knockdowns, Delta rows for PCA {M.shape[0]}, pool check max |diff| {diff:.2e}',
                  flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('cmd', choices=['partners', 'pool'])
    p.add_argument('--filter', default='F0', choices=['F0', 'F1', 'F2'])
    a = p.parse_args()
    {'partners': partners, 'pool': pool}[a.cmd](a)
