"""Submission checklist (Jose, 2026-11-02) on the live entry's cells: config 1' + K0 (Step 18's D/E/F generator outputs) and
the same config on Eval A (15c cells), with Eval A's real cells as truth. Read only; outputs in
results/stage1/2026-11-02_step20/checklist/.

Per file, one chunked pass: per-cell totals, per-group (knockdown / control) gene sums and sums of squares, a hash of
every cell's sparse row, and count shares (top-50 control genes, MT-, RPL/RPS).

  ~/venv/bin/python code/vcc6/checklist20.py stats   per-file statistics (cached .npz)
  ~/venv/bin/python code/vcc6/checklist20.py report  the eight checks
"""
import hashlib
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

V = Path.home() / 'vcc'
OUT = V / 'results/stage1/2026-11-02_step20/checklist'
S18 = V / 'results/stage1/2026-10-30_step18/config1p_K0'
FILES = {'D': S18 / 'D/e2e/salt_none/prediction.h5ad', 'E': S18 / 'E/e2e/salt_none/prediction.h5ad',
         'F': S18 / 'F/e2e/salt_none/prediction.h5ad',
         'EvalA_pred': V / 'results/stage1/2026-10-29_step15/c/config1p_K0/EvalA/e2e/salt_none/prediction.h5ad',
         'EvalA_real': V / 'lib/harness/evalA/real.h5ad'}
CTRL = 'non-targeting'


def labels(f):
    o = f['obs']
    for col in ('target_gene', 'target'):
        if col in o:
            g = o[col]
            if isinstance(g, h5py.Group):
                cats = g['categories'][:].astype(str); return cats[g['codes'][:]]
            return g[:].astype(str)
    raise KeyError('no target column')


def stats(name, chunk=4000):
    out = OUT / f'stats_{name}.npz'
    if out.exists():
        return dict(np.load(out, allow_pickle=True))
    with h5py.File(FILES[name], 'r') as f:
        lab = labels(f); X = f['X']; n, G = X.attrs['shape']
        ip = X['indptr'][:]
        groups, gi = np.unique(lab, return_inverse=True)
        S = np.zeros((len(groups), G)); Q = np.zeros((len(groups), G)); tot = np.zeros(n)
        hashes = np.empty(n, dtype=object)
        for r0 in range(0, n, chunk):
            r1 = min(n, r0 + chunk); a, b = ip[r0], ip[r1]
            d = X['data'][a:b]; ix = X['indices'][a:b]; p = ip[r0:r1 + 1] - a
            M = sp.csr_matrix((d.astype(np.float64), ix, p), shape=(r1 - r0, G))
            tot[r0:r1] = np.asarray(M.sum(1)).ravel()
            Z = sp.csr_matrix((np.ones(r1 - r0), (gi[r0:r1], np.arange(r1 - r0))), shape=(len(groups), r1 - r0))
            S += (Z @ M).toarray(); Q += (Z @ M.multiply(M)).toarray()
            for i in range(r1 - r0):
                hashes[r0 + i] = hashlib.blake2b(ix[p[i]:p[i + 1]].tobytes() + d[p[i]:p[i + 1]].tobytes(), digest_size=16).hexdigest()
            if r0 % 40000 == 0:
                print(name, r0, n, flush=True)
    res = dict(groups=groups, gi=gi, S=S, Q=Q, tot=tot, hashes=hashes.astype(str), n_per_group=np.bincount(gi))
    np.savez(out, **res)
    return res


def report():
    genes = pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str).to_numpy()
    mt = np.char.startswith(genes.astype(str), 'MT-'); rib = pd.Series(genes).str.match(r'^RP[LS]\d').to_numpy()
    R = {}
    st = {k: stats(k) for k in FILES}
    for name in ('D', 'E', 'F', 'EvalA_pred', 'EvalA_real'):
        s = st[name]; g = s['groups']; ci = int(np.flatnonzero(g == CTRL)[0]) if CTRL in g else None
        kd = np.flatnonzero(g != CTRL)
        # controls: the file's own (D/E/F, Eval A pred), or Eval A pred's controls for the real file (the same 18,400 cells?)
        ref = st[name] if ci is not None else st['EvalA_pred']; rci = int(np.flatnonzero(ref['groups'] == CTRL)[0])
        ctrl_sum = ref['S'][rci]; ctrl_n = ref['n_per_group'][rci]
        ctrl_tot = ref['tot'][ref['gi'] == rci]
        ctrl_cpm = 1e6 * ctrl_sum / ctrl_sum.sum()
        gate = ctrl_cpm > 5
        cell_kd = np.isin(s['gi'], kd)
        r = dict(cells=int(len(s['tot'])), knockdowns=int(len(kd)), cells_per_kd_median=float(np.median(s['n_per_group'][kd])),
                 control_cells=int(ctrl_n), gated_genes=int(gate.sum()))
        # 1. gated genes whose knockdown-group mean is exactly 0
        z = (s['S'][kd][:, gate] == 0).sum(1)
        r['c1_gated_genes_all_zero_per_kd'] = dict(median=float(np.median(z)), q90=float(np.quantile(z, 0.9)), max=int(z.max()),
                                                   kd_with_any=int((z > 0).sum()))
        # 3. count shares
        top50 = np.argsort(-ctrl_sum)[:50]
        kd_sum = s['S'][kd].sum(0)
        for nm, m in (('top50', np.isin(np.arange(len(genes)), top50)), ('MT', mt), ('RPL_RPS', rib)):
            r[f'c3_share_{nm}'] = dict(knockdown_cells=float(kd_sum[m].sum() / kd_sum.sum()), controls=float(ctrl_sum[m].sum() / ctrl_sum.sum()))
        # per-knockdown spread of the top-50 share
        sh = s['S'][kd][:, top50].sum(1) / s['S'][kd].sum(1)
        r['c3_top50_share_per_kd'] = dict(q05=float(np.quantile(sh, 0.05)), median=float(np.median(sh)), q95=float(np.quantile(sh, 0.95)))
        # 7, 8. depth and its spread
        t = s['tot'][cell_kd]
        r['c7_8_depth'] = dict(kd_cells=dict(median=float(np.median(t)), q05=float(np.quantile(t, 0.05)), q95=float(np.quantile(t, 0.95)),
                                             cv=float(t.std() / t.mean()), sd_log=float(np.log(t).std())),
                               controls=dict(median=float(np.median(ctrl_tot)), q05=float(np.quantile(ctrl_tot, 0.05)),
                                             q95=float(np.quantile(ctrl_tot, 0.95)), cv=float(ctrl_tot.std() / ctrl_tot.mean()),
                                             sd_log=float(np.log(ctrl_tot).std())))
        # 6. duplicate cells
        h = pd.Series(s['hashes'][cell_kd]); dup = h.duplicated(keep=False)
        hk = pd.DataFrame(dict(h=s['hashes'][cell_kd], g=s['gi'][cell_kd]))
        across = hk[hk.h.isin(hk.h[hk.h.duplicated()])].groupby('h').g.nunique()
        ctrl_h = set(ref['hashes'][ref['gi'] == rci])
        r['c6_duplicates'] = dict(kd_cells_in_a_duplicate=int(dup.sum()), hashes_shared_across_knockdowns=int((across > 1).sum()),
                                  kd_cells_identical_to_a_control=int(h.isin(ctrl_h).sum()) if name != 'EvalA_real' else None)
        # 5. per-gene dispersion within knockdown groups vs the control cells, at matched mean (Fano by mean bin)
        n = s['n_per_group'][kd][:, None]; mu = s['S'][kd] / n; var = s['Q'][kd] / n - mu ** 2
        cm = ctrl_sum / ctrl_n; cv_ = ref['Q'][rci] / ctrl_n - cm ** 2
        bins = [0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30]
        fano = []
        for lo, hi in zip(bins[:-1], bins[1:]):
            m = (mu >= lo) & (mu < hi) & gate[None, :]; mc = (cm >= lo) & (cm < hi) & gate
            if m.sum() > 50 and mc.sum() > 10:
                fano.append(dict(mean_bin=f'{lo}-{hi}', kd_fano_median=float(np.median(var[m] / mu[m])), ctrl_fano_median=float(np.median(cv_[mc] / cm[mc]))))
        r['c5_fano_by_mean'] = fano
        R[name] = r
        # 2. unmeasured genes (no source) in D/E/F: prediction mean vs the control mean
        if name in ('D', 'E', 'F'):
            z_ = np.load(S18 / f'{name}/input.npz', allow_pickle=True)
            order = [list(z_['knockdowns'].astype(str)).index(x) for x in g[kd]]
            un = (z_['n_sources'][order] <= 0) & gate[None, :]
            ratio = (mu[un] + 1e-9) / (cm[np.nonzero(un)[1]] + 1e-9)
            lr = np.log2(ratio)
            r['c2_unmeasured_gated'] = dict(pairs=int(un.sum()), share_of_gated_pairs=float(un.sum() / (len(kd) * gate.sum())),
                                            log2_pred_over_ctrl_median=float(np.median(lr)), abs_log2_q90=float(np.quantile(np.abs(lr), 0.9)),
                                            measured_pairs_abs_log2_median=float(np.median(np.abs(np.log2((mu[~un & gate[None, :]] + 1e-9) / (np.broadcast_to(cm, mu.shape)[~un & gate[None, :]] + 1e-9))))))
            r['c2_genes_unmeasured_for_every_kd'] = int(((z_['n_sources'] <= 0).all(0) & gate).sum())
    # 4. significant genes per knockdown (Eval A: cell-eval's DE tables), ours vs the reference
    H = V / 'results/stage1/2026-10-29_step15/c/config1p_K0/EvalA/e2e/salt_none/harness/evalA'
    cnt = {}
    for side in ('pred', 'real'):
        d = pd.read_parquet(H / f'de_{side}.parquet', columns=['target', 'p_adj'])
        cnt[side] = d[d.p_adj < 0.05].groupby('target').size()
    c = pd.DataFrame(cnt).fillna(0)
    R['c4_EvalA_sig_genes_per_kd'] = dict(pred=dict(median=float(c.pred.median()), q90=float(c.pred.quantile(0.9))),
                                          real=dict(median=float(c.real.median()), q90=float(c.real.quantile(0.9))),
                                          spearman=float(c.pred.rank().corr(c.real.rank())), kd=int(len(c)))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'checklist.json').write_text(json.dumps(R, indent=1, default=float) + '\n')
    print(json.dumps(R, indent=1, default=float))


if __name__ == '__main__':
    OUT.mkdir(parents=True, exist_ok=True)
    if sys.argv[1] == 'stats':
        for k in FILES:
            stats(k); print('done', k, flush=True)
    else:
        report()
