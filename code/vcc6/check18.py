"""Step 18 integrity check (Jose: "nothing shifted / went missing?"). Read-only on Step 18's outputs.

  1. generator reports: 300 knockdowns per context in list order; requested calls = the input arrays' calls (as the
     generator sees them); n_unmeasured; mean-count error; library.
  2. full scans (chunked, backed) of every file: per-cell integer / non-negative / finite counts, library sizes, and the
     per-(context, knockdown) count sums -> pseudobulk.
     files: the three per-context generator outputs, the merged h5ad, the .vcc's own payload (decompressed to a temp
     file, deleted afterwards), and the uploaded C+G merged h5ad (2026-10-28).
  3. identity: per-knockdown count sums equal between generator output, merged file and .vcc payload (exact).
  4. fidelity: generated pseudobulk vs the requested means (input.npz as the generator sees it), per knockdown, in the
     scorer's space; and Delta(generated) vs Delta(requested).
  5. vs the uploaded C+G file: same obs layout; per-knockdown Delta correlation.
-> results/stage1/2026-10-30_step18/tables/check.json
"""
import json
import sys
import tarfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

V = Path.home() / 'vcc'
R = V / 'results/stage1/2026-10-30_step18'
CG = V / 'results/stage1/2026-10-28_DEF/final/prediction_DEF.h5ad'
TG = [t for t in (V / 'results/target_lists/final_DEF_300.txt').read_text().split() if t]
GENES = pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str).tolist()
CTRL = 'non-targeting'


def log(m):
    print(f'[{time.strftime("%H:%M:%S", time.gmtime())}] {m}', flush=True)


def scan(path, context=None, chunk=20000):
    """Per-(context, target) count sums and cell counts, plus per-cell sanity, from a backed h5ad."""
    import anndata as ad
    a = ad.read_h5ad(path, backed='r')
    assert list(a.var_names.astype(str)) == GENES, f'{path}: genes differ from the panel'
    tcol = 'target_gene' if 'target_gene' in a.obs else 'target'
    t = a.obs[tcol].astype(str).to_numpy()
    c = a.obs['context'].astype(str).to_numpy() if 'context' in a.obs else np.full(len(t), context)
    keys = pd.Index(pd.unique(pd.Series([f'{x}|{y}' for x, y in zip(c, t)])))
    code = keys.get_indexer([f'{x}|{y}' for x, y in zip(c, t)])
    S = np.zeros((len(keys), a.n_vars)); n = np.bincount(code, minlength=len(keys))
    lib = np.zeros(a.n_obs); bad = dict(non_integer=0, negative=0, non_finite=0); nnz = 0
    for s in range(0, a.n_obs, chunk):
        X = a.X[s:s + chunk]; X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
        d = X.data
        bad['non_integer'] += int((d != np.round(d)).sum()); bad['negative'] += int((d < 0).sum())
        bad['non_finite'] += int((~np.isfinite(d)).sum()); nnz += X.nnz
        lib[s:s + chunk] = np.asarray(X.sum(1)).ravel()
        cc = code[s:s + chunk]
        M = sp.csr_matrix((np.ones(len(cc)), (cc, np.arange(len(cc)))), shape=(len(keys), len(cc)))
        S += (M @ X).toarray()
    runs = int((np.diff(code) != 0).sum() + 1)                                   # contiguous blocks
    a.file.close()
    return dict(keys=keys, sums=S, n=n, lib=lib, bad=bad, nnz=nnz, blocks=runs, n_obs=len(t))


def log1p_pb(S):
    return np.log1p(5e4 * S / np.maximum(S.sum(1, keepdims=True), 1))


def main():
    out = {}
    # 1. generator reports
    rep = {}
    for x in 'DEF':
        g = pd.read_csv(R / f'config1p_K0/{x}/e2e/salt_none/prediction_generator.csv')
        z = np.load(R / f'config1p_K0/{x}/input.npz', allow_pickle=True)
        un = z['n_sources'] <= 0; sig = z['de_sig'] & ~un
        rep[x] = dict(rows=len(g), order_matches_list=list(g.knockdown) == TG,
                      requested_equals_input=bool((g.sig_requested.to_numpy() == sig.sum(1)).all()),
                      n_unmeasured_equals_input=bool((g.n_unmeasured.to_numpy() == un.sum(1)).all()),
                      mean_count_max_abs_error_max=float(g.mean_count_max_abs_error.max()),
                      realised_minus_requested=dict(min=int((g.sig_realized - g.sig_requested).min()),
                                                    median=float((g.sig_realized - g.sig_requested).median()),
                                                    max=int((g.sig_realized - g.sig_requested).max())),
                      recall_min=float(g.sig_recall.min()), library_median_range=[float(g.library_median.min()), float(g.library_median.max())],
                      seeds_unique=bool(g.seed.is_unique))
    out['generator_reports'] = rep; log(json.dumps(rep))
    # 2. scans
    scans = {}
    for x in 'DEF':
        scans[f'gen_{x}'] = scan(R / f'config1p_K0/{x}/e2e/salt_none/prediction.h5ad', context=x); log(f'scanned gen_{x}')
    scans['merged'] = scan(R / 'final/prediction_DEF.h5ad'); log('scanned merged')
    tmp = R / 'final/_vcc_payload_check'; tmp.mkdir(exist_ok=True)
    try:
        import zstandard as zstd
        with tarfile.open(R / 'final/prediction_DEF.vcc') as tf:
            members = tf.getnames(); out['vcc_members'] = members
            meta = json.load(tf.extractfile('meta.json')) if 'meta.json' in members else None
            out['vcc_meta'] = meta
            with tf.extractfile('pred.h5ad.zst') as src, open(tmp / 'pred.h5ad', 'wb') as dst:
                zstd.ZstdDecompressor().copy_stream(src, dst)
        scans['vcc'] = scan(tmp / 'pred.h5ad'); log('scanned vcc payload')
    finally:
        for f in tmp.glob('*'):
            f.unlink()
        tmp.rmdir()
    scans['cg_uploaded'] = scan(CG); log('scanned uploaded C+G')
    out['files'] = {k: dict(cells=v['n_obs'], groups=len(v['keys']), contiguous_blocks=v['blocks'], nnz=int(v['nnz']),
                            cells_per_group=sorted(set(v['n'].tolist()))[:5], **v['bad'],
                            library=dict(min=float(v['lib'].min()), median=float(np.median(v['lib'])), max=float(v['lib'].max())))
                    for k, v in scans.items()}
    # 3. identity of per-knockdown sums
    ident = {}
    m = scans['merged']; mi = {k: i for i, k in enumerate(m['keys'])}
    for x in 'DEF':
        gsc = scans[f'gen_{x}']; ks = [k for k in gsc['keys'] if not k.endswith('|' + CTRL)]
        gi = [list(gsc['keys']).index(k) for k in ks]
        ident[f'gen_{x}_vs_merged'] = dict(knockdowns=len(ks), all_present=all(k in mi for k in ks),
                                           max_abs_sum_diff=float(np.abs(gsc['sums'][gi] - m['sums'][[mi[k] for k in ks]]).max()),
                                           controls_in_merged=int(sum(k.endswith('|' + CTRL) for k in m['keys'])))
    v = scans['vcc']; vi = {k: i for i, k in enumerate(v['keys'])}
    ident['merged_vs_vcc'] = dict(same_groups=set(vi) == set(mi), same_order=list(v['keys']) == list(m['keys']),
                                  max_abs_sum_diff=float(np.abs(v['sums'][[vi[k] for k in m['keys']]] - m['sums']).max()),
                                  same_nnz=v['nnz'] == m['nnz'])
    out['identity'] = ident; log(json.dumps(ident))
    # 4. fidelity vs requested means, and 5. vs the uploaded C+G
    fid, cgc = {}, {}
    c = scans['cg_uploaded']; ci = {k: i for i, k in enumerate(c['keys'])}
    out['layout_vs_cg'] = dict(same_group_order=list(c['keys']) == list(m['keys']), same_cells=c['n_obs'] == m['n_obs'])
    for x in 'DEF':
        z = np.load(R / f'config1p_K0/{x}/input.npz', allow_pickle=True)
        un = z['n_sources'] <= 0; lfc = np.where(un, 0.0, z['lfc'].astype(float)); m0 = z['m0']; gate = z['gate']
        mu = np.maximum((m0[None, :] + 1e-9) * np.exp2(lfc) - 1e-9, 0.0)
        P_req = np.log1p(5e4 * mu / mu.sum(1, keepdims=True))
        gsc = scans[f'gen_{x}']; ks = list(gsc['keys'])
        P_gen = log1p_pb(gsc['sums'][[ks.index(f'{x}|{t}') for t in TG]])
        P_ctl = log1p_pb(gsc['sums'][[ks.index(f'{x}|{CTRL}')]])[0]
        P_cg = log1p_pb(c['sums'][[ci[f'{x}|{t}'] for t in TG]])
        d_req, d_gen, d_cg = P_req - P_ctl[None, :], P_gen - P_ctl[None, :], P_cg - P_ctl[None, :]
        r_lvl = [np.corrcoef(P_gen[i, gate], P_req[i, gate])[0, 1] for i in range(len(TG))]
        r_d = [np.corrcoef(d_gen[i, gate], d_req[i, gate])[0, 1] for i in range(len(TG))]
        r_cg = [np.corrcoef(d_gen[i, gate], d_cg[i, gate])[0, 1] for i in range(len(TG))]
        # a shifted gene axis would show up as the best correlation at a nonzero offset
        sh = {s: float(np.median([np.corrcoef(np.roll(d_gen[i], s)[gate], d_req[i, gate])[0, 1] for i in range(0, 300, 10)]))
              for s in (-1, 0, 1)}
        fid[x] = dict(level_r_min=float(np.min(r_lvl)), level_r_median=float(np.median(r_lvl)),
                      delta_r_min=float(np.min(r_d)), delta_r_median=float(np.median(r_d)),
                      abs_level_diff_median=float(np.median(np.abs(P_gen - P_req)[:, gate])),
                      abs_level_diff_p999=float(np.quantile(np.abs(P_gen - P_req)[:, gate], 0.999)),
                      knockdown_with_worst_delta_r=TG[int(np.argmin(r_d))], gene_shift_check_median_r=sh)
        cgc[x] = dict(delta_r_vs_cg_min=float(np.min(r_cg)), delta_r_vs_cg_median=float(np.median(r_cg)),
                      delta_norm_ratio_new_over_cg_median=float(np.median(np.linalg.norm(d_gen[:, gate], axis=1)
                                                                          / np.maximum(np.linalg.norm(d_cg[:, gate], axis=1), 1e-12))))
    out['fidelity_vs_requested'] = fid; out['vs_uploaded_cg'] = cgc
    log(json.dumps(dict(fid=fid, cg=cgc), default=float))
    (R / 'tables/check.json').write_text(json.dumps(out, indent=1, default=lambda o: o.tolist() if hasattr(o, 'tolist') else str(o)) + '\n')


if __name__ == '__main__':
    main()
