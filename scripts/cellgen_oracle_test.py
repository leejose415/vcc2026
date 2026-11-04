"""Oracle and null tests of the cell generator (code/cellgen/cell_generator.py), scored like the harness.

Cells: Eval A's (lib/harness/evalA: raw unfiltered H1, thinned to the 20k regime). 10 random knockdowns (seed 2026) from
those both in Eval A and kept in lib/effects_v5/arc_h1.h5ad, with >= 200 Eval A cells; each one's cells split 50/50
(seeded). The scoring reference (results/cellgen_tests/ref/real.h5ad) = half 2 of every knockdown + Eval A's scoring
controls, the same control pool for every arm; its cache, ceiling and bundle (generic baseline = 0, replicate anchor = 1)
are built here with vcc3_generate.OfficialEvaluator, as the harness builds an eval's.
Truth from half 1, v5 definitions against the v5 H1 controls (what a v5 table row is): lfc = log2((m1 + 1e-9) /
(m0 + 1e-9)), m1 = half 1's mean CPM (panel library), m0 = arc_h1_baseline.csv control_mean_cpm; de_sig = BH < 0.05
over the gated genes of the rank test against arc_h1_control_quantiles.npz (the target gene included when it rejects).
Arms (each a 400-cell prediction per knockdown + Eval A's context controls, the control label cell-eval2 needs):
  oracle     generate_cells(true lfc, true de_sig, the v5 quantiles, the v5 control means)
  replicate  real half 1 itself (the ceiling a perfect model can hope for)
  null       generate_cells(lfc = 0, de_sig = False)
Pass: oracle scaled overall >= replicate scaled overall - 0.05.
--extra (after the first run; same reference, cache and bundle) adds diagnostic arms:
  oracle_nhalf     the oracle with n_cells = half 1's size (the truth's call set and the cells in one regime)
  oracle_nhalf_pb  the same, the mean target = half 1's pseudobulk (sum / total, the scorer's expression aggregate)
                   instead of its mean per-cell CPM
  control_paste    each knockdown = a seeded 400 of the context controls (the floor the null arm is compared with)
Outputs in results/cellgen_tests/ only (never lib/harness or lib/effects_v5).

Usage: ~/venv/bin/python scripts/cellgen_oracle_test.py [--threads 4]
"""
import argparse
import json
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import false_discovery_control

V = Path.home() / 'vcc'
sys.path.insert(0, str(V / 'scripts/harness')); sys.path.insert(0, str(V / 'code/cellgen'))
from common import CONTROL, G, eval_dir, panel_genes, read_rows, write_h5ad, log, save_json, seed_for   # noqa: E402
from score import raw_metrics                                                                            # noqa: E402
import cell_generator as C                                                                               # noqa: E402

OUT = V / 'results/cellgen_tests'
V5 = V / 'lib/effects_v5'
N_KD, MIN_CELLS, SEED, TOL = 10, 200, 2026, 0.05


def cpm_of(X):
    lib = np.asarray(X.sum(1)).ravel().astype(float)
    return sp.diags(1e6 / np.where(lib > 0, lib, 1)) @ X.astype(np.float64)


def truth(Xh, cq, m0, alpha=C.FDR):
    """v5 LFC and calls of real cells against the v5 controls (the quantiles and their mean)."""
    P = cpm_of(Xh); m1 = np.asarray(P.mean(0)).ravel()
    lfc = np.where(np.isfinite(m0), np.log2((m1 + C.EPS) / (np.nan_to_num(m0) + C.EPS)), np.nan)
    Pc = P.tocsc(); tested = cq.col
    p = np.array([C.rank_p(np.log1p(Pc[:, c].toarray().ravel()), cq.gene(i), cq.n0) for i, c in enumerate(tested)])
    sig = np.zeros(len(m0), bool); sig[tested[false_discovery_control(p, method='bh') < alpha]] = True
    return lfc, sig


def pred_file(path, blocks, ctx_X, genes):
    labels = [t for t, X in blocks for _ in range(X.shape[0])] + [CONTROL] * ctx_X.shape[0]
    X = sp.vstack([X for _, X in blocks] + [ctx_X]).tocsr()
    write_h5ad(path, X, pd.DataFrame(dict(target=labels, target_gene=labels), index=[f'c{i:07d}' for i in range(len(labels))]), genes)


def evaluator(out, reuse, threads):
    exe = Path(sys.executable).parent / 'cell-eval2'
    a = argparse.Namespace(cell_eval=str(exe) if exe.exists() else 'cell-eval2', eval_profile='vcc2026', num_threads=threads,
                           anchor_splits=5)
    return G.OfficialEvaluator(a, OUT / 'ref/real.h5ad', out, reuse=reuse)


def pseudobulk_cpm(X):
    tot = np.asarray(X.sum(0)).ravel().astype(float); return tot / tot.sum() * 1e6


def main(threads, extra=False):
    t0 = time.time(); genes = panel_genes(); n_panel = len(genes); d = eval_dir('A')
    real = ad.read_h5ad(d / 'real.h5ad', backed='r'); lab = real.obs['target_gene'].astype(str).to_numpy()
    assert real.var_names.astype(str).tolist() == genes; real.file.close()
    v5 = ad.read_h5ad(V5 / 'arc_h1.h5ad', backed='r'); kept = set(v5.obs.target_gene.astype(str)); v5.file.close()
    counts = pd.Series(lab[lab != CONTROL]).value_counts()
    eligible = sorted(t for t, n in counts.items() if t in kept and n >= MIN_CELLS)
    kds = sorted(np.random.default_rng(SEED).choice(eligible, N_KD, replace=False).tolist())
    log(f'{len(eligible)} eligible knockdowns (Eval A, v5-kept, >= {MIN_CELLS} cells); picked {kds}')
    halves = {}
    for t in kds:
        rows = np.random.default_rng(seed_for(SEED, 'cellgen_oracle', t)).permutation(np.flatnonzero(lab == t))
        halves[t] = (np.sort(rows[:len(rows) // 2]), np.sort(rows[len(rows) // 2:]))
    ctrl_rows = np.flatnonzero(lab == CONTROL)
    want = np.sort(np.concatenate([np.concatenate(h) for h in halves.values()] + [ctrl_rows]))
    M = read_rows(d / 'real.h5ad', want, np.arange(n_panel), n_panel); pos = {r: i for i, r in enumerate(want)}
    take = lambda rows: M[[pos[r] for r in rows]]
    # the scoring reference: half 2 + the scoring controls
    OUT.mkdir(parents=True, exist_ok=True)
    ref_labels = [t for t in kds for _ in halves[t][1]] + [CONTROL] * len(ctrl_rows)
    write_h5ad(OUT / 'ref/real.h5ad', sp.vstack([take(halves[t][1]) for t in kds] + [take(ctrl_rows)]).tocsr(),
               pd.DataFrame(dict(target=ref_labels, target_gene=ref_labels), index=[f'r{i:07d}' for i in range(len(ref_labels))]), genes)
    ctx = ad.read_h5ad(d / 'context.h5ad'); ctx_X = G.aligned_X(ctx, genes).astype(np.int32); del ctx
    # truth from half 1 and the three arms
    cq = C.ControlQuantiles(V5 / 'arc_h1_control_quantiles.npz')
    base = pd.read_csv(V5 / 'arc_h1_baseline.csv', index_col=0).reindex(genes); m0 = base.control_mean_cpm.to_numpy(float)
    names = ('oracle_nhalf', 'oracle_nhalf_pb', 'control_paste') if extra else ('oracle', 'replicate', 'null')
    arms = {k: [] for k in names}; reports = {}
    for i, t in enumerate(kds):
        X1 = take(halves[t][0]); lfc, sig = truth(X1, cq, m0)
        seed = int(seed_for(SEED, 'cellgen', t)) % (2 ** 32); n1 = X1.shape[0]; own = genes.index(t) if t in genes else None
        gen = lambda l, s_, n: C.generate_cells(l, s_, cq, m0, n_cells=n, library_size=20000, seed=seed, target=t, genes=genes)
        r = dict(half1_cells=int(n1), half2_cells=int(len(halves[t][1])), true_sig=int(sig.sum()),
                 true_sig_without_target=int(sig.sum() - (sig[own] if own is not None else 0)))
        if extra:
            lfc_pb = np.where(np.isfinite(lfc), np.log2((pseudobulk_cpm(X1) + C.EPS) / (np.nan_to_num(m0) + C.EPS)), np.nan)
            a1, a2 = gen(lfc, sig, n1), gen(lfc_pb, sig, n1)
            pick = np.random.default_rng(seed_for(SEED, 'cellgen_paste', t)).choice(ctx_X.shape[0], 400, replace=False)
            arms['oracle_nhalf'].append((t, a1.X.tocsr())); arms['oracle_nhalf_pb'].append((t, a2.X.tocsr()))
            arms['control_paste'].append((t, ctx_X[np.sort(pick)]))
            r.update(oracle_nhalf=a1.uns['cellgen'], oracle_nhalf_pb=a2.uns['cellgen']); first = a1.uns['cellgen']
        else:
            g_or, g_nu = gen(lfc, sig, 400), gen(np.zeros(n_panel), np.zeros(n_panel, bool), 400)
            arms['oracle'].append((t, g_or.X.tocsr())); arms['null'].append((t, g_nu.X.tocsr())); arms['replicate'].append((t, X1))
            r.update(oracle=g_or.uns['cellgen'], null=g_nu.uns['cellgen']); first = g_or.uns['cellgen']
        reports[t] = r
        log(f'  {i + 1}/{N_KD} {t}: half 1 {n1} cells, {int(sig.sum())} true calls -> realized {first["sig_realized"]} '
            f'(recall {first["sig_recall"]:.3f}, precision {first["sig_precision"]:.3f})')
    for k, blocks in arms.items(): pred_file(OUT / k / 'pred_cells.h5ad', blocks, ctx_X, genes)
    save_json(OUT / ('generation_extra.json' if extra else 'generation.json'), dict(knockdowns=kds, eligible=len(eligible), reports=reports))
    # score: the reference's cache, ceiling and bundle once (kept for --extra), then each arm against it
    if not extra:
        ev = evaluator(OUT / 'ref', None, threads); ev.prepare(); ev.build_bundle()
    res = {}
    for k in arms:
        out = OUT / k; e = evaluator(out, OUT / 'ref', threads)
        e.score(out / 'pred_cells.h5ad', out); scaled = e.scale(out)
        try: G.per_target_table(out).to_csv(out / 'per_target.csv')
        except Exception as ex: log(f'  per-target table skipped: {ex}')
        res[k] = dict(scaled=scaled, raw=raw_metrics(out)); log(f'{k}: scaled overall {scaled["overall"]:.4f}')
    if extra:
        prev = json.loads((OUT / 'summary.json').read_text())['arms']; res = {**prev, **res}
        rows = [dict(arm=k, overall=v['scaled']['overall'], **{m: v['scaled'].get(m) for m in G.VCC2026_METRICS}) for k, v in res.items()]
        pd.DataFrame(rows).to_csv(OUT / 'summary_all.csv', index=False)
        save_json(OUT / 'summary_extra.json', dict(arms=res, minutes=round((time.time() - t0) / 60, 1)))
        print(pd.DataFrame(rows).round(4).to_string(index=False), flush=True)
        return
    gap = res['replicate']['scaled']['overall'] - res['oracle']['scaled']['overall']
    summary = dict(knockdowns=kds, oracle_overall=res['oracle']['scaled']['overall'],
                   replicate_overall=res['replicate']['scaled']['overall'], null_overall=res['null']['scaled']['overall'],
                   gap_replicate_minus_oracle=gap, criterion=f'oracle >= replicate - {TOL}', passed=bool(gap <= TOL),
                   arms=res, minutes=round((time.time() - t0) / 60, 1))
    save_json(OUT / 'summary.json', summary)
    rows = [dict(arm=k, overall=v['scaled']['overall'], **{m: v['scaled'].get(m) for m in G.VCC2026_METRICS}) for k, v in res.items()]
    pd.DataFrame(rows).to_csv(OUT / 'summary.csv', index=False)
    print(json.dumps({k: v for k, v in summary.items() if k != 'arms'}, indent=1), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--threads', type=int, default=4); p.add_argument('--extra', action='store_true')
    a = p.parse_args(); main(a.threads, a.extra)
