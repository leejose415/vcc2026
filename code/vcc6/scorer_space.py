"""Scorer-space MSE (Jose 2026-10-28, results/stage1/2026-10-28_scorer_space/): cell_eval2 0.18.0's
expr_mse_unbiased_capped_norm rebuilt from its own primitives so every term can be read, split and swept.

Per scored knockdown p (G_p = the panel minus p's own gene where it resolves; w_p = 1 / |G_p|):
  num_p = ( ||mu_pred,p - mu_real,p||^2_{G_p} - r * min(C_pred,p, C_real,p) - C_real,p ) / |G_p|
  den_p = ( ||mu_real,p - mu_ctrl||^2_{G_p} - C_real,p - C_ctrl ) / |G_p|
  r     = min(1, B / sum_q w_q min(C_pred,q, C_real,q)),  B = _across_pert_budget(mu_pred, own genes out, w)
  raw   = sum_p num_p / sum_p den_p
mu = bulk_lognorm pseudobulk (log1p(5e4 x group count sum / its total)); C = delete-1 jackknife in that space
(moments.jackknife_correction). The real side (mu_real, mu_ctrl, C_real, C_ctrl) comes from the eval's harness cache
(real_cache/*.moments.npz, built by cell-eval2 itself); the predicted side from the cells via
prep.pseudobulk_bulk_lognorm_with_moments + moments.correction_for, or from a model's means directly (Steps 3-5).
Citations: metrics/delta.py:959 (formula), :43 / :699 (per-row cap), :157 / :858 (budget), :276 (own gene), :1122
(denominator); moments.py:160; prep.py:329-361; catalog.py:345.

  check  ~/venv/bin/python code/vcc6/scorer_space.py check --eval DIR --pred PRED.h5ad --harness SUMMARY.json
"""
import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

V = Path.home() / 'vcc'
CTRL, BTS = 'non-targeting', 50000.0
MSE = 'expr_mse_unbiased_capped_norm'


class Ref:
    """One harness eval's real side, in the scorer's space, plus its MSE anchors."""

    def __init__(self, eval_dir):
        d = Path(eval_dir)
        z = np.load(glob.glob(str(d / 'real_cache' / '*.moments.npz'))[0], allow_pickle=True)
        perts = z['perts'].astype(str); idx = {p: i for i, p in enumerate(perts)}
        self.genes = pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str).tolist()
        gp = {g: j for j, g in enumerate(self.genes)}
        self.targets = sorted(p for p in perts if p != CTRL)
        rows = [idx[t] for t in self.targets]
        self.R = z['means'][rows].astype(np.float64); self.Creal = z['jk'][rows].astype(np.float64)
        self.ctrl = z['means'][idx[CTRL]].astype(np.float64); self.Cctrl = float(z['jk'][idx[CTRL]])
        self.n_real = z['counts'][rows].astype(np.int64)
        self.own = np.array([gp.get(t, -1) for t in self.targets])
        G = len(self.genes); self.nrow = np.where(self.own >= 0, G - 1, G).astype(np.float64)
        self.excl = {i: int(c) for i, c in enumerate(self.own) if c >= 0}
        dr = self.R - self.ctrl[None]
        sq = np.einsum('ij,ij->i', dr, dr) - self.own_sq(self.R, np.broadcast_to(self.ctrl, self.R.shape))
        self.den = (sq - self.Creal - self.Cctrl) / self.nrow
        base = pd.read_csv(d / 'bundle' / 'baseline_agg.csv').set_index('statistic').loc['mean']
        anc = pd.read_parquet(d / 'bundle' / 'anchor_agg.parquet').set_index('metric')['replicate']
        self.b_mse, self.r_mse = float(base[MSE]), float(anc[MSE])

    def own_sq(self, A, B):
        out = np.zeros(A.shape[0])
        for i, c in self.excl.items():
            out[i] = (A[i, c] - B[i, c]) ** 2
        return out

    def parts(self, P, Cpred):
        """Every term per knockdown for predicted means P (T x G, rows = self.targets) and their jackknife Cpred
        (0 for a noise-free prediction). Returns a dict of arrays and the totals."""
        from cell_eval2.metrics.delta import PRED_TRACE_CAP_K, _across_pert_budget
        P = np.asarray(P, np.float64); Cpred = np.asarray(Cpred, np.float64)
        d = P - self.R
        plug = np.einsum('ij,ij->i', d, d) - self.own_sq(P, self.R)
        capped = np.minimum(Cpred, PRED_TRACE_CAP_K * self.Creal)
        w = 1.0 / self.nrow
        claim = float((capped * w).sum())
        B = _across_pert_budget(P, self.excl, w)
        r = min(1.0, B / claim) if claim > 0 else 1.0
        credit = r * capped
        num = (plug - credit - self.Creal) / self.nrow
        raw = float(num.sum() / self.den.sum())
        return dict(plug=plug, credit=credit, capped=capped, num=num, den=self.den, raw=raw, budget=B, claim=claim, r=r,
                    budget_binds=bool(claim > 0 and B < claim), scaled_unclamped=self.scaled(raw),
                    scaled=min(1.0, max(0.0, self.scaled(raw))))

    def scaled(self, raw):
        return (raw - self.b_mse) / (self.r_mse - self.b_mse)


def pred_side(pred_h5ad, targets):
    """Predicted pseudobulk means and jackknife C from generated cells, rows in `targets` order."""
    import anndata as ad
    from cell_eval2.moments import correction_for
    from cell_eval2.prep import pseudobulk_bulk_lognorm_with_moments
    a = ad.read_h5ad(pred_h5ad)
    col = 'target_gene' if 'target_gene' in a.obs else 'target'
    perts, means, mom = pseudobulk_bulk_lognorm_with_moments(a, col, bulk_target_sum=BTS)
    C = correction_for(mom, perts, means, comparator='bulk_lognorm')
    idx = {str(p): i for i, p in enumerate(perts)}
    rows = [idx[t] for t in targets]
    return np.asarray(means)[rows].astype(np.float64), np.asarray(C)[rows].astype(np.float64)


def check(a):
    ref = Ref(a.eval)
    P, C = pred_side(a.pred, ref.targets)
    out = ref.parts(P, C)
    harness = json.loads(Path(a.harness).read_text())['raw'][MSE]
    res = dict(eval=a.eval, pred=a.pred, rebuilt_raw=out['raw'], harness_raw=harness, abs_diff=abs(out['raw'] - harness),
               reproduces_1e6=bool(abs(out['raw'] - harness) < 1e-6), budget=out['budget'], claim=out['claim'], r=out['r'],
               budget_binds=out['budget_binds'], rows_capped=int((C > ref.Creal).sum()), knockdowns=len(ref.targets),
               baseline_raw=ref.b_mse, replicate_raw=ref.r_mse, scaled_unclamped=out['scaled_unclamped'])
    print(json.dumps(res, indent=1))
    if a.out:
        Path(a.out).write_text(json.dumps(res, indent=1) + '\n')
    sys.exit(0 if res['reproduces_1e6'] else 1)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)
    c = sub.add_parser('check'); c.add_argument('--eval', required=True); c.add_argument('--pred', required=True)
    c.add_argument('--harness', required=True); c.add_argument('--out')
    a = p.parse_args()
    {'check': check}[a.cmd](a)
