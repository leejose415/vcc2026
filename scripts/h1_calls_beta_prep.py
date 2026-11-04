"""Gather H1 (Eval A) arrays for the calls and beta figures: truth / C+G realised / C+G requested calls, true and
predicted scorer-space Delta. Writes one npz."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
V = Path.home() / 'vcc'; sys.path.insert(0, str(V / 'code/vcc6'))
import scorer_part2 as P2, scorer_arms as SA, scorer_study as SS, run_line as RL

SCR = Path(sys.argv[1])
ref = P2.refs_for('H1'); h = ref.h
genes = list(h.genes); gp = {g: j for j, g in enumerate(genes)}; tp = {t: i for i, t in enumerate(h.targets)}
E = V / 'results/stage1/2026-10-27_incumbent/evalA/CG'
inp = np.load(E / 'input.npz', allow_pickle=True)
assert list(inp['knockdowns']) == list(h.targets) and list(inp['genes']) == genes

def calls(f):
    d = pd.read_parquet(f, columns=['target', 'feature', 'p_adj', 'log2_fold_change'])
    T = np.zeros((h.T, h.G), bool); L = np.full((h.T, h.G), np.nan, np.float32)
    ii = d['target'].map(tp).to_numpy(); jj = d['feature'].map(gp).to_numpy()
    L[ii, jj] = d['log2_fold_change'].to_numpy(); s = d['p_adj'].to_numpy() < 0.05
    T[ii[s], jj[s]] = True
    return T, L

H = E / 'e2e/salt_none/harness/evalA'
truth, lfc_true = calls(H / 'de_real.parquet'); real, _ = calls(H / 'de_pred.parquet')
tested = ~np.isnan(lfc_true)
X = RL.read_controls(SS.EVAL_DIRS['H1'] / 'context.h5ad', genes)
m0, gate, _, _ = RL.control_profile(X, SCR / 'ctrl_H1')
own = np.zeros((h.T, h.G), bool)
for i in h.rows:
    if h.own[i] >= 0: own[i, h.own[i]] = True
for A in (truth, real): A[own] = False
req = inp['de_sig'] & (inp['n_sources'] > 0) & ~own
lfc = np.where(inp['n_sources'] > 0, inp['lfc'].astype(np.float64), 0.0)
P = SA.bulk(SA.mu_from_lfc(lfc, m0)); ctrl_bulk = np.log1p(SS.BTS * m0 / m0.sum())
chk = np.load(V / 'results/stage1/2026-10-28_scorer_space/tables/pred_cells_CG_H1.npz')
print('P matches step-9 proxy cache:', float(np.abs(P - chk['P']).max()))
ok = gate[None, :] & ~own
np.savez_compressed(SCR / 'h1_arrays.npz', targets=np.array(h.targets), genes=np.array(genes), truth=truth, real=real,
                    req=req, tested=tested & ~own, lfc_true=lfc_true, m0=m0, gate=gate, ok=ok,
                    pred_delta=(P - ctrl_bulk[None, :]).astype(np.float32),
                    true_delta=(h.R - h.ctrl[None, :]).astype(np.float32), n_sources=inp['n_sources'])
print('truth median', np.median(truth.sum(1)), 'real median', np.median(real.sum(1)), 'req median', np.median(req.sum(1)))
