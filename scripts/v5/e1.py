"""E1 loss diagnostics on dev_nopie (no training): gradient-norm share of the three loss terms over one training epoch,
and the share of the (weighted) Huber loss that comes from residuals beyond delta, by tercile of real effect size,
on the training lines (all their rows) and on H1 (zero-shot, 300 knockdowns).

  cd ~/vcc && PYTHONPATH=code:~/v5 ~/venv/bin/python ~/v5/e1.py
Effect size of a row = ||lfc_t - line mean of lfc_t|| over its measured genes, noise-corrected (minus sum se^2);
terciles within each line.
"""
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F

from transformer.data import LINES, kd_keys
from transformer.evaluate import h1_truth, load_run

V = Path.home() / 'vcc'
RUN = V / 'lib/vcc_2026_4/runs/dev_nopie'
OUT = V / 'results/tuning'


def terms(out, b, lc):
    m = b['loss_mask']; mf = m.float(); n = mf.sum().clamp(min=1)
    w = 1.0 / (b['se_t'].square() + lc['se_eps'] ** 2)
    hub = F.huber_loss(out['delta'], b['lfc_t'], reduction='none', delta=lc['huber_delta'])
    huber = (w * hub * mf).sum() / n
    pred = out['delta'] * mf; true = b['lfc_t'] * mf; rows = m.sum(1) >= 2
    cosine = (1 - F.cosine_similarity(pred[rows], true[rows], dim=1, eps=1e-8)).mean()
    y = b['de'].float(); n_pos = (y * mf).sum(); n_neg = n - n_pos
    pw = torch.clamp(n_neg / n_pos.clamp(min=1), max=lc['pos_weight_max']).detach()
    bce = (F.binary_cross_entropy_with_logits(out['p_logit'], y, reduction='none', pos_weight=pw) * mf).sum() / n
    return dict(huber=huber, cosine=lc['cos_weight'] * cosine, bce=lc['bce_weight'] * bce)


def grad_shares(cfg, bank, model):
    model.train(); params = [p for p in model.parameters() if p.requires_grad]
    lc = cfg['loss']; bs = cfg['train']['batch_size']
    order = np.random.default_rng(cfg['seed']).permutation(bank.train_rows); recs = []
    for s in range(0, len(order) - 1, bs):
        b = bank.batch(rows=order[s:s + bs], train=True)
        with torch.autocast('cuda', dtype=torch.bfloat16):
            out = model(b)
        t = terms(out, b, lc); norms = {}
        for k, v in t.items():
            g = torch.autograd.grad(v, params, retain_graph=True, allow_unused=True)
            norms[k] = float(torch.sqrt(sum((x.float() ** 2).sum() for x in g if x is not None)))
        tot = sum(norms.values()); recs.append(dict(**{f'norm_{k}': v for k, v in norms.items()},
                                                    **{f'share_{k}': v / tot for k, v in norms.items()},
                                                    **{f'loss_{k}': float(v) for k, v in t.items()}))
        if len(recs) % 200 == 0: print(f'  grad shares: {len(recs)} batches', flush=True)
    df = pd.DataFrame(recs); model.eval()
    return df.mean().to_dict(), len(df)


def huber_parts(pred, lfc, se, mask, delta, se_eps):
    r = pred - lfc; a = np.abs(r); w = 1 / (se ** 2 + se_eps ** 2)
    hub = np.where(a <= delta, 0.5 * r ** 2, delta * (a - 0.5 * delta)) * w * mask
    beyond = (a > delta) & mask
    return hub.sum(1), (hub * beyond).sum(1), beyond.sum(1), mask.sum(1)


@torch.no_grad()
def training_lines(cfg, bank, model):
    rows = np.arange(len(bank.row_kd)); bs = 64; lc = cfg['loss']
    # pass 1: per-line mean target lfc per gene (over measured entries)
    G = len(bank.gcols); s1 = {}; n = {}
    for s in range(0, len(rows), bs):
        b = bank.batch(rows=rows[s:s + bs]); ln = b['line'].cpu().numpy(); L = b['lfc_t'].cpu().numpy(); M = b['loss_mask'].cpu().numpy()
        for l in np.unique(ln):
            sel = ln == l; s1[l] = s1.get(l, 0) + (L[sel] * M[sel]).sum(0); n[l] = n.get(l, 0) + M[sel].sum(0)
    mean = {l: s1[l] / np.maximum(n[l], 1) for l in s1}
    recs = []
    for s in range(0, len(rows), bs):
        b = bank.batch(rows=rows[s:s + bs])
        with torch.autocast('cuda', dtype=torch.bfloat16):
            out = model(b)
        P = out['delta'].float().cpu().numpy(); L = b['lfc_t'].cpu().numpy(); S = b['se_t'].cpu().numpy()
        M = b['loss_mask'].cpu().numpy(); ln = b['line'].cpu().numpy()
        tot, beyond, nb, nm = huber_parts(P, L, S, M, lc['huber_delta'], lc['se_eps'])
        dm = (L - np.stack([mean[l] for l in ln])) * M
        size = np.sqrt(np.maximum((dm ** 2).sum(1) - ((S ** 2) * M).sum(1), 0))
        for j in range(len(P)):
            recs.append(dict(row=int(rows[s + j]), line=LINES[ln[j]], held_out=bool(rows[s + j] in held), size=size[j],
                             huber=tot[j], huber_beyond=beyond[j], n_beyond=nb[j], n=nm[j]))
    return pd.DataFrame(recs)


def h1_rows(cfg, bank):
    obs, keys, L, de, measured = h1_truth(bank)
    with h5py.File(Path(cfg['paths']['effects']) / 'arc_h1.h5ad', 'r') as f:
        S = f['layers/se'][:][:, bank.gcols]
    S = np.where(np.isfinite(S), S, 1.0)
    z = np.load(RUN / 'h1_pred_g3/raw_predictions.npz'); P = z['delta']
    assert np.array_equal(z['gene_cols'], bank.gcols)
    M = measured.astype(float); mean = (L * M).sum(0) / np.maximum(M.sum(0), 1)
    dm = (L - mean) * M; size = np.sqrt(np.maximum((dm ** 2).sum(1) - ((S ** 2) * M).sum(1), 0))
    tot, beyond, nb, nm = huber_parts(P, L, S, measured, cfg['loss']['huber_delta'], cfg['loss']['se_eps'])
    return pd.DataFrame(dict(target=obs['target_gene'].astype(str).to_numpy(), line='H1', size=size, huber=tot,
                             huber_beyond=beyond, n_beyond=nb, n=nm))


def tercile_table(df):
    df = df.copy(); df['tercile'] = df.groupby('line')['size'].transform(lambda x: np.argsort(np.argsort(x.to_numpy())) * 3 // len(x))
    df['tercile'] = df.tercile.map({0: 'weak', 1: 'mid', 2: 'strong'})
    g = df.groupby('tercile').agg(rows=('size', 'size'), median_size=('size', 'median'), huber=('huber', 'sum'),
                                  huber_beyond=('huber_beyond', 'sum'), n_beyond=('n_beyond', 'sum'), n=('n', 'sum'))
    g['share_of_huber_loss_beyond_delta'] = g.huber_beyond / g.huber; g['share_of_entries_beyond_delta'] = g.n_beyond / g.n
    return g.loc[['weak', 'mid', 'strong']]


if __name__ == '__main__':
    cfg, bank, model = load_run(RUN, device='cuda', ckpt='best.pt')
    held = set(bank.val_rows.tolist())
    res = {}
    tr = training_lines(cfg, bank, model); tr.to_csv(OUT / 'v5_E1_train_rows.csv', index=False)
    h1 = h1_rows(cfg, bank); h1.to_csv(OUT / 'v5_E1_h1_rows.csv', index=False)
    pd.set_option('display.width', 220)
    for nm, df in (('training_lines_all_rows', tr), ('training_lines_heldout_genes', tr[tr.held_out]), ('H1_300', h1)):
        t = tercile_table(df); res[f'huber_beyond_{nm}'] = t.reset_index().to_dict('records'); print(nm); print(t.round(4).to_string(), flush=True)
    shares, nb = grad_shares(cfg, bank, model)
    res['grad_shares_one_epoch'] = dict(batches=nb, **shares); print(json.dumps(res['grad_shares_one_epoch'], indent=1))
    (OUT / 'v5_E1_gpu.json').write_text(json.dumps(res, indent=1, default=float))
