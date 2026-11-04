"""Per-epoch evaluation for hyperparameter sweeps (validation mode only).

Three evaluation sets, the model in eval mode:
  train    a fixed random subset of training rows, as many as the held-out rows (seed split_seed + 1, drawn from the
           full training rows before any data-fraction subsetting, so it is identical across runs)
  heldout  the 20 % held-out knockdowns (same lines, unseen knockdowns; fixed by split_seed)
  h1       all 300 H1 knockdowns (an unseen line; H1 is never an input or a target in validation mode)
Metrics per set, over genes with loss_mask = 1 (measured, finite, not the knocked-down gene's own column, measured by
some source): per-knockdown R^2 = 1 - SS_res / SS_tot (SS_tot around that knockdown's own mean), mean and median;
DE R^2 (genes with DE label 1, knockdowns with >= 3 of them); per-knockdown Pearson r; the pooled calibration slope of
observed on predicted (1 = right size, > 1 = predictions too small); per-knockdown average precision of p_DE against the
DE labels (median and mean; "AUPRC").
Baselines, once per set: B1 = the precision-weighted mean of the evidence (Model 2; score |z|); B2 = the average knockdown
response: for train / heldout the mean lfc over that line's training knockdowns, for H1 the mean over H1's own 300
knockdowns ("reference only: uses H1 answers"); score |B2|.
"""
import numpy as np
import torch

from .data import LINE, LINES, SLOT, SOURCES
from .evaluate import h1_truth

METRICS = ['r2_mean', 'r2_median', 'de_r2_mean', 'de_r2_median', 'pearson_mean', 'pearson_median', 'slope',
           'auprc_median', 'auprc_mean', 'n_kd']


class Acc:
    """Accumulates per-knockdown statistics on the GPU, batch by batch."""

    def __init__(self):
        self.r2, self.de_r2, self.r, self.ap = [], [], [], []
        self.sx = self.sy = self.sxy = self.sxx = self.n = 0.0

    def add(self, pred, obs, mask, de, score):
        m = mask.float(); n = m.sum(1)
        ok = n >= 3
        mean_y = (obs * m).sum(1) / n.clamp(min=1); mean_x = (pred * m).sum(1) / n.clamp(min=1)
        ss_tot = (((obs - mean_y[:, None]) ** 2) * m).sum(1); ss_res = (((obs - pred) ** 2) * m).sum(1)
        r2 = 1 - ss_res / ss_tot.clamp(min=1e-12)
        cov = ((obs - mean_y[:, None]) * (pred - mean_x[:, None]) * m).sum(1)
        vx = (((pred - mean_x[:, None]) ** 2) * m).sum(1)
        r = cov / (vx * ss_tot).sqrt().clamp(min=1e-12)
        self.r2.append(r2[ok & (ss_tot > 0)]); self.r.append(r[ok & (vx > 0) & (ss_tot > 0)])
        d = (mask & de).float(); nd = d.sum(1)
        mdy = (obs * d).sum(1) / nd.clamp(min=1)
        dtot = (((obs - mdy[:, None]) ** 2) * d).sum(1); dres = (((obs - pred) ** 2) * d).sum(1)
        self.de_r2.append((1 - dres / dtot.clamp(min=1e-12))[(nd >= 3) & (dtot > 0)])
        self.sx += float((pred * m).sum()); self.sy += float((obs * m).sum())
        self.sxy += float((pred * obs * m).sum()); self.sxx += float((pred * pred * m).sum()); self.n += float(m.sum())
        # average precision per knockdown: unmasked genes sorted last with label 0, so they never count
        s = torch.where(mask, score, torch.full_like(score, -torch.inf))
        o = torch.argsort(s, dim=1, descending=True); lab = torch.gather((mask & de).float(), 1, o)
        tp = lab.cumsum(1); prec = tp / torch.arange(1, lab.shape[1] + 1, device=lab.device)[None]
        npos = lab.sum(1); ap = (prec * lab).sum(1) / npos.clamp(min=1)
        self.ap.append(ap[(npos > 0) & (npos < mask.sum(1))])

    def result(self):
        cat = lambda xs: torch.cat(xs).float().cpu().numpy() if xs else np.zeros(0)
        r2, de_r2, r, ap = cat(self.r2), cat(self.de_r2), cat(self.r), cat(self.ap)
        num = self.sxy - self.sx * self.sy / max(self.n, 1); den = self.sxx - self.sx ** 2 / max(self.n, 1)
        f = lambda a, fn: float(fn(a)) if len(a) else float('nan')
        return dict(r2_mean=f(r2, np.mean), r2_median=f(r2, np.median), de_r2_mean=f(de_r2, np.mean),
                    de_r2_median=f(de_r2, np.median), pearson_mean=f(r, np.mean), pearson_median=f(r, np.median),
                    slope=num / den if den > 0 else float('nan'), auprc_median=f(ap, np.median), auprc_mean=f(ap, np.mean),
                    n_kd=float(len(r2)))


class EvalSuite:
    def __init__(self, bank, cfg, bs=64):
        if bank.mode != 'validation': raise ValueError('sweep evaluation is validation mode only')
        self.bank, self.bs = bank, bs
        seed = cfg.get('split_seed', cfg['seed'])
        full_train = getattr(bank, 'train_rows_full', bank.train_rows)
        rng = np.random.default_rng(seed + 1)
        self.sets = dict(train=np.sort(rng.choice(full_train, len(bank.val_rows), replace=False)), heldout=bank.val_rows)
        obs, keys, L, de, measured = h1_truth(bank)
        self.h1_rows = bank.rows_for('H1', keys)
        dev = bank.device
        self.h1 = dict(obs=torch.as_tensor(L, device=dev), de=torch.as_tensor(de, device=dev),
                       mask=torch.as_tensor(measured, device=dev) & bank.has_data[None, :])
        self.h1['de'] &= self.h1['mask']
        # B2: average knockdown response per target line over its training knockdowns; H1 over its own knockdowns
        self.b2 = {}
        for line in bank.target_lines:
            s = SLOT[line]; tr = bank.row_t[np.intersect1d(np.flatnonzero(bank.row_line == LINE[line]), full_train)]
            X = bank.lfc[s][torch.as_tensor(tr, device=dev)].float()
            self.b2[LINE[line]] = torch.nanmean(X, 0).nan_to_num(0.0)
        Lh = torch.where(self.h1['mask'], self.h1['obs'], torch.nan)
        self.b2[LINE['H1']] = torch.nanmean(Lh, 0).nan_to_num(0.0)

    def _batches(self, name):
        bank = self.bank
        if name == 'h1':
            kd, line, trow = self.h1_rows
            for s in range(0, len(kd), self.bs):
                b = bank.batch(kd=kd[s:s + self.bs], line=line[s:s + self.bs], trow=trow[s:s + self.bs])
                sl = slice(s, s + self.bs)
                b.update(lfc_t=self.h1['obs'][sl], de=self.h1['de'][sl], loss_mask=self.h1['mask'][sl])
                yield b
        else:
            rows = self.sets[name]
            for s in range(0, len(rows), self.bs):
                yield bank.batch(rows=rows[s:s + self.bs])

    @torch.no_grad()
    def baselines(self):
        out = {}
        for name in ('train', 'heldout', 'h1'):
            a1, a2 = Acc(), Acc()
            for b in self._batches(name):
                a1.add(b['pw_mean'], b['lfc_t'], b['loss_mask'], b['de'], b['pw_z'].abs())
                m2 = torch.stack([self.b2[int(l)] for l in b['line']])
                a2.add(m2, b['lfc_t'], b['loss_mask'], b['de'], m2.abs())
            out[name] = dict(B1=a1.result(), B2=a2.result())
        return out

    @torch.no_grad()
    def evaluate(self, model, bf16=True):
        was = model.training; model.eval(); out = {}
        for name in ('train', 'heldout', 'h1'):
            acc = Acc()
            for b in self._batches(name):
                with torch.autocast('cuda', dtype=torch.bfloat16, enabled=bf16):
                    o = model(b)
                acc.add(o['delta'].float(), b['lfc_t'], b['loss_mask'], b['de'], torch.sigmoid(o['p_logit'].float()))
            out[name] = acc.result()
        model.train(was)
        return out


def trust_weights(model):
    """Per-source trust (the learned bias of the source attention) and its softmax over all slots."""
    t = model.evidence.trust.detach().float().cpu()
    names = [line for _, line in SOURCES] + ['PIE']
    return dict(raw={n: float(v) for n, v in zip(names, t)}, softmax={n: float(v) for n, v in zip(names, torch.softmax(t, 0))})
