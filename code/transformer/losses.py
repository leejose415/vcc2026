"""Losses, averaged over the measured genes of the whole batch (loss_mask):
Huber (delta 1) on delta_hat vs lfc_t weighted 1 / (se_t^2 + se_eps^2); + cos_weight x (1 - cosine) per row over its
measured genes; + bce_weight x BCE on p_de vs the DE label, pos_weight = min(n_neg / n_pos, pos_weight_max)."""
import torch
import torch.nn.functional as F


def losses(out, b, cfg):
    lc = cfg['loss']; m = b['loss_mask']; mf = m.float(); n = mf.sum().clamp(min=1)
    w = 1.0 / (b['se_t'].square() + lc['se_eps'] ** 2)
    hub = F.huber_loss(out['delta'], b['lfc_t'], reduction='none', delta=lc['huber_delta'])
    huber = (w * hub * mf).sum() / n
    pred = out['delta'] * mf; true = b['lfc_t'] * mf
    rows = m.sum(1) >= 2
    cos = F.cosine_similarity(pred[rows], true[rows], dim=1, eps=1e-8)
    cosine = (1 - cos).mean() if rows.any() else out['delta'].sum() * 0
    y = b['de'].float(); n_pos = (y * mf).sum(); n_neg = n - n_pos
    pos_weight = torch.clamp(n_neg / n_pos.clamp(min=1), max=lc['pos_weight_max']).detach()
    bce_all = F.binary_cross_entropy_with_logits(out['p_logit'], y, reduction='none', pos_weight=pos_weight)
    bce = (bce_all * mf).sum() / n
    total = huber + lc['cos_weight'] * cosine + lc['bce_weight'] * bce
    return total, dict(loss=total.item(), huber=huber.item(), cosine=cosine.item(), bce=bce.item(),
                       pos_weight=float(pos_weight), measured=int(n.item()))
