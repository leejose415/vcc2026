"""The vcc_2026_4 transformer: one token per gene plus a knockdown token, evidence pooled per gene by source attention.

Gene embedding e_g: per PIE feature source (esm2, string_space, depmap_gene_effect, gene_text) Linear(d_src -> d), summed
over the sources that cover the gene, plus a Linear of the coverage flags (gene_features='pie'); or a learned embedding
only (gene_features='none'). Evidence encoder: MLP per (gene, source) slot -> d, plus a source embedding. Source attention:
score = learned per-source trust b_s + <q(e_g), k_s>, softmax over available slots, weighted sum; a learned "no evidence"
vector when no slot is available. Gene token = e_g + MLP(basal) + pooled evidence; knockdown token = e_k + MLP(knockdown
features) + type embedding, prepended. TransformerEncoder (pre-norm, batch-first) on bf16 autocast.
Heads: delta_hat = Linear(d, 1) + gamma * precision-weighted mean evidence lfc (gamma learned, starts at 1);
p_de logit = Linear(d, 1).
"""
import json
from pathlib import Path

import numpy as np
import torch
from torch import nn

from .data import BASAL_FEATURES, KD_FEATURES, N_SLOTS, SLOT_FEATURES

PIE_SOURCES = ['esm2', 'string_space', 'depmap_gene_effect', 'gene_text']


def mlp(i, d, o=None):
    return nn.Sequential(nn.Linear(i, d), nn.GELU(), nn.Linear(d, o or d))


class GeneEmbedding(nn.Module):
    def __init__(self, n_genes, d, features='pie', feature_dir=None):
        super().__init__()
        self.features = features
        self.unknown = nn.Parameter(torch.zeros(d))            # knockdown genes off the panel
        if features == 'pie':
            meta = json.loads((Path(feature_dir) / 'sources.json').read_text())
            assert meta['columns'] == PIE_SOURCES, meta['columns']
            mask = np.load(Path(feature_dir) / 'coverage_mask.npy')
            self.register_buffer('cover', torch.as_tensor(mask, dtype=torch.float32), persistent=False)
            self.proj = nn.ModuleList()
            for k, name in enumerate(PIE_SOURCES):
                X = torch.as_tensor(np.load(Path(feature_dir) / f'{name}.npy'))      # float16 (18,533 x d_src)
                self.register_buffer(f'feat_{k}', X, persistent=False)
                self.proj.append(nn.Linear(X.shape[1], d))
            self.flags = nn.Linear(len(PIE_SOURCES), d)
        elif features == 'none':
            self.table = nn.Embedding(n_genes, d)
        else:
            raise ValueError(features)

    def forward(self, cols):
        """cols: panel column per token (-1 = not on the panel)."""
        c = cols.clamp(min=0)
        if self.features == 'pie':
            cover = self.cover[c]
            e = self.flags(cover)
            for k, p in enumerate(self.proj):
                e = e + p(getattr(self, f'feat_{k}')[c].float()) * cover[..., k:k + 1]
        else:
            e = self.table(c)
        return torch.where((cols >= 0)[..., None], e, self.unknown)


class EvidencePool(nn.Module):
    def __init__(self, d, key_dim):
        super().__init__()
        self.enc = mlp(SLOT_FEATURES, d)
        self.src = nn.Embedding(N_SLOTS, d)
        self.trust = nn.Parameter(torch.zeros(N_SLOTS))
        self.keys = nn.Parameter(torch.randn(N_SLOTS, key_dim) * 0.02)
        self.query = nn.Linear(d, key_dim)
        self.none = nn.Parameter(torch.zeros(d))

    def forward(self, feats, mask, e_g):
        """feats (B, G, S, F), mask (B, G, S), e_g (G, d) -> pooled (B, G, d), attention (B, G, S)."""
        h = self.enc(feats) + self.src.weight                       # (B, G, S, d)
        score = self.trust + (self.query(e_g) @ self.keys.T)[None]  # (B?, G, S)
        score = score.expand(mask.shape).masked_fill(~mask, float('-inf'))
        anyv = mask.any(-1, keepdim=True)
        a = torch.softmax(torch.where(anyv, score, 0.0), -1) * mask
        pooled = (a[..., None] * h).sum(-2)
        return torch.where(anyv, pooled, self.none), a


class VCCTransformer(nn.Module):
    def __init__(self, cfg, n_genes):
        super().__init__()
        m = cfg['model']; d = m['d']
        self.genes = GeneEmbedding(n_genes, d, cfg['gene_features'], cfg['paths']['pie_features'])
        self.evidence = EvidencePool(d, m['key_dim'])
        self.basal = mlp(BASAL_FEATURES, d)
        self.kd = mlp(KD_FEATURES, d)
        self.kd_type = nn.Parameter(torch.zeros(d))
        layer = nn.TransformerEncoderLayer(d, m['heads'], m['ff'], m['dropout'], activation='gelu', batch_first=True,
                                           norm_first=True)
        self.encoder = nn.TransformerEncoder(layer, m['layers'], enable_nested_tensor=False)
        self.norm = nn.LayerNorm(d)
        self.head_delta = nn.Linear(d, 1)
        self.head_p = nn.Linear(d, 1)
        self.gamma = nn.Parameter(torch.ones(()))

    def forward(self, b):
        e_g = self.genes(b['gene_cols'])                             # (G, d)
        pooled, attn = self.evidence(b['slot_feats'], b['slot_mask'], e_g)
        tokens = e_g[None] + self.basal(b['basal']) + pooled          # (B, G, d)
        kd_tok = self.genes(b['kd_col']) + self.kd(b['kd_feats']) + self.kd_type
        x = torch.cat([kd_tok[:, None], tokens], 1)
        x = self.norm(self.encoder(x))[:, 1:]
        delta = self.head_delta(x).squeeze(-1).float() + self.gamma * b['pw_mean']
        p_logit = self.head_p(x).squeeze(-1).float()
        return dict(delta=delta, p_logit=p_logit, attn=attn)
