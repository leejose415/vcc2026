"""Configuration for the vcc_2026_4 transformer: YAML files in code/transformer/configs/, merged over DEFAULTS."""
import copy
from pathlib import Path

import yaml

V = Path.home() / 'vcc'

DEFAULTS = dict(
    name='dev',
    mode='validation',              # validation: targets = the six cell lines, no H1 anywhere; final: + H1
    gene_set='dev',                 # dev: our 300 targets + top genes up to n_genes; full: all 18,533 panel genes
    n_genes=5000,
    seed=2026,
    holdout_frac=0.2,               # share of target genes held out (validation rows), split by gene
    pie='auto',                     # on / off / auto (on when the PIE tables exist)
    gene_features='pie',            # pie: PIE feature projections + coverage flags; none: learned embedding only
    model=dict(d=128, layers=4, heads=4, ff=512, dropout=0.1, key_dim=16),
    train=dict(batch_size=32, lr=3e-4, weight_decay=0.01, warmup=500, epochs=20, patience=3, grad_clip=1.0,
               eff_noise=0.2, steps_per_epoch=None, bf16=True, eval_batches=None),
    loss=dict(huber_delta=1.0, se_eps=0.1, cos_weight=0.1, bce_weight=1.0, pos_weight_max=50.0),
    paths=dict(effects='lib/effects_v4', pie_features='lib/pie_features', cache='lib/transformer_cache',
               runs='lib/vcc_2026_4/runs'),
)


def _merge(a, b):
    out = copy.deepcopy(a)
    for k, v in (b or {}).items():
        out[k] = _merge(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else v
    return out


def load(path=None, **overrides):
    cfg = _merge(DEFAULTS, yaml.safe_load(Path(path).read_text()) if path else {})
    cfg = _merge(cfg, overrides)
    if cfg['mode'] not in ('validation', 'final'):
        raise ValueError(f"mode must be validation or final, got {cfg['mode']!r}")
    if cfg['gene_set'] not in ('dev', 'full'):
        raise ValueError(f"gene_set must be dev or full, got {cfg['gene_set']!r}")
    for k, v in cfg['paths'].items():
        cfg['paths'][k] = str(V / v) if not str(v).startswith('/') else str(v)
    return cfg
