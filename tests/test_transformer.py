"""Leakage and shape tests for the vcc_2026_4 transformer data loader and model, on the real effects_v4 tables
(tiny config: 500 genes). Run from ~/vcc: PYTHONPATH=code python -m pytest tests/test_transformer.py -q
Skipped when the effects_v4 tables are not on this machine."""
from pathlib import Path

import numpy as np
import pytest
import torch

V = Path.home() / 'vcc'
pytestmark = pytest.mark.skipif(not (V / 'lib/effects_v4/replogle_k562_gw.h5ad').exists(), reason='effects_v4 not here')

from transformer import config as C                                    # noqa: E402
from transformer.data import LINE, LINES, N_SLOTS, PIE_LEAKED, PIE_SLOT, SLOT, SLOT_FEATURES, Bank   # noqa: E402
from transformer.losses import losses                                  # noqa: E402
from transformer.model import VCCTransformer                           # noqa: E402

DEV = 'cuda' if torch.cuda.is_available() else 'cpu'


@pytest.fixture(scope='module')
def bank():
    return Bank(C.load(V / 'code/transformer/configs/tiny.yaml'), device=DEV, verbose=False)


@pytest.fixture(scope='module')
def batch(bank):
    rng = np.random.default_rng(0)
    rows = np.concatenate([rng.choice(np.flatnonzero(bank.row_line == LINE[l]), 4, replace=False)
                           for l in bank.target_lines])            # every target line in one batch
    return rows, bank.batch(rows=rows, train=True)


def test_validation_mode_has_no_h1(bank, batch):
    rows, b = batch
    assert bank.mode == 'validation' and 'H1' not in bank.loaded
    assert bank.lfc[SLOT['H1']] is None and bank.se[SLOT['H1']] is None
    assert not np.any(bank.row_line == LINE['H1'])
    assert torch.all(b['slot_feats'][:, :, SLOT['H1']] == 0)
    assert torch.all(torch.isnan(bank.eff[:, SLOT['H1']]))           # no H1 efficiency either


def test_target_line_never_in_its_own_evidence(bank, batch):
    rows, b = batch
    for i, r in enumerate(rows):
        s = SLOT[LINES[bank.row_line[r]]]
        assert torch.all(b['slot_feats'][i, :, s] == 0), LINES[bank.row_line[r]]
        assert not b['slot_mask'][i, :, s].any()


def test_efficiency_never_uses_target_line(bank):
    r = np.flatnonzero(bank.row_line == LINE['K562'])[:16]
    b = bank.batch(rows=r)
    e = bank.eff[bank._rows_t[0][torch.as_tensor(r, device=bank.device)]].clone()
    e[:, SLOT['K562']] = torch.nan
    want = torch.nanmedian(e, 1).values
    got = b['kd_feats'][:, -2]
    known = torch.isfinite(want)
    assert torch.allclose(got[known], want[known]) and torch.all(b['kd_feats'][:, -1] == known.float())


def test_no_pie_for_leaked_targets(bank):
    for l in PIE_LEAKED:
        assert int(bank.line_pie[LINE[l]]) == -1
    for l in ('HCT116', 'HEK293T'):
        r = np.flatnonzero(bank.row_line == LINE[l])[:8]
        assert torch.all(bank.batch(rows=r)['slot_feats'][:, :, PIE_SLOT] == 0)


def test_heldout_genes_never_training_targets(bank):
    train_keys = set(bank.kd_vocab[bank.row_kd[bank.train_rows]])
    val_keys = set(bank.kd_vocab[bank.row_kd[bank.val_rows]])
    assert not (train_keys & bank.heldout_keys)
    assert val_keys <= bank.heldout_keys and len(val_keys) > 0
    assert len(bank.train_rows) + len(bank.val_rows) == len(bank.row_kd)


def test_tensors_nan_free_and_shapes(bank, batch):
    rows, b = batch
    B, G = len(rows), len(bank.gcols)
    for k, v in b.items():
        if torch.is_floating_point(v): assert torch.isfinite(v).all(), k
    assert b['slot_feats'].shape == (B, G, N_SLOTS, SLOT_FEATURES)
    assert b['slot_mask'].shape == (B, G, N_SLOTS)
    assert b['basal'].shape == (B, G, 5) and b['kd_feats'].shape == (B, 7)
    for k in ('lfc_t', 'se_t', 'de', 'loss_mask', 'pw_mean'): assert b[k].shape == (B, G), k
    assert not (b['de'] & ~b['loss_mask']).any()
    own = b['gene_cols'][None, :] == b['kd_col'][:, None]
    assert not (b['loss_mask'] & own).any()


def test_de_labels_match_table_counts(bank):
    for line, (mine, table) in bank.n_de_bh.items():
        assert np.array_equal(mine, table), line


def test_model_forward_backward(bank, batch):
    cfg = bank.cfg; rows, b = batch
    model = VCCTransformer(cfg, len(bank.genes)).to(bank.device)
    with torch.autocast(DEV, dtype=torch.bfloat16, enabled=DEV == 'cuda'):
        out = model(b)
    assert out['delta'].shape == b['lfc_t'].shape and torch.isfinite(out['delta']).all()
    loss, parts = losses(out, b, cfg)
    loss.backward()
    assert torch.isfinite(loss) and model.gamma.grad is not None
