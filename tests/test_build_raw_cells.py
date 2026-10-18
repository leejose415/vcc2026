"""Checks for scripts/build_raw_cells.py: labelling rules, the UMI cut and the streamed CSR files."""
import importlib
from pathlib import Path
import sys

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture()
def brc(tmp_path, monkeypatch):
    sys.path.insert(0, str(ROOT / 'scripts'))
    try:
        module = importlib.import_module('build_raw_cells')
        monkeypatch.setattr(module, 'OUT', tmp_path / 'cells')
        monkeypatch.setattr(module, 'CHUNK', 7)                       # several chunks on tiny data
        yield module
    finally:
        sys.path.remove(str(ROOT / 'scripts')); sys.modules.pop('build_raw_cells', None)


def test_labels(brc):
    assert brc.standard_target('non-targeting', {'non-targeting'}) == ('non-targeting', 0)
    assert brc.standard_target('TP53', {'non-targeting'}) == ('TP53', 1)
    assert brc.pb_target('CTRL') == ('non-targeting', 0) and brc.pb_target('None') is None
    assert brc.pb_target('TP53,RB1') == ('RB1+TP53', 2) and brc.pb_target('FDPS_2') == ('FDPS', 1)
    assert brc.pb_target('NegCtrl3,HMGCR') == ('HMGCR', 1) and brc.pb_target('NKX2-1') == ('NKX2-1', 1)


def test_replogle_layout_uses_library_umis(brc, tmp_path, monkeypatch):
    rng = np.random.default_rng(0); n, g = 30, 6
    X = rng.poisson(3, (n, g)).astype(np.float32)
    obs = pd.DataFrame(dict(gene=['non-targeting'] * 10 + ['A'] * 10 + ['B'] * 10, gene_id='x', sgID_AB='g', gem_group=1,
                            UMI_count=np.r_[np.full(5, 500.0), np.full(25, 5000.0)]), index=[f'c{i}' for i in range(n)])
    var = pd.DataFrame(dict(gene_name=['G0', 'G1', 'G1', 'G3', 'G4', 'G5']), index=[f'E{i}' for i in range(g)])
    path = tmp_path / 'replogle.h5ad'; ad.AnnData(X, obs=obs, var=var).write_h5ad(path)
    brc.DATASETS['toy'] = ('replogle', path, dict(cell_line='Toy', dataset='Test'))
    row = brc.build('toy')                                                      # default: no cell dropped
    a = ad.read_h5ad(brc.OUT / 'toy.h5ad')
    assert a.n_obs == 30 and row['dropped_low_umi'] == 0 and row['cells_below_1000_umi'] == 5
    assert np.array_equal(a.X.toarray(), X.astype(np.int32)) and a.obs.total_umi.tolist()[:5] == [500.0] * 5
    monkeypatch.setattr(brc, 'MIN_UMI', 1000)                                     # --min-umi 1000
    row = brc.build('toy')
    a = ad.read_h5ad(brc.OUT / 'toy.h5ad')
    assert a.n_obs == 25 and row['dropped_low_umi'] == 5 and (a.obs.total_umi >= 1000).all()   # the library total decides
    assert a.X.dtype == np.int32 and np.array_equal(a.X.toarray(), X[5:].astype(np.int32))
    assert list(a.var_names) == ['G0', 'G1', 'G1-1', 'G3', 'G4', 'G5'] and a.var.gene_symbol.tolist()[2] == 'G1'
    assert a.var.gene_id.tolist()[0] == 'E0' and set(a.obs.cell_line) == {'Toy'}
    assert (a.obs.n_genes_detected.to_numpy() == (X[5:] > 0).sum(1)).all()
    assert pd.read_csv(brc.OUT / 'manifest.csv').min_umi.tolist() == [1000]
    assert pd.read_csv(brc.OUT / 'manifest.csv').name.tolist() == ['toy']


def test_perturbase_layout_keeps_and_labels_every_cell(brc, tmp_path, monkeypatch):
    rng = np.random.default_rng(1); files = {}
    for cond in ('Dox', 'Mock'):
        X = sp.csr_matrix(rng.poisson(60, (12, 25)).astype(np.float32)); X[0] = 0
        labels = ['CTRL', 'None', 'TP53', 'TP53,RB1', 'FDPS_2', 'NegCtrl1'] * 2
        p = tmp_path / cond / 'raw.h5ad'; p.parent.mkdir()
        ad.AnnData(X, obs=pd.DataFrame(dict(gene=labels), index=[f'{cond}{i}' for i in range(12)]),
                   var=pd.DataFrame(index=[f'G{i}' for i in range(25)])).write_h5ad(p)
        files[cond] = p
    brc.DATASETS['toy_pb'] = ('perturbase', files, dict(cell_line='MCF10A', dataset='Test'))
    row = brc.build('toy_pb')
    a = ad.read_h5ad(brc.OUT / 'toy_pb.h5ad')
    # per condition: row 0 (CTRL) has 0 UMIs and rows 1 and 7 have no guide ('None'); all are kept and labelled
    assert a.n_obs == 24 and row['dropped_low_umi'] == 0 and row['unassigned_cells'] == 4 and row['cells_below_1000_umi'] == 2
    assert set(a.obs.condition) == {'Dox', 'Mock'}
    assert set(a.obs.target_gene) == {'non-targeting', 'TP53', 'RB1+TP53', 'FDPS', 'unassigned'}
    assert a.obs.loc[a.obs.target_gene == 'unassigned', 'n_targets'].eq(-1).all()
    assert a.obs.loc[a.obs.target_gene == 'RB1+TP53', 'n_targets'].eq(2).all() and row['pairs'] == 1 and row['knockdowns'] == 2
