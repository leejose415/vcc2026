"""Behavioral checks for vcc_2026_3b: Model 2's sign arrays and top-N calling with gate_from."""
import argparse
import ast
import importlib
from pathlib import Path
import subprocess
import sys

import nbformat
import numpy as np
import pytest
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
NAMES = ['vcc2_common', 'vcc2_meta', 'vcc3_ridge', 'vcc3_generate', 'vcc3_screen']


def embedded(notebook):
    return {c.source.split('\n', 1)[0].split('/')[-1]: c.source.split('\n', 1)[1]
            for c in nbformat.read(ROOT / notebook, as_version=4).cells if c.source.startswith('%%writefile')}


@pytest.fixture(scope='module')
def scripts_v3b(tmp_path_factory):
    directory = tmp_path_factory.mktemp('v3b_scripts')
    for name, body in embedded('vcc_2026_3b.ipynb').items():
        (directory / name).write_text(body)
    for name in NAMES: sys.modules.pop(name, None)
    sys.path.insert(0, str(directory))
    try:
        yield {name: importlib.import_module(name) for name in NAMES}
    finally:
        sys.path.remove(str(directory))
        for name in NAMES: sys.modules.pop(name, None)


def test_notebook_structure_and_commands(scripts_v3b):
    from IPython.core.inputtransformer2 import TransformerManager
    notebook = nbformat.read(ROOT / 'vcc_2026_3b.ipynb', as_version=4)
    nbformat.validate(notebook)
    transform = TransformerManager()
    written = [c.source.split('\n', 1)[0] for c in notebook.cells if c.source.startswith('%%writefile')]
    assert written == [f'%%writefile ~/vcc/code/vcc_2026_3b/{n}.py' for n in NAMES[:2] + ['vcc3_ridge', 'vcc3_generate', 'vcc3_screen']]
    for i, cell in enumerate(notebook.cells):
        if cell.cell_type != 'code': continue
        assert cell.execution_count is None and cell.outputs == []
        compile(transform.transform_cell(cell.source), f'cell_{i}', 'exec')
        if not cell.source.startswith('%%writefile'):
            for node in ast.walk(ast.parse(cell.source)):
                if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id in ['run', 'bg']:
                    value = node.args[-1]
                    if isinstance(value, ast.Constant) and isinstance(value.value, str):
                        result = subprocess.run(['bash', '-n'], input=value.value, text=True, capture_output=True)
                        assert result.returncode == 0, result.stderr
    for name in ['vcc3_ridge', 'vcc3_generate', 'vcc3_screen', 'vcc2_meta']:
        result = subprocess.run([sys.executable, scripts_v3b[name].__file__, '--help'], capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
    # Only the three intended scripts differ from vcc_2026_3, and the schemas stay, so vcc_2026_3's fits load.
    v3, v3b = embedded('vcc_2026_3.ipynb'), embedded('vcc_2026_3b.ipynb')
    assert v3b['vcc2_common.py'] == v3['vcc2_common.py'] and v3b['vcc3_screen.py'] == v3['vcc3_screen.py']
    for name in ['vcc2_meta.py', 'vcc3_ridge.py', 'vcc3_generate.py']:
        assert v3b[name] != v3[name]
    assert scripts_v3b['vcc2_meta'].SCHEMA == 'vcc_2026_2_meta_v1' and scripts_v3b['vcc3_generate'].SCHEMA == 'vcc_2026_3_ridge_v1'
    headings = [c.source.split('\n', 1)[0] for c in notebook.cells if c.cell_type == 'markdown' and c.source.startswith('## ')]
    numbers = [h.split('.')[0] for h in headings]
    assert len(numbers) == len(set(numbers)), headings


def test_combine_sign_agreement_keeps_re_and_se(scripts_v3b):
    meta = scripts_v3b['vcc2_meta']
    lfc = np.array([[[1.0, -1.0, np.nan]], [[0.5, 2.0, np.nan]], [[-0.1, np.nan, np.nan]]])   # 3 sources x 1 row x 3 genes
    se = np.full_like(lfc, 0.2)
    re, re_se, k, agree = meta.combine(lfc, se, np.zeros(3))
    w = 1 / 0.04
    np.testing.assert_allclose(re[0, :2], [(1.0 + 0.5 - 0.1) / 3, (-1.0 + 2.0) / 2], rtol=1e-6)
    np.testing.assert_allclose(re_se[0, :2], [np.sqrt(1 / (3 * w)), np.sqrt(1 / (2 * w))], rtol=1e-6)
    np.testing.assert_array_equal(k[0], [3, 2, 0])
    np.testing.assert_allclose(agree[0, :2], [2 / 3, 1 / 2])                  # gene 0: two of three positive; gene 1: one of two
    assert np.isnan(re[0, 2]) and np.isnan(agree[0, 2])


def test_top_n_gate_from_ranks_by_gate_and_keeps_values(scripts_v3b):
    c, g = scripts_v3b['vcc2_common'], scripts_v3b['vcc3_generate']
    line = g.Line(sp.csr_matrix([[100, 100, 100, 100, 100]] * 2), None, ['KD', 'G1', 'G2', 'G3', 'G4'])
    line.gate = c.expression_gate([10, 10, 10, 10, 10])
    D = np.array([[-2.0, 0.3, 0.9, -0.4, 0.8]])                                # corrected values: G2 and G4 largest
    info = dict(gate_delta=np.array([[-2.0, 1.0, 0.1, -0.9, 0.2]]), gate_se=np.full((1, 5), 0.1),
                n_bg=np.array([[3.0, 3.0, 3.0, 3.0, 1.0]]), sign_agree=np.array([[1.0, 1.0, 1.0, 2 / 3, np.nan]]))
    top2 = dict(call='topn', top_n=2, alpha=1.0, gate_from='re')
    # Model 2 ranks G1 (z 10) and G3 (z 9) highest; the called genes keep the corrected values; the own gene never counts.
    np.testing.assert_allclose(g.calibrate(D, line, ['KD'], g.Cfg(argparse.Namespace(**top2)), info=info), [[0, 0.3, 0, -0.4, 0]])
    # The sign rule drops G3 (3 lines, not all agreeing); G4 is measured in one line only, so it stays a candidate.
    np.testing.assert_allclose(g.calibrate(D, line, ['KD'], g.Cfg(argparse.Namespace(sign_rule=True, **top2)), info=info),
                               [[0, 0.3, 0, 0, 0.8]])
    # Without gate_from the array's own z ranks (here supplied as re_z would be).
    own = dict(z=np.abs(D) / 0.1, n_bg=info['n_bg'], sign_agree=info['sign_agree'])
    np.testing.assert_allclose(g.calibrate(D, line, ['KD'], g.Cfg(argparse.Namespace(call='topn', top_n=2)), info=own), [[0, 0, 0.9, 0, 0.8]])
    with pytest.raises(ValueError, match='invalid generator settings'):
        g.Cfg(argparse.Namespace(gate_from='re', z=0.0))                         # the z call still needs z* > 0 to gate


def test_min_bg_drops_one_line_genes_and_keeps_old_settings(scripts_v3b):
    c, g = scripts_v3b['vcc2_common'], scripts_v3b['vcc3_generate']
    line = g.Line(sp.csr_matrix([[100, 100, 100, 100]] * 2), None, ['KD', 'G1', 'G2', 'G3'])
    line.gate = c.expression_gate([10, 10, 10, 10])
    D = np.array([[-2.0, 0.9, 0.5, -0.3]])
    info = dict(z=np.abs(D) / 0.1, n_bg=np.array([[3.0, 1.0, 2.0, np.nan]]), sign_agree=np.array([[1.0, 1.0, 0.5, np.nan]]))
    cfg = g.Cfg(argparse.Namespace(call='topn', top_n=2, min_bg=2))
    # G1 is measured in one line only and G3 in none: the best remaining candidate is G2, even though its lines disagree.
    np.testing.assert_allclose(g.calibrate(D, line, ['KD'], cfg, info=info), [[0, 0, 0.5, 0]])
    with pytest.raises(ValueError, match='min_bg needs n_bg'):
        g.calibrate(D, line, ['KD'], cfg, info=dict(z=info['z']))
    with pytest.raises(ValueError, match='invalid generator settings'):
        g.Cfg(argparse.Namespace(min_bg=0))
    assert 'min_bg' not in g.Cfg(argparse.Namespace()).as_dict() and cfg.as_dict()['min_bg'] == 2
    assert 'min_bg' in g.ARM_KEYS


def test_ridge_carries_sign_arrays(scripts_v3b):
    source = embedded('vcc_2026_3b.ipynb')['vcc3_ridge.py']
    assert "for name in ['re_z', 'n_bg', 'sign_agree']:" in source
