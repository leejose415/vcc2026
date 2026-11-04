"""Checks for vcc_2026_3c: notebook structure and the generator's multi-arm, scaled oracle mode."""
import argparse
import ast
import importlib
import json
from pathlib import Path
import subprocess
import sys

import nbformat
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
NAMES = ['vcc2_common', 'vcc2_meta', 'vcc3_ridge', 'vcc3_generate', 'vcc3_screen']


def embedded(notebook):
    return {c.source.split('\n', 1)[0].split('/')[-1]: c.source.split('\n', 1)[1]
            for c in nbformat.read(ROOT / notebook, as_version=4).cells if c.source.startswith('%%writefile')}


@pytest.fixture(scope='module')
def scripts_v3c(tmp_path_factory):
    directory = tmp_path_factory.mktemp('v3c_scripts')
    for name, body in embedded('vcc_2026_3c.ipynb').items():
        (directory / name).write_text(body)
    for name in NAMES: sys.modules.pop(name, None)
    sys.path.insert(0, str(directory))
    try:
        yield {name: importlib.import_module(name) for name in NAMES}
    finally:
        sys.path.remove(str(directory))
        for name in NAMES: sys.modules.pop(name, None)


def test_notebook_structure_and_commands(scripts_v3c):
    from IPython.core.inputtransformer2 import TransformerManager
    notebook = nbformat.read(ROOT / 'vcc_2026_3c.ipynb', as_version=4)
    nbformat.validate(notebook)
    transform = TransformerManager()
    written = [c.source.split('\n', 1)[0] for c in notebook.cells if c.source.startswith('%%writefile')]
    assert written == [f'%%writefile ~/vcc/code/vcc_2026_3c/{n}.py' for n in NAMES]
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
        result = subprocess.run([sys.executable, scripts_v3c[name].__file__, '--help'], capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
    oracle_help = subprocess.run([sys.executable, scripts_v3c['vcc3_generate'].__file__, 'oracle', '--help'], capture_output=True, text=True)
    assert '--scale' in oracle_help.stdout and '--variants' in oracle_help.stdout
    # Only the generator differs from vcc_2026_3b, and the schema stays, so vcc_2026_3's predictions load.
    v3b, v3c = embedded('vcc_2026_3b.ipynb'), embedded('vcc_2026_3c.ipynb')
    assert [n for n in v3b if v3b[n] != v3c[n]] == ['vcc3_generate.py']
    assert scripts_v3c['vcc3_generate'].SCHEMA == 'vcc_2026_3_ridge_v1'
    headings = [c.source.split('\n', 1)[0] for c in notebook.cells if c.cell_type == 'markdown' and c.source.startswith('## ')]
    numbers = [h.split('.')[0] for h in headings]
    assert len(numbers) == len(set(numbers)), headings


def oracle_args(tmp_path, variants, **extra):
    genes = tmp_path / 'genes.csv'; pd.DataFrame({'gene_name': ['KD', 'G1']}).to_csv(genes, index=False)
    path = tmp_path / 'variants.json'; path.write_text(json.dumps(variants))
    values = dict(mode='oracle', genes=str(genes), out=str(tmp_path / 'out'), variants=str(path), reference=None, resume=False)
    return argparse.Namespace(**(values | extra))


def test_oracle_arms_are_oracle_or_nochange(scripts_v3c, tmp_path):
    g = scripts_v3c['vcc3_generate']
    with pytest.raises(ValueError, match='oracle arms use method oracle or nochange'):
        g.cmd_benchmark(oracle_args(tmp_path, [dict(method='oracle'), dict(method='re')]))


def test_scaled_oracle_end_to_end(scripts_v3c, tmp_path):
    """Real CPU evaluator on a small synthetic H1: oracle arms and true no change on one split, leaderboard-scaled."""
    import anndata as ad
    import numpy as np
    import scipy.sparse as sp
    binary = Path(sys.executable).parent / 'cell-eval2'
    if not binary.exists(): pytest.skip('official cell-eval2 is not installed')
    g = scripts_v3c['vcc3_generate']; genes = np.array(['K0', 'K1', 'K2'] + [f'G{i}' for i in range(37)])
    rng = np.random.default_rng(2); baseline = rng.integers(30, 100, len(genes))
    controls = rng.poisson(baseline, (160, len(genes))); groups = {}
    for i, kd in enumerate(genes[:3]):
        factors = np.ones(len(genes)); factors[3 + i::3] = 3; factors[i] = 0.2
        groups[kd] = rng.poisson(baseline * factors, (120, len(genes)))
    labels = ['non-targeting'] * len(controls) + [t for t, x in groups.items() for _ in range(len(x))]
    ad.AnnData(sp.csr_matrix(np.vstack([controls] + list(groups.values()))),
               obs=pd.DataFrame({'target_gene': labels}, index=[str(i) for i in range(len(labels))]),
               var=pd.DataFrame(index=genes)).write_h5ad(tmp_path / 'h1.h5ad')
    pd.DataFrame({'gene_name': genes}).to_csv(tmp_path / 'genes.csv', index=False)
    pd.DataFrame({'target_gene': list(groups)}).to_csv(tmp_path / 'targets.csv', index=False)
    arms = [dict(method='oracle', alpha=1.0, cell_rule='matched', label='oracle_matched'),
            dict(method='oracle', alpha=0.5, cell_rule='fixed', label='oracle_half'),
            dict(method='nochange', alpha=1.0, cell_rule='fixed', label='nochange')]
    (tmp_path / 'arms.json').write_text(json.dumps(arms))
    out = tmp_path / 'oracle'
    result = subprocess.run([sys.executable, g.__file__, 'oracle', '--cells', str(tmp_path / 'h1.h5ad'), '--genes', str(tmp_path / 'genes.csv'),
                             '--targets', str(tmp_path / 'targets.csv'), '--cell-eval', str(binary), '--num-threads', '1', '--n-cells', '60',
                             '--min-cells', '50', '--anchor-splits', '2', '--scale', '--variants', str(tmp_path / 'arms.json'), '--out', str(out)],
                            capture_output=True, text=True, timeout=900)
    assert result.returncode == 0, result.stdout[-3000:] + result.stderr[-3000:]
    scaled = pd.read_csv(out / 'scaled_scores.csv', index_col='arm')
    assert set(scaled.index) == {a['label'] for a in arms} and scaled['overall'].notna().all()
    assert (out / 'bundle/manifest.json').exists() and not (out / 'selected_generator.json').exists()
    # The true effects tell the knockdowns apart; no change cannot.
    assert scaled.loc['oracle_matched', 'pds_cosine'] > scaled.loc['nochange', 'pds_cosine']
    assert scaled.loc['oracle_matched', 'overall'] > scaled.loc['nochange', 'overall']


def test_oracle_needs_its_own_split(scripts_v3c, tmp_path):
    g = scripts_v3c['vcc3_generate']
    arms = [dict(method='oracle', cell_rule='matched'), dict(method='nochange')]
    for extra in (dict(reference=str(tmp_path)), dict(resume=True)):
        with pytest.raises(ValueError, match='no --reference or --resume'):
            g.cmd_benchmark(oracle_args(tmp_path, arms, **extra))
