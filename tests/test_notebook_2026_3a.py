"""Structure checks for vcc_2026_3a: vcc_2026_3's code and Δ tables, generated with the top-1,000 arm."""
import ast
from pathlib import Path
import subprocess

import nbformat

ROOT = Path(__file__).resolve().parents[1]


def test_notebook_structure_and_commands():
    from IPython.core.inputtransformer2 import TransformerManager
    notebook = nbformat.read(ROOT / 'vcc_2026_3a.ipynb', as_version=4)
    nbformat.validate(notebook)
    transform = TransformerManager()
    sources = [c.source for c in notebook.cells if c.cell_type == 'code']
    assert not any(s.startswith('%%writefile') for s in sources)                # reuses vcc_2026_3's scripts
    assert "CODE=str(V / 'code' / 'vcc_2026_3')" in sources[0] and 'TOP_N = 1000' in sources[0]
    assert "set(changed) != {'top_n'}" in sources[1]                              # only N may differ from submission 3
    for i, cell in enumerate(notebook.cells):
        if cell.cell_type != 'code': continue
        assert cell.execution_count is None and cell.outputs == []
        compile(transform.transform_cell(cell.source), f'cell_{i}', 'exec')
        for node in ast.walk(ast.parse(cell.source)):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id in ['run', 'bg']:
                value = node.args[-1]
                if isinstance(value, ast.Constant) and isinstance(value.value, str):
                    result = subprocess.run(['bash', '-n'], input=value.value, text=True, capture_output=True)
                    assert result.returncode == 0, result.stderr
