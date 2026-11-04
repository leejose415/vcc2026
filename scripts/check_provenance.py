"""Provenance checks for the vcc_2026_4 inputs. Read only; exits 1 if any check fails.

Per table in lib/effects_v4/ (every <source>.h5ad):
  1. uns['effects_v4'] exists, with se_floor 0.05 and fdr 0.05;
  2. built from lib/qc_cells/ (cell sources: uns cell_file) or, for the table sources (Orion pseudobulk, Marson CD4 T),
     from the source table named in uns (lib/orion/, lib/tcell/), as the v4 builder documents;
  3. on-target threshold -0.5 and MIN_CELLS 30: read from uns when recorded there, otherwise from perturbation_qc.json
     next to the knockdown_qc.csv that uns names (the builder records the QC files, the QC records its parameters);
  4. rows = kept knockdowns (knockdown_keep) in that knockdown_qc.csv.
Then: the old effects path lib/effects/ must not appear in code/transformer/, the v4 scripts or vcc_2026_4.ipynb
(other mentions are listed for information). Writes results/qc/provenance_check.json.
"""
import ast
import json
import re
import sys
from pathlib import Path

import h5py
import pandas as pd
from anndata.io import read_elem

V = Path.home() / 'vcc'
E = V / 'lib/effects_v4'
WANT = dict(se_floor=0.05, fdr=0.05, MAX_ONTARGET_LFC=-0.5, MIN_CELLS=30)
STRICT = ['code/transformer', 'scripts/build_effects_v4.py', 'scripts/check_provenance.py', 'scripts/build_baseline_v2.py',
          'scripts/k562_benchmark.py', 'scripts/export_pie_features.py', 'vcc_2026_4.ipynb']


def local(p):
    """uns paths are absolute on the VM that built them; resolve them inside this ~/vcc."""
    p = str(p); i = p.find('/vcc/')
    return V / p[i + 5:] if i >= 0 else Path(p)


def check_table(path):
    with h5py.File(path, 'r') as f:
        u = read_elem(f['uns']).get('effects_v4', {}); n_rows = f['layers/lfc'].shape[0]
    r = dict(source=path.stem, rows=n_rows, fails=[])
    if not u:
        r['fails'].append('no uns effects_v4'); return r
    for k in ('se_floor', 'fdr'):
        r[k] = float(u.get(k, 'nan'))
        if r[k] != WANT[k]: r['fails'].append(f'{k}={u.get(k)}')
    if 'cell_file' in u:
        r['built_from'] = str(u['cell_file'])
        if '/lib/qc_cells/' not in r['built_from']: r['fails'].append('cell_file not in lib/qc_cells/')
    elif 'table' in u:
        r['built_from'] = str(u['table'])
        if not re.search(r'/lib/(orion|tcell)/', r['built_from']): r['fails'].append('table source outside lib/orion, lib/tcell')
    else:
        r['fails'].append('no cell_file or table in uns')
    kq = local(u.get('knockdown_qc', ''))
    params = {}
    for k in ('MAX_ONTARGET_LFC', 'MIN_CELLS'):
        if k in u or k.lower() in u:
            params[k] = (float(u.get(k, u.get(k.lower()))), 'uns')
    pj = kq.parent / 'perturbation_qc.json'
    if len(params) < 2:
        if not pj.exists():
            r['fails'].append(f'no {pj}')
        else:
            p = json.loads(pj.read_text()).get('parameters', {})
            for k in ('MAX_ONTARGET_LFC', 'MIN_CELLS'):
                if k not in params and k in p: params[k] = (float(p[k]), 'perturbation_qc.json')
    for k in ('MAX_ONTARGET_LFC', 'MIN_CELLS'):
        if k not in params: r['fails'].append(f'{k} not recorded')
        else:
            r[k], r[f'{k}_from'] = params[k]
            if params[k][0] != WANT[k]: r['fails'].append(f'{k}={params[k][0]}')
    if kq.exists():
        kept = int(pd.read_csv(kq)['knockdown_keep'].astype(bool).sum()); r['kept_in_qc'] = kept
        if kept != n_rows: r['fails'].append(f'rows {n_rows} != kept {kept}')
    else:
        r['fails'].append(f'missing {kq}')
    r['ok'] = not r['fails']
    return r


def code_strings(src):
    """(line, text) of every string constant in Python source, docstrings excluded (comments are not tokens); a bare
    'effects' counts only when joined onto a path with / (a dict key named 'effects' is not a path)."""
    try:
        tree = ast.parse(src)
    except SyntaxError:
        return [(i, l) for i, l in enumerate(src.splitlines(), 1) if not l.lstrip().startswith(('#', '%', '!'))]
    docs = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) and node.body:
            first = node.body[0]
            if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant) and isinstance(first.value.value, str):
                docs.add(id(first.value))
    joined = {id(n.right) for n in ast.walk(tree) if isinstance(n, ast.BinOp) and isinstance(n.op, ast.Div)}
    return [(n.lineno, n.value if not (n.value == 'effects' and id(n) not in joined) else '')
            for n in ast.walk(tree) if isinstance(n, ast.Constant) and isinstance(n.value, str) and id(n) not in docs]


def grep_old():
    """Uses of the old effects folder in code: a string naming lib/effects/ or 'effects' joined onto lib."""
    hits = []
    files = [p for p in list((V / 'code').rglob('*.py')) + list((V / 'scripts').glob('*.py')) + list(V.glob('*.ipynb'))
             if '__pycache__' not in str(p) and p.name != 'check_provenance.py']
    for p in files:
        if p.suffix == '.ipynb':
            cells = [''.join(c['source']) for c in json.loads(p.read_text()).get('cells', []) if c['cell_type'] == 'code']
            strings = [(i, t) for i, c in enumerate(cells) for _, t in code_strings(c)]
        else:
            strings = code_strings(p.read_text(errors='replace'))
        rel = str(p.relative_to(V))
        for line, text in strings:
            if 'lib/effects/' in text or re.fullmatch(r'lib/effects', text.rstrip('/')) or text == 'effects':
                hits.append(dict(file=rel, line=line, text=text.strip()[:160], strict=any(rel.startswith(x) for x in STRICT)))
    return hits


def main():
    tables = sorted(p for p in E.glob('*.h5ad'))
    rows = [check_table(p) for p in tables]
    hits = grep_old()
    out = dict(tables=rows, old_path_hits=hits, ok=all(r['ok'] for r in rows) and not any(h['strict'] for h in hits))
    (V / 'results/qc').mkdir(parents=True, exist_ok=True)
    (V / 'results/qc/provenance_check.json').write_text(json.dumps(out, indent=1) + '\n')
    t = pd.DataFrame(rows)
    cols = [c for c in ['source', 'rows', 'kept_in_qc', 'se_floor', 'fdr', 'MAX_ONTARGET_LFC', 'MAX_ONTARGET_LFC_from',
                        'MIN_CELLS', 'MIN_CELLS_from', 'built_from', 'ok'] if c in t]
    pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 70)
    print(t[cols].to_string(index=False))
    for r in rows:
        if r['fails']: print('FAIL', r['source'], r['fails'])
    print(f'\nold path lib/effects/ mentions: {len(hits)} ({sum(h["strict"] for h in hits)} in current v4 code)')
    for h in hits: print(f"  {'STRICT ' if h['strict'] else ''}{h['file']}:{h['line']}  {h['text']}")
    print('PROVENANCE OK' if out['ok'] else 'PROVENANCE FAILED')
    sys.exit(0 if out['ok'] else 1)


if __name__ == '__main__':
    main()
