"""Config 1 is frozen (code/vcc6/CONFIG1_FROZEN.md): its values, the v7 notebook's data-loader and model cells, the
stage-1 model, run_line and the table set must match the recorded hashes."""
import hashlib
import json
import re
import sys
from pathlib import Path

V = Path.home() / 'vcc'
RECORD = (V / 'code/vcc6/CONFIG1_FROZEN.md').read_text()
WANT = dict(s_c=0.5, m=0.5, kappa=1.0, c1=1.0)


def recorded(name):
    return re.search(rf'sha256 {name} ([0-9a-f]{{64}})', RECORD).group(1)


def v7_code_cells():
    return [''.join(c['source']) for c in json.load(open(V / 'vcc_2026_v7.ipynb'))['cells'] if c['cell_type'] == 'code']


def test_config1_values():
    got = re.search(r'^CONFIG = dict\(([^)]*)\)', v7_code_cells()[0], re.M).group(1)
    assert {k.strip(): float(v) for k, v in (kv.split('=') for kv in got.split(','))} == WANT
    sys.path.insert(0, str(V / 'code/vcc6'))
    import scorer_step16 as S16
    assert S16.CONFIG1['arm'] == 'i' and all(p == WANT for p in S16.CONFIG1['params'].values())


def test_config1_is_frozen():
    cells = v7_code_cells()
    assert hashlib.sha256((cells[1] + '\n' + cells[2]).encode()).hexdigest() == recorded('cells')
    for name in ('stage1', 'run_line'):
        assert hashlib.sha256((V / f'code/vcc6/{name}.py').read_bytes()).hexdigest() == recorded(name)


def test_config1_tables():
    sys.path.insert(0, str(V / 'code/vcc6'))
    import stage1 as S
    stems = sorted(p.stem for p in S.discover_tables([str(V / 'lib/effects_v5'), str(V / 'lib/effects_v5_extra')]))
    assert len(stems) == 117
    assert hashlib.sha256('\n'.join(stems).encode()).hexdigest() == recorded('tables')
