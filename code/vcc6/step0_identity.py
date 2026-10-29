"""Step 0a check: predictions.npz on 20 RPE1 knockdowns, before (the notebook's model cells, verbatim) and after
(code/vcc6/stage1.py with the salt and calibration_stats changes). Byte comparison of the files, and of every array.

  ~/venv/bin/python code/vcc6/step0_identity.py before <verbatim copy of the notebook cells as a .py>
  ~/venv/bin/python code/vcc6/step0_identity.py after
  ~/venv/bin/python code/vcc6/step0_identity.py compare
"""
import hashlib
import importlib.util
import json
import sys
import zipfile
from pathlib import Path

import numpy as np

V = Path.home() / 'vcc'
OUT = V / 'results/stage1/2026-10-27_rpe1/step0'
DIRS = [str(V / 'lib/effects_v5'), str(V / 'lib/effects_v5_extra')]


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def run(tag, before_py):
    """One run per process (each loads all 117 tables, ~15 GB)."""
    OUT.mkdir(parents=True, exist_ok=True)
    kd = OUT / 'rpe1_first20.txt'
    kd.write_text('\n'.join((V / 'results/target_lists/rpe1_dev_300.txt').read_text().split()[:20]) + '\n')
    m = load(before_py, 'stage1_before') if tag == 'before' else load(V / 'code/vcc6/stage1.py', 'stage1_after')
    print(f'== {tag}: {m.__file__}', flush=True)
    m.run_predict(DIRS, 'replogle_rpe1', str(OUT / tag), knockdowns_file=str(kd))


def compare():
    a, b = OUT / 'before/predictions.npz', OUT / 'after/predictions.npz'
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    rep = dict(file_bytes_identical=a.read_bytes() == b.read_bytes(), sha_before=sha(a), sha_after=sha(b))
    with zipfile.ZipFile(a) as za, zipfile.ZipFile(b) as zb:
        rep['member_bytes_identical'] = {n: za.read(n) == zb.read(n) for n in za.namelist()}
        rep['zip_timestamps'] = {n: [list(za.getinfo(n).date_time), list(zb.getinfo(n).date_time)] for n in za.namelist()}
    A, B = np.load(a, allow_pickle=True), np.load(b, allow_pickle=True)
    rep['arrays_equal'] = {k: bool(np.array_equal(A[k], B[k])) for k in A.files}
    sb = json.loads((OUT / 'before/summary.json').read_text()); sa = json.loads((OUT / 'after/summary.json').read_text())
    rep['summary_keys_added'] = sorted(set(sa) - set(sb))
    rep['summary_equal_except_fit'] = {k: sb[k] == sa[k] for k in sb}
    rep['fit'] = sa.get('fit')
    (OUT / 'identity.json').write_text(json.dumps(rep, indent=1) + '\n')
    print(json.dumps({k: rep[k] for k in ('file_bytes_identical', 'member_bytes_identical', 'arrays_equal', 'summary_keys_added')}, indent=1))


if __name__ == '__main__':
    if sys.argv[1] == 'compare':
        compare()
    else:
        run(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
