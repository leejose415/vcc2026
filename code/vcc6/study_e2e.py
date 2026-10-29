"""End-to-end inputs for the shrinkage study's four harness runs (Jose 2026-10-27, after analyze): baseline and
Stage C on HCT116final (the HCT116/293T fold: HCT116 and HEK293T both held out) and on H1 Eval A (the H1 fold).

  baseline  cache/base_<line>__<fold>.npz  the current model (C+G) with the fold's tables removed
  stageC    stageC_<line>.npz              c x B_k x p^A (shrink_analyze), same knockdowns
Both keep the baseline's DE calls and n_sources (the study changes magnitude only). Writes
<study>/e2e/inputs/<line>_<arm>.npz in run_line's --lfc-npz format (knockdowns, genes = the panel, lfc, de_sig,
n_sources); run_line predict --lfc-npz then generates cells, run_line score scores them on the eval.

  ~/venv/bin/python code/vcc6/study_e2e.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

V = Path.home() / 'vcc'
OUT = V / 'results/stage1/2026-10-27_shrinkage_study'
RUNS = {'HCT116': ('HCT116+HEK293T', V / 'results/target_lists/final_HCT116.txt'),
        'H1': ('H1', V / 'results/target_lists/evalA_h1_292.txt')}


def main():
    panel = pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str).tolist()
    d = OUT / 'e2e/inputs'; d.mkdir(parents=True, exist_ok=True)
    for line, (fold, lst) in RUNS.items():
        want = [k for k in lst.read_text().split() if k]
        b = np.load(OUT / f'cache/base_{line}__{fold}.npz', allow_pickle=True)
        c = np.load(OUT / f'stageC_{line}.npz', allow_pickle=True)
        for z in (b, c):
            assert list(map(str, z['genes'])) == panel, f'{line}: genes are not the panel order'
        bi = {str(k): i for i, k in enumerate(b['knockdowns'])}; ci = {str(k): i for i, k in enumerate(c['knockdowns'])}
        missing = [k for k in want if k not in bi or k not in ci]
        assert not missing, f'{line}: no prediction for {missing[:10]}'
        ib = [bi[k] for k in want]; ic = [ci[k] for k in want]
        sig, n_src = b['de_sig'][ib], b['n_sources'][ib]
        for arm, lfc in (('baseline', b['lfc'][ib]), ('stageC', c['lfc'][ic])):
            np.savez_compressed(d / f'{line}_{arm}.npz', knockdowns=np.array(want, dtype=object),
                                genes=np.array(panel, dtype=object), lfc=lfc.astype(np.float32), de_sig=sig, n_sources=n_src)
            print(f'{line} {arm}: {len(want)} knockdowns, |lfc| median {np.median(np.abs(lfc[lfc != 0])):.4f}, '
                  f'calls median {np.median(sig.sum(1)):.0f}')


if __name__ == '__main__':
    main()
