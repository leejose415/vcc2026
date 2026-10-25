"""Merge the public 2025 H1 split files into one raw-count h5ad with ONE copy of the controls.

Each 2025 split (Training, Validation, Test) carries the same non-targeting control cells.
Naive concatenation triple-counts them, and copies of one cell then land in both halves of
the H1 benchmark's generator/scoring control split. Controls come from the first file only;
the other files' controls are compared with it. Rows are streamed in blocks (low memory).

Example: python vcc_merge_h1.py --inputs adata_Training.h5ad adata_Validation.h5ad adata_Test.h5ad \
             --out H1_full_2025.h5ad
"""
import argparse
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

from vcc_common import CONTROL, int64_indptr, save_json, validate_counts

try:
    from anndata.io import sparse_dataset, write_elem
except ImportError:
    from anndata.experimental import sparse_dataset, write_elem


def merge(inputs, output, obs_col="target_gene", control=CONTROL, block=4096):
    output = Path(output)
    if output.exists(): raise FileExistsError(f"{output} exists; delete it or choose another path")
    plans, var, control_names, seen, report = [], None, None, set(), dict(inputs=[])
    for path in inputs:
        h = ad.read_h5ad(path, backed="r")
        try:
            if var is None: var = h.var.copy()
            elif not h.var_names.equals(var.index): raise ValueError(f"{path}: gene axis differs from {inputs[0]}")
            is_control = h.obs[obs_col].astype(str).values == control
            names = set(h.obs_names[is_control].astype(str))
            if control_names is None:
                if not names: raise ValueError(f"{path}: no {control!r} cells in obs[{obs_col!r}]")
                control_names, keep, same = names, np.ones(h.n_obs, bool), True
            else:
                keep, same = ~is_control, names == control_names
                if not same: print(f"WARNING {path}: its controls differ from {inputs[0]}'s; only the first file's are kept")
            rows = np.flatnonzero(keep)
            obs = h.obs.iloc[rows].copy(); obs.index = obs.index.astype(str)
            duplicated = seen & set(obs.index)
            if duplicated: raise ValueError(f"{path}: {len(duplicated)} cell names already merged from an earlier file")
            seen |= set(obs.index); plans.append((path, rows, obs))
            report['inputs'].append(dict(path=str(path), cells_kept=int(len(rows)), controls_in_file=int(is_control.sum()),
                                         controls_match_first=bool(same)))
        finally:
            h.file.close()
    columns = [col for col in plans[0][2].columns if all(col in p[2].columns for p in plans)]
    obs = pd.concat([p[2][columns] for p in plans])
    for col in obs.columns:
        if obs[col].dtype == object or isinstance(obs[col].dtype, pd.CategoricalDtype):
            obs[col] = obs[col].astype(str).astype("category")
    partial = output.with_name(output.name + ".partial")
    output.parent.mkdir(parents=True, exist_ok=True); ds = None
    with h5py.File(partial, "w") as f:
        f.attrs["encoding-type"] = "anndata"; f.attrs["encoding-version"] = "0.1.0"
        write_elem(f, "obs", obs); write_elem(f, "var", var)
        for path, rows, _ in plans:
            h = ad.read_h5ad(path, backed="r")
            try:
                for start in range(0, len(rows), block):
                    m = h.X[rows[start:start + block]]
                    m = (m.tocsr() if sp.issparse(m) else sp.csr_matrix(np.asarray(m))).astype(np.float32)
                    validate_counts(m, str(path))
                    if ds is None: write_elem(f, "X", int64_indptr(m)); ds = sparse_dataset(f["X"])
                    else: ds.append(m)
            finally:
                h.file.close()
            print(f"  merged {path}: {len(rows):,} cells", flush=True)
    partial.rename(output)
    labels = obs[obs_col].astype(str)
    report.update(output=str(output), n_cells=int(len(obs)), n_controls=int((labels == control).sum()),
                  n_targets=int(labels[labels != control].nunique()), n_genes=int(len(var)))
    save_json(output.with_name(output.stem + "_merge_report.json"), report)
    print(f"wrote {output}: {report['n_cells']:,} cells, {report['n_controls']:,} controls, {report['n_targets']} targets")
    return report


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--inputs", nargs="+", required=True, help="split files; controls are taken from the FIRST")
    p.add_argument("--out", required=True); p.add_argument("--obs-col", default="target_gene")
    p.add_argument("--control", default=CONTROL); p.add_argument("--block", type=int, default=4096)
    a = p.parse_args()
    merge(a.inputs, a.out, a.obs_col, a.control, a.block)
