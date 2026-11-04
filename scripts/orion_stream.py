"""X-Atlas/Orion per-cell subset for a list of targets, streamed from Hugging Face (Xaira-Therapeutics/X-Atlas-Orion,
CC-BY-NC-SA-4.0), in the format of lib/orion/<LINE>_cells_targets.h5ad (the A-C subset of 2026-10-08):
  obs   target_gene ('non-targeting' for the controls), context (= the line); index '<parquet file>:<cell barcode>'
  var   the 18,533 challenge genes (data/gene_names.csv order), measured = the gene is in Orion's gene metadata
  X     raw integer counts, CSR (Orion token -> gene symbol -> challenge column, summed)
Cells: every cell of every listed target (the release is already guide-filtered: pass_guide_filter = 1, two guides),
plus N_CTRL controls ('Non-Targeting'): the N_CTRL with the smallest sha256(salt + barcode), so the draw does not depend
on the order the batches finish in. Each parquet batch is downloaded, read and deleted (about 47 / 80 GB in total).

  ~/venv/bin/python scripts/orion_stream.py --line HCT116 --targets LIST [LIST ...] --out lib/orion/HCT116_cells_targets.h5ad
"""
import argparse
import hashlib
import heapq
import json
import subprocess
import sys
import tempfile
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

V = Path.home() / 'vcc'
REPO = 'https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion'
CTRL_SRC, CONTROL = 'Non-Targeting', 'non-targeting'
N_CTRL, SALT = 20000, 'orion-controls-2026'


def log(msg):
    print(f'[{time.strftime("%H:%M:%S")}] {msg}', flush=True)


def batch_files(line):
    with urllib.request.urlopen(f'https://huggingface.co/api/datasets/Xaira-Therapeutics/X-Atlas-Orion/tree/main/data') as r:
        tree = json.load(r)
    return sorted(f['path'] for f in tree if f.get('type') == 'file' and f['path'].startswith(f'data/{line}_Batch'))


def fetch(path, dest):
    subprocess.run(['curl', '-sSfL', '--retry', '5', '-o', str(dest), f'{REPO}/resolve/main/{path}'], check=True)


def ctrl_key(barcode):
    return int.from_bytes(hashlib.sha256((SALT + barcode).encode()).digest()[:8], 'big')


def read_batch(path, tmp, targets, col_of_token, n_genes):
    """One parquet batch -> (target rows as (name, label, data, cols), control rows as (key, name, data, cols))."""
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    local = Path(tmp) / Path(path).name
    fetch(path, local)
    fname = Path(path).name
    want = pa.array(sorted(targets | {CTRL_SRC}))
    tgt, ctl, n = [], [], 0
    # small record batches, filtered before anything is materialised: a whole batch file decompressed is several GB
    for rb in pq.ParquetFile(local).iter_batches(batch_size=1024, columns=['cell_barcode', 'gene_target', 'gene_token_id',
                                                                          'gene_expression']):
        n += rb.num_rows
        rb = rb.filter(pc.is_in(rb.column('gene_target'), value_set=want))
        if rb.num_rows == 0:
            continue
        gt = rb.column('gene_target').to_pylist(); bc = rb.column('cell_barcode').to_pylist()
        tok = rb.column('gene_token_id'); ex = rb.column('gene_expression')
        for i in range(rb.num_rows):
            cols = col_of_token[np.asarray(tok[i].as_py(), dtype=np.int64)]
            vals = np.asarray(ex[i].as_py(), dtype=np.float64)
            ok = cols >= 0
            c, v = cols[ok], vals[ok]
            if c.size:                                               # sum tokens that map to one challenge gene
                order = np.argsort(c, kind='stable'); c, v = c[order], v[order]
                u, start = np.unique(c, return_index=True); v = np.add.reduceat(v, start); c = u
            name = f'{fname}:{bc[i]}'
            if gt[i] == CTRL_SRC:
                ctl.append((ctrl_key(bc[i]), name, v.astype(np.int32), c.astype(np.int32)))
            else:
                tgt.append((name, gt[i], v.astype(np.int32), c.astype(np.int32)))
    local.unlink()
    return fname, n, tgt, ctl


def main(line, lists, out, workers, n_batches=None):
    import pyarrow.parquet as pq
    import anndata as ad
    genes = pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str).tolist()
    targets = set()
    for f in lists:
        if f.endswith('.csv'):
            targets |= set(pd.read_csv(f)['target_gene'].astype(str))
        else:
            targets |= {t for t in Path(f).read_text().split() if t}
    targets.discard(CONTROL)
    with tempfile.TemporaryDirectory(dir=str(Path(out).parent)) as tmp:
        meta = Path(tmp) / 'gene_metadata.parquet'; fetch('metadata/gene_metadata.parquet', meta)
        m = pq.read_table(meta).to_pandas()
        col = {g: j for j, g in enumerate(genes)}
        col_of_token = np.full(int(m.gene_token_id.max()) + 1, -1, dtype=np.int64)
        for tok, sym in zip(m.gene_token_id, m.gene_name):
            col_of_token[int(tok)] = col.get(str(sym), -1)
        measured = np.zeros(len(genes), bool); measured[col_of_token[col_of_token >= 0]] = True
        files = batch_files(line)
        if n_batches:
            files = files[:n_batches]                               # a test on the first batches only
        log(f'{line}: {len(files)} batches, {len(targets)} targets, {int(measured.sum())} of {len(genes)} challenge genes measured')
        rows, heap, seen, n_cells = [], [], {}, 0                    # heap: the N_CTRL smallest control keys (as -key)
        with ThreadPoolExecutor(workers) as pool:
            futs = [pool.submit(read_batch, f, tmp, targets, col_of_token, len(genes)) for f in files]
            for k, fu in enumerate(as_completed(futs), 1):
                fname, n, tgt, ctl = fu.result()
                n_cells += n; rows.extend(tgt)
                for key, name, v, c in ctl:
                    if len(heap) < N_CTRL:
                        heapq.heappush(heap, (-key, name, v, c))
                    elif -key > heap[0][0]:
                        heapq.heapreplace(heap, (-key, name, v, c))
                if k % 10 == 0 or k == len(files):
                    log(f'  {k}/{len(files)} batches, {n_cells:,} cells read, {len(rows):,} target cells kept')
    ctrl = sorted(heap, key=lambda x: x[1])
    rows.sort(key=lambda x: x[0])
    allr = [(n, lab, v, c) for n, lab, v, c in rows] + [(n, CONTROL, v, c) for _, n, v, c in ctrl]
    indptr = np.zeros(len(allr) + 1, np.int64)
    for i, r in enumerate(allr): indptr[i + 1] = indptr[i] + r[3].size
    data = np.concatenate([r[2] for r in allr]).astype(np.int32); idx = np.concatenate([r[3] for r in allr])
    X = sp.csr_matrix((data, idx, indptr), shape=(len(allr), len(genes)))
    obs = pd.DataFrame(dict(target_gene=pd.Categorical([r[1] for r in allr]), context=pd.Categorical([line] * len(allr))),
                       index=pd.Index([r[0] for r in allr], name='cell'))
    var = pd.DataFrame(dict(measured=measured), index=pd.Index(genes, name='gene_name'))
    a = ad.AnnData(X=X, obs=obs, var=var)
    a.uns['orion_stream'] = dict(line=line, targets_requested=len(targets), targets_found=int(obs.target_gene.nunique() - 1),
                                 missing=sorted(targets - set(obs.target_gene.astype(str))), controls=len(ctrl),
                                 control_rule=f'the {N_CTRL} smallest sha256("{SALT}" + barcode)', batches=len(files),
                                 cells_read=n_cells, source=REPO, license='CC-BY-NC-SA-4.0', lists=[str(f) for f in lists])
    a.write_h5ad(out, compression='gzip')
    log(f'wrote {out}: {a.n_obs:,} cells ({len(rows):,} target, {len(ctrl):,} control), {a.uns["orion_stream"]["targets_found"]} targets; '
        f'missing {a.uns["orion_stream"]["missing"][:10]}')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--line', required=True, choices=['HCT116', 'HEK293T'])
    p.add_argument('--targets', nargs='+', required=True, help='target lists: .txt (one per line) or a pert_counts .csv')
    p.add_argument('--out', required=True)
    p.add_argument('--workers', type=int, default=4)
    p.add_argument('--batches', type=int, default=None, help='first N batches only (a test)')
    a = p.parse_args()
    main(a.line, a.targets, a.out, a.workers, a.batches)
