"""Raw single-cell files for every perturbation source, one per cell line, with one labelling scheme.

For each dataset: stream the raw counts and write lib/cells/<NAME>.h5ad (CSR, integer counts, gzip) with the same obs
columns everywhere. No cell is dropped by default: total_umi is stored so a UMI cut can be applied later
(--min-umi N applies one at build time).

  cell_line        the cell line / type (K562, H1, HFF, ...)
  dataset          the source study (Replogle2022_GWPS, Arc2025_H1, Hein2021, ...)
  condition        a study condition where there is one; else ''
  target_gene      the knocked-down gene's symbol; 'non-targeting' for every control; 'A+B' for a pair
  n_targets        0 for controls, 1 for a single knockdown, 2 for a pair, -1 for a cell with no assigned guide
                   (PerturBase 'None'; target_gene 'unassigned')
  perturbation     the source's own label, unchanged
  guide            the source's guide identifier where it has one
  batch            the source's batch / gem group where it has one
  total_umi        UMIs per cell (see umi_source)
  umi_source       'library' = the source's per-cell UMI total over ALL genes (Replogle / Nadig release only
                   ~8,250 expressed genes, so the matrix sum would undercount); 'matrix' = sum over the stored genes
  n_genes_detected genes with >= 1 count in the stored matrix

var: index = gene symbol made unique (a repeat gets -1, -2), gene_symbol = the source symbol, gene_id = Ensembl ID
where the source has one. Genes are the source's own, not
re-indexed to the challenge's 18,533.

Usage:  ~/venv/bin/python scripts/build_raw_cells.py --datasets K562 [--fetch] [--delete-raw] [--min-umi 0] [--limit 50000]
Writes: lib/cells/<NAME>.h5ad and a row per dataset in lib/cells/manifest.csv
"""
import argparse
import re
import json
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem, write_elem

V = Path.home() / 'vcc'
LIB = V / 'lib'
OUT = LIB / 'cells'
PB = LIB / 'downloads/perturbase'
MIN_UMI = 0                     # no cell is dropped unless --min-umi is given
REPORT_UMI = 1000               # the manifest counts cells below this, without dropping them
UNASSIGNED = 'unassigned'
CONTROL = 'non-targeting'
CHUNK = 20000

# name -> (reader, file, labels)
DATASETS = {
    'K562':    ('replogle', LIB / 'replogle/K562_gwps_raw_singlecell_01.h5ad', dict(cell_line='K562', dataset='Replogle2022_GWPS')),
    'RPE1':    ('replogle', LIB / 'downloads/rpe1/rpe1_raw_singlecell_01.h5ad', dict(cell_line='RPE1', dataset='Replogle2022_essential')),
    'HepG2':   ('replogle', LIB / 'downloads/nadig_hepg2/GSE264667_hepg2_raw_singlecell_01.h5ad', dict(cell_line='HepG2', dataset='Nadig2025')),
    'Jurkat':  ('replogle', LIB / 'downloads/nadig_jurkat/GSE264667_jurkat_raw_singlecell_01.h5ad', dict(cell_line='Jurkat', dataset='Nadig2025')),
    'H1':      ('h1', LIB / 'h1/H1_full_2025.h5ad', dict(cell_line='H1', dataset='Arc2025_H1')),
    'HCT116_targets':  ('orion', LIB / 'orion/HCT116_cells_targets.h5ad', dict(cell_line='HCT116', dataset='XAtlasOrion2025_subset')),
    'HEK293T_targets': ('orion', LIB / 'orion/HEK293T_cells_targets.h5ad', dict(cell_line='HEK293T', dataset='XAtlasOrion2025_subset')),
    # PerturBase (perturbase.cn) raw releases: one raw.h5ad per study, obs 'gene' only ('CTRL' = control,
    # 'A,B' = a pair, 'None' = no guide detected -> dropped), gene symbols only, full transcriptome.
    'hESC':            ('perturbase', {'': PB / 'PRJNA524331/raw.h5ad'}, dict(cell_line='hESC', dataset='Genga2019',
                        note='CRISPRi during definitive-endoderm differentiation; arrayed screen')),
    'HFF':             ('perturbase', {'': PB / 'CRISPRi_perturb_host/raw.h5ad'}, dict(cell_line='HFF', dataset='Hein2021',
                        note='human foreskin fibroblasts; host-factor knockdowns; the file holds the uninfected 00h cells only (knockdowns and controls)')),
}
GEO = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/suppl/'
PB_URL = 'http://www.perturbase.cn/static/{0}/{0}.raw.tar.gz'      # the site refuses https
SOURCES = {                                                         # where --fetch gets each missing raw file
    'K562': 'gs://vcc2026/replogle/K562_gwps_raw_singlecell_01.h5ad',
    'RPE1': 'https://ndownloader.figshare.com/files/35775606',
    'HepG2': GEO + 'GSE264667_hepg2_raw_singlecell_01.h5ad', 'Jurkat': GEO + 'GSE264667_jurkat_raw_singlecell_01.h5ad',
    'H1': 'gs://vcc2026/lib/h1/H1_full_2025.h5ad',
    'HCT116_targets': 'gs://vcc2026/orion/HCT116_cells_targets.h5ad', 'HEK293T_targets': 'gs://vcc2026/orion/HEK293T_cells_targets.h5ad',
}
PB_CONTROLS = {'CTRL', 'control', 'Control', 'NT', 'non-targeting'}
PB_UNASSIGNED = {'None', 'nan', '', 'multiple', 'Multiple'}


# ---------------------------------------------------------------- labels
def standard_target(raw, control_labels, sep=None):
    """Source label -> (target_gene, n_targets). Controls become 'non-targeting'; a pair 'A,B' becomes 'A+B'."""
    s = str(raw).strip()
    if s in control_labels or s.lower().startswith(('non-targeting', 'nontargeting', 'negctrl', 'neg_ctrl', 'safe-targeting')):
        return CONTROL, 0
    parts = [p.strip() for p in (s.split(sep) if sep else [s]) if p.strip()]
    parts = [p for p in parts if p not in control_labels and not p.lower().startswith(('negctrl', 'non-targeting'))]
    if not parts: return CONTROL, 0
    return '+'.join(parts), len(parts)


def labelled_obs(n, labels, perturbation, target, n_targets, guide=None, batch=None, condition=''):
    return pd.DataFrame(dict(
        cell_line=labels['cell_line'], dataset=labels['dataset'], condition=condition,
        target_gene=np.asarray(target, object), n_targets=np.asarray(n_targets, np.int8),
        perturbation=np.asarray(perturbation, object).astype(str),
        guide=np.asarray(guide if guide is not None else [''] * n, object).astype(str),
        batch=np.asarray(batch if batch is not None else [''] * n, object).astype(str)))


# ---------------------------------------------------------------- streaming writer
class CSRWriter:
    """Appends integer CSR rows to an h5ad X group (anndata encoding), gzip-compressed, then writes obs / var."""
    def __init__(self, path, n_vars, dtype=np.int32):
        self.path = Path(path); self.tmp = self.path.with_suffix('.partial.h5ad')
        self.f = h5py.File(self.tmp, 'w'); self.n_vars = n_vars; self.rows = 0; self.nnz = 0
        g = self.f.create_group('X'); g.attrs['encoding-type'] = 'csr_matrix'; g.attrs['encoding-version'] = '0.1.0'
        kw = dict(maxshape=(None,), chunks=(1 << 20,), compression='gzip', compression_opts=4, shuffle=True)
        self.data = g.create_dataset('data', (0,), dtype=dtype, **kw)
        self.indices = g.create_dataset('indices', (0,), dtype=np.int32, **kw)
        self.indptr = [np.zeros(1, np.int64)]
        self.dtype = np.dtype(dtype)

    def append(self, block):
        block = sp.csr_matrix(block)
        if block.nnz:
            if block.data.max() > np.iinfo(self.dtype).max or not np.all(block.data == np.round(block.data)):
                raise ValueError('counts are not integers within the stored dtype')
            n = block.nnz
            for ds, values in ((self.data, block.data.astype(self.dtype)), (self.indices, block.indices.astype(np.int32))):
                ds.resize((self.nnz + n,)); ds[self.nnz:] = values
        self.indptr.append(block.indptr[1:].astype(np.int64) + self.nnz)
        self.nnz += block.nnz; self.rows += block.shape[0]

    def close(self, obs, var, uns):
        g = self.f['X']; g.create_dataset('indptr', data=np.concatenate(self.indptr))
        g.attrs['shape'] = np.array([self.rows, self.n_vars])
        if len(obs) != self.rows: raise ValueError(f'obs has {len(obs)} rows, X has {self.rows}')
        obs = obs.copy(); obs.index = obs.index.astype(str)
        for k, v in dict(obs=obs, var=var, uns=uns).items(): write_elem(self.f, k, v)
        for k in ('obsm', 'varm', 'obsp', 'varp', 'layers'): write_elem(self.f, k, {})
        self.f.attrs['encoding-type'] = 'anndata'; self.f.attrs['encoding-version'] = '0.1.0'
        self.f.close(); self.tmp.replace(self.path)


# ---------------------------------------------------------------- readers
def read_replogle(path, labels, limit):
    """Replogle / Nadig raw release: dense float counts, obs gene / gene_id / sgID_AB / gem_group / UMI_count."""
    f = h5py.File(path, 'r'); obs = read_elem(f['obs']); var = read_elem(f['var'])
    n = len(obs) if not limit else min(limit, len(obs))
    obs = obs.iloc[:n]
    total = obs['UMI_count'].to_numpy(float); keep = total >= MIN_UMI
    target, k = zip(*(standard_target(x, {'non-targeting'}) for x in obs['gene'].astype(str)))
    out = labelled_obs(n, labels, obs['gene'], target, k, guide=obs.get('sgID_AB'), batch=obs.get('gem_group'))
    out.index = obs.index.astype(str); out['total_umi'] = total; out['umi_source'] = 'library'
    genes = pd.DataFrame({'gene_id': var.index.astype(str)}, index=pd.Index(var['gene_name'].astype(str), name='gene'))
    def blocks():
        for start in range(0, n, CHUNK):
            stop = min(start + CHUNK, n); m = keep[start:stop]
            yield f['X'][start:stop][m], m
    return out, keep, genes, blocks, f


def read_csr(path, labels, limit, kind):
    """CSR h5ad (H1, Orion subsets): matrix-sum UMIs, labels from obs."""
    f = h5py.File(path, 'r'); obs = read_elem(f['obs']); var = read_elem(f['var'])
    n = len(obs) if not limit else min(limit, len(obs))
    obs = obs.iloc[:n]; X = f['X']; indptr = X['indptr'][:n + 1]
    total = np.zeros(n)
    for start in range(0, n, CHUNK * 5):                       # one pass for the per-cell totals
        stop = min(start + CHUNK * 5, n); a, b = indptr[start], indptr[stop]
        d = X['data'][a:b]; rows = np.repeat(np.arange(stop - start), np.diff(indptr[start:stop + 1]))
        total[start:stop] = np.bincount(rows, weights=d, minlength=stop - start)
    keep = total >= MIN_UMI
    tcol = 'target_gene'
    target, k = zip(*(standard_target(x, {'non-targeting'}) for x in obs[tcol].astype(str)))
    out = labelled_obs(n, labels, obs[tcol], target, k, guide=obs.get('guide_id'), batch=obs.get('batch'))
    out.index = obs.index.astype(str); out['total_umi'] = total; out['umi_source'] = 'matrix'
    gid = var['gene_id'].astype(str).to_numpy() if 'gene_id' in var else np.full(len(var), '', object)
    genes = pd.DataFrame({'gene_id': gid}, index=pd.Index(var.index.astype(str), name='gene'))
    def blocks():
        for start in range(0, n, CHUNK):
            stop = min(start + CHUNK, n); a, b = indptr[start], indptr[stop]
            block = sp.csr_matrix((X['data'][a:b], X['indices'][a:b], indptr[start:stop + 1] - a), shape=(stop - start, len(var)))
            m = keep[start:stop]; yield block[m], m
    return out, keep, genes, blocks, f


class _Closed:
    def close(self): pass


def pb_target(raw):
    """PerturBase label -> (target_gene, n_targets), or None for a cell with no assigned guide. A trailing guide
    number ('FDPS_2') is dropped; 'CTRL' and negative-control guides are controls."""
    s = str(raw).strip()
    if s in PB_UNASSIGNED: return None
    parts = [re.sub(r'_\d+$', '', p.strip()) for p in s.split(',') if p.strip()]
    parts = [p for p in parts if p not in PB_CONTROLS and not p.lower().startswith(('negctrl', 'non-targeting', 'nt_', 'safe'))]
    return ('+'.join(sorted(parts)), len(parts)) if parts else (CONTROL, 0)


def read_perturbase(files, labels, limit):
    """One or more PerturBase raw.h5ad files (one per condition), read whole: they are small."""
    import anndata as ad
    obs_parts, mats, genes = [], [], None
    for condition, path in files.items():
        a = ad.read_h5ad(path)
        if genes is None: genes = a.var_names.astype(str)
        elif list(a.var_names.astype(str)) != list(genes): raise ValueError(f'{path}: genes differ between conditions')
        X = sp.csr_matrix(a.X); n = a.n_obs if not limit else min(limit, a.n_obs); X = X[:n]
        lab = [pb_target(x) for x in a.obs['gene'].astype(str)[:n]]
        assigned = np.array([x is not None for x in lab])
        target = [x[0] if x else UNASSIGNED for x in lab]; k = [x[1] if x else -1 for x in lab]
        o = labelled_obs(n, labels, a.obs['gene'].astype(str)[:n], target, k, condition=condition)
        o.index = (a.obs_names[:n].astype(str) + (f'-{condition}' if condition else '')); o['total_umi'] = np.asarray(X.sum(1)).ravel()
        o['umi_source'] = 'matrix'; o['assigned'] = assigned
        obs_parts.append(o); mats.append(X)
    obs = pd.concat(obs_parts); X = sp.vstack(mats).tocsr()
    obs.pop('assigned'); keep = obs['total_umi'].to_numpy() >= MIN_UMI
    obs['n_targets'] = obs['n_targets'].astype(np.int8)
    gene_table = pd.DataFrame({'gene_id': np.full(len(genes), '', object)}, index=pd.Index(genes, name='gene'))
    def blocks():
        for start in range(0, X.shape[0], CHUNK):
            stop = min(start + CHUNK, X.shape[0]); m = keep[start:stop]; yield X[start:stop][m], m
    return obs, keep, gene_table, blocks, _Closed()


def build(name, limit=0):
    reader, path, labels = DATASETS[name]
    paths = list(path.values()) if isinstance(path, dict) else [path]
    missing = [str(p) for p in paths if not Path(p).exists()]
    if missing: raise FileNotFoundError(f'{name}: {missing} missing; download first')
    t0 = time.time()
    if reader == 'replogle': obs, keep, genes, blocks, f = read_replogle(path, labels, limit)
    elif reader == 'perturbase': obs, keep, genes, blocks, f = read_perturbase(path, labels, limit)
    else: obs, keep, genes, blocks, f = read_csr(path, labels, limit, reader)
    OUT.mkdir(parents=True, exist_ok=True)
    dest = OUT / f'{name}{"_sample" if limit else ""}.h5ad'
    w = CSRWriter(dest, len(genes)); detected = []
    for i, (block, _) in enumerate(blocks()):
        block = sp.csr_matrix(block); w.append(block); detected.append(np.diff(block.indptr))
        if i % 10 == 0: print(f'  {name}: {w.rows:,} cells written ({time.time() - t0:.0f}s)', flush=True)
    kept = obs[keep].copy(); kept['n_genes_detected'] = np.concatenate(detected) if detected else np.zeros(0, int)
    # Unique gene names (a repeated symbol gets '-1', '-2', ...); the source's symbol and Ensembl ID stay as columns.
    symbols = genes.index.astype(str).to_numpy(); unique = pd.Series(symbols)
    dup = unique.duplicated(keep='first'); unique[dup] = unique[dup] + '-' + unique[dup].groupby(unique[dup]).cumcount().add(1).astype(str)
    genes = genes.assign(gene_symbol=symbols); genes.index = pd.Index(unique.to_numpy(), name='gene')
    for c in ('cell_line', 'dataset', 'condition', 'target_gene', 'perturbation', 'guide', 'batch', 'umi_source'):
        kept[c] = kept[c].astype('category')
    low = obs['total_umi'].to_numpy() < MIN_UMI; unassigned = kept['n_targets'].to_numpy() < 0
    below_report = int((kept['total_umi'].to_numpy() < REPORT_UMI).sum())
    uns = dict(source_file=json.dumps({k: str(v) for k, v in path.items()}) if isinstance(path, dict) else str(path), reader=reader,
               min_umi=MIN_UMI, cells_before=int(len(obs)), cells_after=int(keep.sum()), dropped_low_umi=int(low.sum()),
               unassigned_cells=int(unassigned.sum()), cells_below_1000_umi=below_report, control_label=CONTROL, note=labels.get('note', ''),
               built=time.strftime('%Y-%m-%d %H:%M'))
    w.close(kept, genes, uns); f.close()
    row = dict(name=name, cell_line=labels['cell_line'], dataset=labels['dataset'], cells_before=len(obs), cells_after=int(keep.sum()),
               min_umi=MIN_UMI, dropped_low_umi=int(low.sum()), cells_below_1000_umi=below_report, unassigned_cells=int(unassigned.sum()),
               knockdowns=int(kept.loc[kept.n_targets == 1, 'target_gene'].nunique()),
               pairs=int(kept.loc[kept.n_targets > 1, 'target_gene'].nunique()), control_cells=int((kept.n_targets == 0).sum()),
               genes=len(genes), median_umi=float(kept.total_umi.median()), umi_source=kept.umi_source.iloc[0] if len(kept) else '',
               file=str(dest), gb=round(dest.stat().st_size / 1e9, 2), minutes=round((time.time() - t0) / 60, 1))
    print(json.dumps(row), flush=True)
    if not limit:
        import fcntl
        manifest = OUT / 'manifest.csv'
        with open(OUT / '.manifest.lock', 'w') as lock:          # two builds may finish at once
            fcntl.flock(lock, fcntl.LOCK_EX)
            old = pd.read_csv(manifest) if manifest.exists() else pd.DataFrame()
            pd.concat([old[old.name != name] if len(old) else old, pd.DataFrame([row])]).to_csv(manifest, index=False)
    return row


def fetch(name):
    """Download any missing raw file of `name`; returns the files it fetched (so they can be deleted after the build)."""
    import subprocess, tarfile
    reader, path, _ = DATASETS[name]
    fetched = []
    for dest in (path.values() if isinstance(path, dict) else [path]):
        dest = Path(dest)
        if dest.exists(): continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        if reader == 'perturbase':
            study = dest.parent.name; archive = dest.parent / f'{study}.raw.tar.gz'
            subprocess.run(['curl', '-sSL', '--retry', '5', '-o', str(archive), PB_URL.format(study)], check=True)
            with tarfile.open(archive) as t:
                member = t.getmember('raw.h5ad')                      # the only file used; nothing else is extracted
                with t.extractfile(member) as src, open(dest, 'wb') as out:
                    while chunk := src.read(1 << 24): out.write(chunk)
            archive.unlink()
        else:
            url = SOURCES[name]; part = dest.with_suffix('.part')
            cmd = ['gcloud', 'storage', 'cp', url, str(part)] if url.startswith('gs://') else ['curl', '-sSL', '--retry', '5', '-o', str(part), url]
            subprocess.run(cmd, check=True); part.replace(dest)
        print(f'  fetched {dest} ({dest.stat().st_size / 1e9:.1f} GB)', flush=True); fetched.append(dest)
    return fetched


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--datasets', nargs='+', required=True, choices=list(DATASETS))
    p.add_argument('--limit', type=int, default=0, help='only the first N cells, written to <NAME>_sample.h5ad (a size test)')
    p.add_argument('--fetch', action='store_true', help='download missing raw files first (bucket, GEO, figshare, PerturBase)')
    p.add_argument('--delete-raw', action='store_true', help='after a successful build, delete the raw files this run fetched')
    p.add_argument('--min-umi', type=int, default=0, help='drop cells with fewer total UMIs (default 0: keep every cell)')
    a = p.parse_args(); MIN_UMI = a.min_umi
    for name in a.datasets:
        fetched = fetch(name) if a.fetch else []
        build(name, a.limit)
        if a.delete_raw and not a.limit:
            for f in fetched: f.unlink(); print(f'  deleted fetched raw {f}', flush=True)
