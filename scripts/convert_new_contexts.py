"""Convert the new CRISPRi context datasets into raw single-cell files shaped like the existing QC inputs.

Writes lib/cells/<name>.h5ad: X = integer UMI counts (CSR int32, gzip), every gene the source released; var index =
gene symbol made unique (a repeat gets -1, -2), gene_symbol, gene_id = Ensembl ID; obs = target_gene (the knocked-down
gene's symbol, 'non-targeting' for every control), guide_id (the construct), batch (the lane the UMI cutoffs use), plus
source columns kept under author_* / descriptive names. Only single-target CRISPRi cells are written; the source files
are never modified. These files then go through the original QC sequence (qc_replogle_umi -> qc_replogle_apply ->
write_replogle_qc -> qc_perturbations) via the 'h1' layout (UMI = the matrix row sum, target IDs from var['gene_id']).

Scope (user decisions, 2026-10-26): CRISPRi only; single-target cells only; Datlinger 2021, Gasperini, Jiang, Dixit,
Datlinger 2017, Frangieh and Norman excluded.
  kolf      KOLF2.1J hiPSC (author-selected strong perturbations): counts from layers['counts'] (X is normalized);
            target gene_target ('NTC' -> control), guide gRNA, lane channel. Every cell carries one guide.
  adamson   Adamson 2016 K562, one file per experiment (10X001, 10X005, 10X010): construct = the perturbation label;
            target = the label before '_p<plasmid>'; controls 62(mod), 63(mod), Gal4-4(mod) and 3x_neg_ctrl (non-targeting
            / Gal4-targeting); 10X005 keeps only ATF6_only, IRE1_only (ERN1) and PERK_only (EIF2AK3) of its combinations;
            unassigned '*' dropped. No lane is recorded: one lane ('all').
  feng      Feng HipSci iPSC CRISPRi (gene x cell CSV.gz, streamed by scripts/feng_csv_to_csc.cpp). Assigned cells only
            (Guide_Call != 'unassigned'); NonTarget_* -> control; target = Guide_Call before its last '_'.
            Targeted screen: one file per donor line (Day3), lane = the 10x inlet (Batch).
            Genome-wide fitness / non-fitness: one file per screen, all donor lines and days pooled, lane = donor line;
            donor_line and day kept in obs (they have ~0-2 non-targeting cells per line: QC uses the rest-of-screen
            reference).
RNA columns map to metadata after stripping a leading 'MP-' / 'PC-' (and, in the targeted screen, '-D3_' -> '_').

Usage: ~/venv/bin/python scripts/convert_new_contexts.py kolf | adamson | feng_targeted | feng_gw [--force]
"""
import argparse
import gzip
import json
import subprocess
import sys
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem, write_elem

V = Path.home() / 'vcc'
OUT = V / 'lib/cells'
CONTROL = 'non-targeting'
FENG = V / 'lib/ipsc/feng_counts'
STREAMER_SRC = V / 'scripts/feng_csv_to_csc.cpp'
RESERVED = {'source', 'background', 'condition', 'batch', 'target_gene', 'is_control', 'target_gene_id', 'guide_id',
            'n_umi', 'qc_low_ref'}        # set by write_replogle_qc: only target_gene, guide_id and batch are written here


def log(msg):
    print(f'[{time.strftime("%H:%M:%S")}] {msg}', flush=True)


def unique_symbols(sym):
    seen, out = {}, []
    for s in map(str, sym):
        k = seen.get(s, 0); out.append(s if k == 0 else f'{s}-{k}'); seen[s] = k + 1
    return out


def write_raw(name, X, obs, symbols, gene_ids, meta):
    """Integer counts CSR int32 + standard obs / var, gzip; refuses fractional or negative counts."""
    X = X.tocsr(); X.sort_indices()
    d = X.data
    if d.size and (d.min() < 0 or not np.all(d == np.rint(d))): raise ValueError(f'{name}: counts are not non-negative integers')
    if d.size and d.max() > np.iinfo(np.int32).max: raise ValueError(f'{name}: count overflow')
    X = sp.csr_matrix((d.astype(np.int32), X.indices.astype(np.int32), X.indptr.astype(np.int64)), shape=X.shape)
    X.eliminate_zeros()
    bad = RESERVED & (set(obs.columns) - {'target_gene', 'guide_id', 'batch'})
    if bad: raise ValueError(f'{name}: obs columns reserved for write_replogle_qc: {bad}')
    obs = obs.copy()
    for c in obs.columns:
        if pd.api.types.is_string_dtype(obs[c]) or obs[c].dtype == object: obs[c] = obs[c].astype(str).astype('category')
    var = pd.DataFrame(dict(gene_symbol=list(map(str, symbols)), gene_id=list(map(str, gene_ids))),
                       index=pd.Index(unique_symbols(symbols), name='gene'))
    a = ad.AnnData(X=X, obs=obs, var=var)
    lab = obs['target_gene'].astype(str)
    a.uns['conversion'] = dict(meta, converted=time.strftime('%Y-%m-%d %H:%M'), script='scripts/convert_new_contexts.py',
                               cells=int(a.n_obs), genes=int(a.n_vars), controls=int((lab == CONTROL).sum()),
                               targets=int(lab[lab != CONTROL].nunique()), nnz=int(X.nnz))
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f'{name}.h5ad'; tmp = path.with_suffix('.h5ad.part')
    a.write_h5ad(tmp, compression='gzip'); tmp.rename(path)
    row = np.asarray(X.sum(1)).ravel()
    log(f'wrote {path}: {a.n_obs:,} cells x {a.n_vars:,} genes, {a.uns["conversion"]["controls"]:,} controls, '
        f'{a.uns["conversion"]["targets"]:,} targets, median {np.median(row):.0f} UMI')
    return path


def read_sparse(h, key, n_obs, n_vars):
    """An h5ad sparse group (CSC or CSR) fully into memory as a scipy matrix."""
    g = h[key]; enc = g.attrs['encoding-type']
    M = (sp.csc_matrix if enc == 'csc_matrix' else sp.csr_matrix)((g['data'][:], g['indices'][:], g['indptr'][:]), shape=(n_obs, n_vars))
    return M


# ---------------------------------------------------------------- KOLF
def kolf(force=False):
    name = 'kolf_kolf21j'
    if (OUT / f'{name}.h5ad').exists() and not force: return log(f'{name}: exists, skipped')
    path = V / 'lib/ipsc/kolf/KOLF_Strong_Perturbations.h5ad'
    with h5py.File(path, 'r') as h:
        obs = read_elem(h['obs']); var = read_elem(h['var'])
        log(f'KOLF: reading layers/counts ({h["layers/counts/data"].shape[0]:,} stored values)')
        X = read_sparse(h, 'layers/counts', len(obs), len(var))
    if not (obs['n_gRNA'].to_numpy() == 1).all(): raise ValueError('KOLF: cells with more than one guide')
    target = obs['gene_target'].astype(str).where(obs['gene_target'].astype(str) != 'NTC', CONTROL)
    o = pd.DataFrame(dict(target_gene=target.to_numpy(), guide_id=obs['gRNA'].astype(str).to_numpy(),
                          batch=obs['channel'].astype(str).to_numpy(), source_cell_id=obs.index.astype(str),
                          author_perturbation=obs['perturbation'].astype(str).to_numpy(),
                          author_target_knockdown=obs['target_knockdown'].to_numpy(),
                          author_perturbation_significant=obs['perturbation_significant'].to_numpy(),
                          author_batch=obs['batch'].astype(str).to_numpy()), index=obs.index.astype(str))
    write_raw(name, X, o, var.index.astype(str), var['gene_ids'].astype(str),
              dict(source=str(path), dataset='KOLF strong perturbations (KRAB-ZIM3-dCas9 KOLF2.1J hiPSC)', modality='CRISPRi',
                   background='KOLF2.1J', counts="layers['counts'] (X is normalized)", lane='channel', controls='NTC',
                   note='author-selected strong perturbations: not an unbiased sample of guide performance'))


# ---------------------------------------------------------------- Adamson
ADAMSON = {'10X001': 'GSM2406675_10X001', '10X005': 'GSM2406677_10X005', '10X010': 'GSM2406681_10X010'}
ADAMSON_CONTROL = ('62(mod)', '63(mod)', 'Gal4-4(mod)', '3x_neg_ctrl')
ADAMSON_ALIAS = {'ATF6_only': 'ATF6', 'IRE1_only': 'ERN1', 'PERK_only': 'EIF2AK3'}


def adamson_target(label):
    """Construct label -> target symbol, CONTROL, or None (dropped: combination, unassigned)."""
    if label in ('*', 'nan', ''): return None
    stem = label.rsplit('_p', 1)[0]
    if stem.startswith(ADAMSON_CONTROL): return CONTROL
    if stem in ADAMSON_ALIAS: return ADAMSON_ALIAS[stem]
    if '_' in stem: return None                                  # ATF6_PERK_IRE1, PERK_IRE1, ...: combinations
    return stem


def adamson(force=False):
    for exp, file in ADAMSON.items():
        name = f'adamson_k562_{exp.lower()}'
        if (OUT / f'{name}.h5ad').exists() and not force: log(f'{name}: exists, skipped'); continue
        path = V / f'lib/scperturb/AdamsonWeissman2016_{file}.h5ad'
        with h5py.File(path, 'r') as h:
            obs = read_elem(h['obs']); var = read_elem(h['var'])
            X = read_sparse(h, 'X', len(obs), len(var))
        label = np.array([x if isinstance(x, str) else '' for x in obs['perturbation']], dtype=object)   # missing -> '' (dropped)
        target = np.array([adamson_target(x) for x in label], dtype=object)
        keep = target != None                                                       # noqa: E711
        dropped = pd.Series(label[~keep]).value_counts().to_dict()
        o = pd.DataFrame(dict(target_gene=target[keep].astype(str), guide_id=label[keep], batch='all',
                              source_cell_id=obs.index.astype(str)[keep]), index=obs.index.astype(str)[keep])
        write_raw(name, X.tocsr()[np.flatnonzero(keep)], o, var.index.astype(str), var['ensembl_id'].astype(str),
                  dict(source=str(path), dataset=f'Adamson 2016 K562 CRISPRi ({exp})', modality='CRISPRi', background='K562',
                       controls=', '.join(ADAMSON_CONTROL), lane="none recorded: one lane ('all')",
                       dropped_cells=json.dumps(dropped)))


# ---------------------------------------------------------------- Feng
def feng_metadata(screen):
    m = pd.read_csv(FENG / f'{screen}_Cell-Metadata.tsv.gz', sep='\t', dtype=str)
    norm = m.Cell_ID.str.replace(r'^MP-', '', regex=True).str.replace(r'^PC-', '', regex=True)
    if screen == 'TargetedScreen': norm = norm.str.replace(r'-D3_', '_', regex=True)
    if norm.duplicated().any(): raise ValueError(f'{screen}: duplicated normalized cell IDs')
    m.index = norm
    with gzip.open(FENG / f'{screen}_RNA-UMI-Counts.csv.gz', 'rt') as f: columns = f.readline().rstrip('\n').split(',')[1:]
    if len(columns) != len(set(columns)): raise ValueError(f'{screen}: duplicated RNA columns')
    m = m.reindex(columns)                                     # RNA column order; a column without metadata -> NaN row
    log(f'{screen}: {len(columns):,} RNA columns, {int(m.Cell_ID.isna().sum())} without metadata')
    guide = m.Guide_Call.fillna('unassigned').astype(str)
    assigned = (guide != 'unassigned') & m.Cell_ID.notna() & m.Cell_Line.notna() & (m.Cell_Line.astype(str) != 'unassigned_cell_line')
    if guide.str.contains('[;,]').any(): raise ValueError(f'{screen}: multi-guide calls')
    target = guide.str.rsplit('_', n=1).str[0]
    target = target.where(~target.eq('NonTarget'), CONTROL)
    m['target_gene'] = target; m['guide_id'] = guide; m['assigned'] = assigned
    m['day'] = m.Batch.str.extract(r'(Day\d+)', expand=False)
    return m, columns


def run_streamer(screen, m, groups, work):
    """groups: list of (source index, boolean mask over RNA columns, row order) -> CSC blocks per source."""
    exe = work / 'feng_csv_to_csc'
    subprocess.run(['g++', '-O3', '-march=native', '-o', str(exe), str(STREAMER_SRC), '-lz'], check=True)
    n = len(m); src = np.full(n, -1, np.int32); row = np.zeros(n, np.int32)
    for k, mask in enumerate(groups):
        idx = np.flatnonzero(mask); src[idx] = k; row[idx] = np.arange(len(idx), dtype=np.int32)
    (work / 'assign.bin').write_bytes(src.tobytes() + row.tobytes())
    log(f'{screen}: streaming {FENG / f"{screen}_RNA-UMI-Counts.csv.gz"} into {len(groups)} source(s)')
    t0 = time.time()
    subprocess.run([str(exe), str(FENG / f'{screen}_RNA-UMI-Counts.csv.gz'), str(work / 'assign.bin'), str(len(groups)), str(work)], check=True)
    log(f'{screen}: streamed in {(time.time() - t0) / 60:.1f} min')
    feats = (work / 'genes.txt').read_text().splitlines()
    parts = [f.split(':') for f in feats]
    gene_ids = [p[0] for p in parts]; symbols = [p[1] if len(p) > 1 else p[0] for p in parts]
    return gene_ids, symbols


def load_block(work, k, n_cells, n_genes):
    data = np.fromfile(work / f'src{k}.data', dtype=np.uint32); rows = np.fromfile(work / f'src{k}.rows', dtype=np.uint32)
    ptr = np.fromfile(work / f'src{k}.colptr', dtype=np.uint64)
    if len(ptr) != n_genes + 1 or int(ptr[-1]) != len(data): raise ValueError(f'source {k}: inconsistent CSC block')
    return sp.csc_matrix((data.astype(np.int32), rows.astype(np.int32), ptr.astype(np.int64)), shape=(n_cells, n_genes))


def feng_targeted(force=False):
    screen = 'TargetedScreen'; m, cols = feng_metadata(screen)
    lines = sorted(m.loc[m.assigned, 'Cell_Line'].astype(str).unique())
    names = [f'feng_targeted_{l}' for l in lines]
    todo = [i for i, nm in enumerate(names) if force or not (OUT / f'{nm}.h5ad').exists()]
    if not todo: return log('feng targeted: all lines exist, skipped')
    work = OUT / '_feng_work_targeted'; work.mkdir(parents=True, exist_ok=True)
    masks = [(m.assigned & (m.Cell_Line.astype(str) == lines[i])).to_numpy() for i in todo]
    gene_ids, symbols = run_streamer(screen, m, masks, work)
    for k, i in enumerate(todo):
        sel = m[masks[k]]
        X = load_block(work, k, len(sel), len(gene_ids))
        o = pd.DataFrame(dict(target_gene=sel.target_gene.to_numpy(), guide_id=sel.guide_id.to_numpy(), batch=sel.Batch.astype(str).to_numpy(),
                              donor_line=sel.Cell_Line.astype(str).to_numpy(), day=sel.day.astype(str).to_numpy(),
                              source_cell_id=sel.Cell_ID.astype(str).to_numpy()), index=sel.index.astype(str))
        write_raw(names[i], X, o, symbols, gene_ids,
                  dict(source=str(FENG / f'{screen}_RNA-UMI-Counts.csv.gz'), dataset=f'Feng HipSci iPSC CRISPRi targeted screen, donor line {lines[i]}',
                       modality='CRISPRi', background=lines[i], lane='10x inlet (Batch)', controls='NonTarget_*',
                       note='assigned cells only (unassigned are not controls)'))
        for f in ('data', 'rows', 'colptr'): (work / f'src{k}.{f}').unlink()


def feng_gw(force=False):
    for screen, name in (('GenomeWideScreen_FitnessGenes', 'feng_gw_fitness'), ('GenomeWideScreen_NonFitnessGenes', 'feng_gw_nonfitness')):
        if (OUT / f'{name}.h5ad').exists() and not force: log(f'{name}: exists, skipped'); continue
        m, cols = feng_metadata(screen)
        work = OUT / f'_feng_work_{name}'; work.mkdir(parents=True, exist_ok=True)
        mask = m.assigned.to_numpy()
        gene_ids, symbols = run_streamer(screen, m, [mask], work)
        sel = m[mask]; X = load_block(work, 0, len(sel), len(gene_ids))
        o = pd.DataFrame(dict(target_gene=sel.target_gene.to_numpy(), guide_id=sel.guide_id.to_numpy(),
                              batch=sel.Cell_Line.astype(str).to_numpy(), donor_line=sel.Cell_Line.astype(str).to_numpy(),
                              day=sel.day.astype(str).to_numpy(), inlet=sel.Batch.astype(str).to_numpy(),
                              source_cell_id=sel.Cell_ID.astype(str).to_numpy()), index=sel.index.astype(str))
        write_raw(name, X, o, symbols, gene_ids,
                  dict(source=str(FENG / f'{screen}_RNA-UMI-Counts.csv.gz'), dataset=f'Feng HipSci iPSC CRISPRi {screen}, donor lines and days pooled',
                       modality='CRISPRi', background='HipSci iPSC (pooled donor lines)', lane='donor line',
                       controls='NonTarget_* (only a handful: QC uses the rest-of-screen reference)',
                       note='assigned cells only; donor_line and day kept in obs'))
        for f in ('data', 'rows', 'colptr'): (work / f'src0.{f}').unlink()


# ---------------------------------------------------------------- Tahoe-100M (PIE-relevant extra source; 2026-10-26)
TAHOE = V / 'lib/tahoe'
TAHOE_WORK = OUT / '_tahoe_work'
TAHOE_URL = 'https://huggingface.co/datasets/tahoebio/Tahoe-100M/resolve/main/data/train-{:05d}-of-03388.parquet'
TAHOE_SHARDS, TAHOE_DMSO_CAP, TAHOE_COND_CAP, TAHOE_SEED = 3388, 18400, 450, 2026


def tahoe_name(cell_name):
    return 'tahoe_' + ''.join(c.lower() if c.isalnum() else '_' for c in str(cell_name)).strip('_')


def tahoe_single_target_drugs():
    """Drugs with exactly one annotated protein target -> that target's symbol."""
    d = pd.read_parquet(TAHOE / 'metadata/drug_metadata.parquet'); t = d.targets.fillna('').astype(str).str.strip()
    parts = t.apply(lambda s: [x.strip() for x in s.replace(';', ',').split(',') if x.strip()])
    one = parts.apply(len) == 1
    return dict(zip(d.drug[one].astype(str).str.strip(), parts[one].str[0]))


def tahoe_select(force=False):
    """Choose the cells to stream, from the 100M-row obs table (no expression read): per cell line, a seeded 18,400 of
    its DMSO_TF vehicle cells (context basals + the reference), and for every single-target drug at its HIGHEST dose
    (user decision), a seeded 450 of its cells (the effect builder samples 400 of those surviving cell QC)."""
    out = TAHOE_WORK / 'selection.parquet'
    if out.exists() and not force: return log(f'{out}: exists, skipped')
    TAHOE_WORK.mkdir(parents=True, exist_ok=True)
    import ast
    import pyarrow.parquet as pq
    drugs = tahoe_single_target_drugs()
    # the highest dose of each drug, from the 1,344-row sample table (doses are all uM)
    sm = pd.read_parquet(TAHOE / 'metadata/sample_metadata.parquet')
    parsed = sm.drugname_drugconc.astype(str).map(lambda s: ast.literal_eval(s)[0] if s.startswith('[(') else None)
    if {p[2] for p in parsed.dropna()} != {'uM'}: raise ValueError('Tahoe dose units are not all uM')
    sm['drug_s'] = sm.drug.astype(str).str.strip(); sm['conc'] = parsed.map(lambda p: float(p[1]) if p else np.nan)
    top = sm[sm.drug_s.isin(drugs)].groupby('drug_s').conc.max()
    keep_conds = set(sm.loc[sm.drug_s.isin(drugs) & sm.conc.eq(sm.drug_s.map(top)), 'drugname_drugconc'].astype(str))
    conc_of = dict(zip(sm.drugname_drugconc.astype(str), sm.conc))
    # fixed categories, so the per-row-group pieces concatenate as categoricals (~2 GB, not ~10 GB of strings)
    cats = dict(plate=sorted(sm.plate.astype(str).unique()), sample=sorted(sm['sample'].astype(str).unique()),
                drug=sorted(set(drugs) | {'DMSO_TF'}), drugname_drugconc=sorted(keep_conds | set(sm.loc[sm.drug_s.eq('DMSO_TF'), 'drugname_drugconc'].astype(str))),
                cell_name=sorted(pd.read_parquet(TAHOE / 'metadata/cell_line_metadata.parquet').cell_name.astype(str).unique()),
                pass_filter=['full', 'minimal'])
    pf = pq.ParquetFile(TAHOE / 'raw/obs_metadata.parquet'); parts = []
    cols = ['BARCODE_SUB_LIB_ID', 'plate', 'sample', 'drug', 'drugname_drugconc', 'cell_name', 'tscp_count', 'pass_filter']
    for i in range(pf.metadata.num_row_groups):
        t = pf.read_row_group(i, columns=cols).to_pandas()
        dr = t.drug.astype(str).str.strip(); dc = t.drugname_drugconc.astype(str)
        keep = (dr.eq('DMSO_TF') | dc.isin(keep_conds)).to_numpy()
        t = t[keep].copy(); t['drug'] = dr[keep].to_numpy()
        for c, cv in cats.items():
            v = t[c].astype(str)
            if not v.isin(cv).all(): raise ValueError(f'row group {i}: {c} values outside the metadata, e.g. {v[~v.isin(cv)].unique()[:3]}')
            t[c] = pd.Categorical(v, categories=cv)
        parts.append(t)
        if i % 16 == 0: log(f'  obs row group {i + 1} / {pf.metadata.num_row_groups}')
    o = pd.concat(parts, ignore_index=True); del parts
    o['conc_uM'] = o.drugname_drugconc.astype(str).map(conc_of).astype(float)
    is_dmso = o.drug.eq('DMSO_TF')
    rng = np.random.default_rng(TAHOE_SEED); picked = []
    for (line, drug), idx in o.groupby(['cell_name', 'drug'], sort=True, observed=True).indices.items():
        cap = TAHOE_DMSO_CAP if drug == 'DMSO_TF' else TAHOE_COND_CAP
        idx = np.sort(idx); picked.append(np.sort(rng.choice(idx, cap, replace=False)) if len(idx) > cap else idx)
    s = o.iloc[np.concatenate(picked)].copy()
    for c in cats: s[c] = s[c].astype(str)
    s['target_gene'] = np.where(s.drug.eq('DMSO_TF'), CONTROL, s.drug.map(drugs).fillna(''))
    s['guide_id'] = np.where(s.drug.eq('DMSO_TF'), 'DMSO_TF', s.drug + '_' + s.conc_uM.astype(str) + 'uM')
    if not s.BARCODE_SUB_LIB_ID.is_unique: raise ValueError('Tahoe barcodes are not unique')
    s.to_parquet(out, index=False)
    g = s.groupby('cell_name').drug
    log(f'selected {len(s):,} cells over {s.cell_name.nunique()} lines: {int(s.drug.eq("DMSO_TF").sum()):,} DMSO_TF, '
        f'{int((~s.drug.eq("DMSO_TF")).sum()):,} drug cells in {s[~s.drug.eq("DMSO_TF")].groupby(["cell_name", "drug"]).ngroups:,} line x drug conditions '
        f'({s.drug[~s.drug.eq("DMSO_TF")].nunique()} drugs at their highest dose); -> {out}')


def tahoe_stream(force=False, workers=6):
    """Stream every expression shard once (download -> keep the selected cells -> delete), appending each line's cells to
    its own CSR pieces; resumable (done.jsonl records each finished shard with the per-line file sizes). Then one raw
    file per cell line: lib/cells/tahoe_<line>.h5ad (target_gene = the drug's target, 'non-targeting' for DMSO_TF;
    guide_id = drug_<dose>uM; batch = plate)."""
    import concurrent.futures as cf
    import pyarrow.parquet as pq
    sel = pd.read_parquet(TAHOE_WORK / 'selection.parquet')
    lines = sorted(sel.cell_name.unique()); line_idx = {l: i for i, l in enumerate(lines)}
    want = dict(zip(sel.BARCODE_SUB_LIB_ID, sel.cell_name.map(line_idx)))
    gm = pd.read_parquet(TAHOE / 'metadata/gene_metadata.parquet').sort_values('token_id').reset_index(drop=True)
    tok2col = np.full(int(gm.token_id.max()) + 1, -1, np.int64); tok2col[gm.token_id.to_numpy()] = np.arange(len(gm))
    parts = TAHOE_WORK / 'parts'; parts.mkdir(parents=True, exist_ok=True); shards = TAHOE_WORK / 'shards'; shards.mkdir(exist_ok=True)
    done_path = TAHOE_WORK / 'done.jsonl'
    done, last = {}, None
    if done_path.exists():
        for line in done_path.read_text().splitlines():
            if line.strip(): r = json.loads(line); done[r['shard']] = r['sizes']; last = r['sizes']   # file order = append order
    files = {k: [parts / f'{i:02d}.{k}' for i in range(len(lines))] for k in ('data', 'idx', 'len', 'bc')}
    # crash-safe resume: cut every piece back to its size after the last shard recorded as finished
    for k, paths in files.items():
        for i, p in enumerate(paths):
            size = (last or {}).get(f'{i}.{k}', 0)
            with open(p, 'ab') as fh: fh.truncate(size)
    todo = [s for s in range(TAHOE_SHARDS) if s not in done]
    log(f'Tahoe stream: {len(want):,} selected cells in {len(lines)} lines; {len(todo)} of {TAHOE_SHARDS} shards to go')

    def fetch(s):
        dest = shards / f'{s:05d}.parquet'
        for attempt in range(5):
            r = subprocess.run(['curl', '-sL', '--fail', '--retry', '3', '-o', str(dest), TAHOE_URL.format(s)])
            if r.returncode == 0: return s, dest
            time.sleep(10 * (attempt + 1))
        raise RuntimeError(f'shard {s}: download failed')
    t0 = time.time(); n_cells = 0
    handles = {k: [open(p, 'a' if k == 'bc' else 'ab') for p in ps] for k, ps in files.items()}
    with cf.ThreadPoolExecutor(workers) as ex:
        futs = {ex.submit(fetch, s): s for s in todo[:workers * 2]}; nxt = workers * 2
        while futs:
            f = next(cf.as_completed(futs)); futs.pop(f); s, dest = f.result()
            if nxt < len(todo): futs[ex.submit(fetch, todo[nxt])] = todo[nxt]; nxt += 1
            t = pq.read_table(dest, columns=['genes', 'expressions', 'BARCODE_SUB_LIB_ID']); dest.unlink()
            bc = t['BARCODE_SUB_LIB_ID'].to_pylist(); li = np.array([want.get(b, -1) for b in bc])
            rows = np.flatnonzero(li >= 0)
            if len(rows):
                g = t['genes'].combine_chunks(); e = t['expressions'].combine_chunks()
                go = g.offsets.to_numpy(); gv = g.values.to_numpy(); ev = e.values.to_numpy()
                if not np.array_equal(go, e.offsets.to_numpy()): raise ValueError(f'shard {s}: genes / expressions misaligned')
                lens = np.diff(go); rid = np.repeat(np.arange(len(lens)), lens)        # row of every stored value, once
                # selected cells' values; tokens 0-2 are not genes (token 1 opens every cell with the value -2): dropped
                keep = (li[rid] >= 0) & (gv >= 3)
                rv = rid[keep]; cols = tok2col[gv[keep]]; vals = ev[keep]
                if (cols < 0).any(): raise ValueError(f'shard {s}: gene tokens outside gene_metadata')
                if not np.all(vals == np.rint(vals)) or (vals < 0).any(): raise ValueError(f'shard {s}: non-integer counts')
                nl = np.bincount(rv, minlength=len(lens))                                # kept values per cell
                order = np.argsort(li[rv], kind='stable')                                # by line; rows stay in order
                rv, cols, vals = rv[order], cols[order], vals[order]
                bounds = np.searchsorted(li[rv], np.arange(len(lines) + 1))
                for i in np.unique(li[rows]):
                    a, b = bounds[i], bounds[i + 1]; r = rows[li[rows] == i]
                    handles['data'][i].write(vals[a:b].astype(np.int32).tobytes())
                    handles['idx'][i].write(cols[a:b].astype(np.int32).tobytes())
                    handles['len'][i].write(nl[r].astype(np.int64).tobytes())
                    handles['bc'][i].write('\n'.join(bc[j] for j in r) + '\n')
                n_cells += len(rows)
            for hs in handles.values():
                for h in hs: h.flush()
            sizes = {f'{i}.{k}': p.stat().st_size for k, ps in files.items() for i, p in enumerate(ps) if p.exists()}
            with open(done_path, 'a') as fh: fh.write(json.dumps(dict(shard=s, sizes=sizes)) + '\n')
            done[s] = sizes
            if len(done) % 50 == 0:
                log(f'  {len(done):,} / {TAHOE_SHARDS} shards, {n_cells:,} cells kept this run ({(time.time() - t0) / 60:.0f} min)')
    for hs in handles.values():
        for h in hs: h.close()
    log(f'Tahoe stream complete: {len(done):,} shards')


def tahoe_finalize(force=False):
    """One raw file per cell line from the streamed pieces."""
    sel = pd.read_parquet(TAHOE_WORK / 'selection.parquet').set_index('BARCODE_SUB_LIB_ID')
    lines = sorted(sel.cell_name.unique())
    gm = pd.read_parquet(TAHOE / 'metadata/gene_metadata.parquet').sort_values('token_id').reset_index(drop=True)
    parts = TAHOE_WORK / 'parts'
    for i, line in enumerate(lines):
        name = tahoe_name(line)
        if (OUT / f'{name}.h5ad').exists() and not force: log(f'{name}: exists, skipped'); continue
        bc = (parts / f'{i:02d}.bc').read_text().split()
        lens = np.fromfile(parts / f'{i:02d}.len', dtype=np.int64)
        data = np.fromfile(parts / f'{i:02d}.data', dtype=np.int32); idx = np.fromfile(parts / f'{i:02d}.idx', dtype=np.int32)
        if len(bc) != len(lens) or int(lens.sum()) != len(data): raise ValueError(f'{line}: inconsistent streamed pieces')
        X = sp.csr_matrix((data, idx, np.concatenate([[0], np.cumsum(lens)])), shape=(len(bc), len(gm)))
        s = sel.loc[bc]
        missing = int((sel.cell_name == line).sum()) - len(bc)
        o = pd.DataFrame(dict(target_gene=s.target_gene.to_numpy(), guide_id=s.guide_id.to_numpy(), batch=s.plate.to_numpy(),
                              drug=s.drug.to_numpy(), conc_uM=s.conc_uM.to_numpy(), tahoe_sample=s['sample'].to_numpy(),
                              author_pass_filter=s.pass_filter.to_numpy(), author_tscp_count=s.tscp_count.to_numpy(),
                              source_cell_id=s.index.astype(str)), index=s.index.astype(str))
        write_raw(name, X, o, gm.gene_symbol.astype(str), gm.ensembl_id.astype(str),
                  dict(source='hf://datasets/tahoebio/Tahoe-100M (data/train-*.parquet, metadata/obs_metadata.parquet)',
                       dataset=f'Tahoe-100M {line}: DMSO_TF vehicle + single-target drugs at their highest dose', modality='drug',
                       background=str(line), lane='plate', controls='DMSO_TF (plate-matched vehicle)',
                       selection=f'per line <= {TAHOE_DMSO_CAP:,} DMSO_TF and <= {TAHOE_COND_CAP} cells per drug (seed {TAHOE_SEED})',
                       selected_not_found_in_shards=missing,
                       note='target_gene = the drug\'s single annotated target; guide_id = drug_<dose>uM; authors\' pass_filter kept, not applied'))
        for k in ('data', 'idx', 'len', 'bc'): (parts / f'{i:02d}.{k}').unlink()       # this line's pieces are in the file now


# ---------------------------------------------------------------- Jiang (re-included 2026-10-26: CRISPRko, on-target filtered)
JIANG_WORK = OUT / '_jiang_work'
JIANG_PATHWAYS = ('IFNB', 'IFNG', 'INS', 'TGFB', 'TNFA')


def jiang_symbol_table(force=False):
    """symbol -> Ensembl ID (one per symbol, '' when ambiguous) from every gene table here, for the conversion on vcc2027."""
    lookup = symbol_ids(); JIANG_WORK.mkdir(parents=True, exist_ok=True)
    rows = [(s, next(iter(ids)) if len(ids) == 1 else '') for s, ids in lookup.items()]
    pd.DataFrame(rows, columns=['symbol', 'gene_id']).to_csv(JIANG_WORK / 'symbol_ids.csv', index=False)
    log(f'{JIANG_WORK / "symbol_ids.csv"}: {len(rows):,} symbols ({sum(1 for _, g in rows if g):,} unambiguous)')


def jiang(force=False):
    """Per pathway (one Seurat object, exported by scripts/extract_jiang_counts.R into _jiang_work/<pathway>/), one raw file
    per cell line: lib/cells/jiang_<pathway>_<line>.h5ad. target_gene = metadata 'gene' ('NT' -> control), guide_id =
    'guide', batch = 'sample_ID' (the Parse sample well). Gene IDs from _jiang_work/symbol_ids.csv (symbols only in the source)."""
    ids = pd.read_csv(JIANG_WORK / 'symbol_ids.csv', keep_default_na=False)
    sym2id = dict(zip(ids.symbol.astype(str), ids.gene_id.astype(str)))
    for pw in JIANG_PATHWAYS:
        d = JIANG_WORK / pw
        if not (d / 'complete.txt').exists(): log(f'{pw}: not exported, skipped'); continue
        n_cells, n_genes, nnz = map(int, (d / 'shape.txt').read_text().split())
        genes = (d / 'genes.txt').read_text().splitlines()
        x = np.fromfile(d / 'x.bin', dtype=np.int32); i = np.fromfile(d / 'i.bin', dtype=np.int32); p = np.fromfile(d / 'p.bin', dtype=np.float64).astype(np.int64)
        if len(x) != nnz or len(i) != nnz or len(p) != n_cells + 1 or len(genes) != n_genes: raise ValueError(f'{pw}: inconsistent export')
        X = sp.csr_matrix((x, i, p), shape=(n_cells, n_genes)); del x, i
        m = pd.read_csv(d / 'metadata.csv.gz', index_col=0, low_memory=False)
        gene_ids = [sym2id.get(g, '') for g in genes]
        target = m['gene'].astype(str).where(m['gene'].astype(str) != 'NT', CONTROL)
        for line in sorted(m['cell_type'].astype(str).unique()):
            name = f'jiang_{pw.lower()}_{line.lower()}'
            if (OUT / f'{name}.h5ad').exists() and not force: log(f'{name}: exists, skipped'); continue
            rows = np.flatnonzero(m['cell_type'].astype(str).to_numpy() == line); s = m.iloc[rows]
            o = pd.DataFrame(dict(target_gene=target.iloc[rows].to_numpy(), guide_id=s['guide'].astype(str).to_numpy(),
                                  batch=s['sample_ID'].astype(str).to_numpy(), pathway=s['pathway'].astype(str).to_numpy(),
                                  author_replicate=s['Batch_info'].astype(str).to_numpy(), author_mixscale_score=s['mixscale_score'].to_numpy(),
                                  source_cell_id=s.index.astype(str)), index=s.index.astype(str))
            write_raw(name, X[rows], o, genes, gene_ids,
                      dict(source=str(V / f'lib/jiang/Seurat_object_{pw}_Perturb_seq.rds'), dataset=f'Jiang 2025 pathway Perturb-seq {pw}, {line}',
                           modality='CRISPRko', background=line, lane="sample_ID (Parse sample well)", controls="NT (non-targeting)",
                           genes_without_ensembl=int(sum(1 for g in gene_ids if not g)),
                           note='raw RNA counts from the Seurat object (read through its slots, no Seurat package); Parse split-pool, shallow (~3.9k UMI)'))


# ---------------------------------------------------------------- target symbols
def symbol_ids():
    """Symbol -> Ensembl IDs from every gene table we hold (older annotations included): Replogle / Nadig / H1 QC'd
    files and the scPerturb downloads. Resolves previous symbols (ATP5B, DGCR14, ...) through their Ensembl ID."""
    out = {}
    def add(sym, ids):
        for s, g in zip(map(str, sym), map(str, ids)):
            if g.startswith('ENSG'): out.setdefault(s, set()).add(g.split('.')[0])
    for n in ('replogle_k562_gw', 'replogle_rpe1', 'nadig_hepg2', 'nadig_jurkat', 'arc_h1'):
        p = V / f'lib/qc_cells/{n}_qc.h5ad'
        if not p.exists(): continue
        with h5py.File(p, 'r') as f: v = read_elem(f['var'])
        add(v['gene_name'], v.index) if 'gene_name' in v else add(v.index, v['gene_id'])
    for p in sorted((V / 'lib/scperturb').glob('*.h5ad')):
        with h5py.File(p, 'r') as f: v = read_elem(f['var'])
        c = next((c for c in ('ensembl_id', 'gene_id', 'ensemble_id') if c in v), None)
        if c: add(v.index, v[c])
    return out


def rebuilt_conversion(name, obs, f):
    """The uns['conversion'] record write_raw wrote, rebuilt from the file (only feng_gw_fitness lost it, 26 Oct)."""
    screen = {'feng_gw_fitness': 'GenomeWideScreen_FitnessGenes', 'feng_gw_nonfitness': 'GenomeWideScreen_NonFitnessGenes'}.get(name)
    if screen is None: raise ValueError(f'{name}: no conversion record and no rebuild rule')
    lab = obs['author_target'].astype(str) if 'author_target' in obs else obs['target_gene'].astype(str)
    X = f['X']
    return dict(source=str(FENG / f'{screen}_RNA-UMI-Counts.csv.gz'), dataset=f'Feng HipSci iPSC CRISPRi {screen}, donor lines and days pooled',
                modality='CRISPRi', background='HipSci iPSC (pooled donor lines)', lane='donor line',
                note='assigned cells only; donor_line and day kept in obs', converted='2026-10-26 02:06',
                script='scripts/convert_new_contexts.py', cells=int(len(obs)), genes=int(X.attrs['shape'][1]),
                controls=int((lab == CONTROL).sum()), targets=int(lab[lab != CONTROL].nunique()), nnz=int(X['indptr'][-1]),
                rebuilt='uns/conversion rebuilt 2026-10-26 after a failed symbol-fix write deleted it')


def fix_symbols(force=False):
    """Feng targets named with previous / non-HGNC symbols get the file's current symbol, so write_replogle_qc finds their
    Ensembl ID: 'HLA_A' -> 'HLA-A'; a previous symbol -> the file's symbol for the same Ensembl ID. Rewrites obs only (the
    original label stays in author_target); anything unresolved keeps its name and gets its own 'unresolved:' ID at QC."""
    lookup = symbol_ids(); log(f'symbol table: {len(lookup):,} symbols')
    for p in sorted(OUT.glob('feng_*.h5ad')):
        with h5py.File(p, 'r+') as f:
            obs = read_elem(f['obs']); var = read_elem(f['var']); uns = read_elem(f['uns'])
            if 'symbol_fixes' in uns.get('conversion', {}) and not force: log(f'{p.name}: already fixed, skipped'); continue
            sym = set(var.index.astype(str)); id2sym = {g.split('.')[0]: s for s, g in zip(var.index.astype(str), var['gene_id'].astype(str))}
            orig = obs['author_target'].astype(str) if 'author_target' in obs else obs['target_gene'].astype(str)
            fixes, unresolved = {}, []
            for t in sorted(set(orig) - {CONTROL} - sym):
                cand = t.replace('_', '-')
                if cand in sym: fixes[t] = cand; continue
                hit = sorted({id2sym[g] for g in lookup.get(t, ()) if g in id2sym})
                if len(hit) == 1: fixes[t] = hit[0]
                else: unresolved.append(t)
            obs['author_target'] = orig.to_numpy()
            obs['target_gene'] = orig.map(lambda t: fixes.get(t, t)).to_numpy()
            for c in ('target_gene', 'author_target'): obs[c] = obs[c].astype(str).astype('category')
            del f['obs']; write_elem(f, 'obs', obs)
            conv = dict(uns.get('conversion', {}))
            if not conv: conv = rebuilt_conversion(p.stem, obs, f)            # lost by the 26 Oct failed first attempt
            conv['symbol_fixes'] = json.dumps(fixes); conv['symbols_unresolved'] = json.dumps(unresolved)
            if 'conversion' in f['uns']: del f['uns']['conversion']
            write_elem(f['uns'], 'conversion', conv)
        log(f'{p.name}: {len(fixes)} targets renamed {dict(list(fixes.items())[:6])}; {len(unresolved)} unresolved {unresolved[:8]}')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('what', choices=['kolf', 'adamson', 'feng_targeted', 'feng_gw', 'fix_symbols', 'tahoe_select', 'tahoe_stream', 'jiang_symbol_table', 'jiang',
                                    'tahoe_finalize'])
    p.add_argument('--force', action='store_true')
    a = p.parse_args()
    dict(kolf=kolf, adamson=adamson, feng_targeted=feng_targeted, feng_gw=feng_gw, fix_symbols=fix_symbols, tahoe_select=tahoe_select,
         tahoe_stream=tahoe_stream, tahoe_finalize=tahoe_finalize, jiang_symbol_table=jiang_symbol_table, jiang=jiang)[a.what](a.force)
