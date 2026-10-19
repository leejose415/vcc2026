"""PIE gene features on the challenge's 18,533-gene panel (lib/pie_features/). Reads the PIE_sources files, changes nothing.

Sources (arcinstitute/PIE_sources, the commit pinned by PIE's replogle_xdataset checkpoint): esm2, string_space,
depmap_gene_effect, gene_text. Each panel gene is matched to a source row by, in order:
  1. symbol: the panel symbol is a source key;
  2. alias:  the source's aliases.yaml maps the panel symbol to a key (alias -> current symbol);
  3. ensembl: another symbol with the same Ensembl ID is a key, directly or through aliases.yaml. Symbol <-> ID pairs come
     from the gene tables of lib/qc_cells/* and lib/effects_v4 and, with --hgnc, from the HGNC complete set (current
     symbol per ID; a panel symbol with no ID gets one through HGNC previous symbols when that is unambiguous).
Unmatched genes get a zero row and mask False. Outputs per source: <source>.npy (18,533 x d, float16); plus
coverage_mask.npy (18,533 x 4 bool, column order in sources.json), coverage.csv (per gene, the match per source)
and sources.json (commit, dims, counts).

Usage: python scripts/export_pie_features.py --sources <PIE_sources dir> [--hgnc hgnc_complete_set.txt]
"""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import yaml
from anndata.io import read_elem

V = Path.home() / 'vcc'
OUT = V / 'lib/pie_features'
SOURCES = ['esm2', 'string_space', 'depmap_gene_effect', 'gene_text']
COMMIT = 'fb624a5158a7f37dfa13a015cdacadb985549d26'


def symbol_ids():
    """(symbol -> Ensembl ID, Ensembl ID -> set of symbols) from every local gene table that has both."""
    s2i, i2s = {}, {}
    files = sorted((V / 'lib/qc_cells').glob('*.h5ad')) + [V / 'lib/effects_v4/replogle_k562_gw.h5ad']
    for p in files:
        with h5py.File(p, 'r') as f:
            var = read_elem(f['var'])
        if 'gene_name' in var:
            pairs = zip(var['gene_name'].astype(str), var.index.astype(str))
        elif 'gene_id' in var:
            pairs = zip(var.index.astype(str), var['gene_id'].astype(str))
        else:
            continue
        for s, g in pairs:
            if not g.startswith('ENSG'):
                continue
            s2i.setdefault(s, g); i2s.setdefault(g, set()).add(s)
    return s2i, i2s


def add_hgnc(path, s2i, i2s):
    """HGNC: ID -> current symbol (and previous symbols); previous symbol -> ID when it names one gene only."""
    h = pd.read_csv(path, sep='\t', usecols=['symbol', 'prev_symbol', 'ensembl_gene_id'], dtype=str).dropna(subset=['ensembl_gene_id'])
    prev = {}
    for sym, prevs, gid in h.itertuples(index=False):
        i2s.setdefault(gid, set()).add(sym); s2i.setdefault(sym, gid)
        for p in str(prevs).split('|') if isinstance(prevs, str) else []:
            i2s[gid].add(p); prev.setdefault(p, set()).add(gid)
    for p, gids in prev.items():
        if len(gids) == 1: s2i.setdefault(p, next(iter(gids)))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--sources', required=True, type=Path); ap.add_argument('--hgnc', type=Path)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    panel = pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str).to_numpy()
    targets = set(pd.read_csv(V / 'data/pert_counts.csv')['target_gene'].astype(str))
    with h5py.File(V / 'lib/effects_v4/replogle_k562_gw.h5ad', 'r') as f:
        panel_id = read_elem(f['var'])['gene_id'].astype(str).to_numpy()        # '' where no source maps the symbol
    s2i, i2s = symbol_ids()
    if a.hgnc: add_hgnc(a.hgnc, s2i, i2s)
    mask = np.zeros((len(panel), len(SOURCES)), bool); how = {}; info = {}
    for k, name in enumerate(SOURCES):
        d = a.sources / name
        meta = json.loads((d / 'meta.json').read_text()); keys = {s: i for i, s in enumerate(meta['keys'])}
        aliases = yaml.safe_load((d / 'aliases.yaml').read_text()) if (d / 'aliases.yaml').exists() else {}
        aliases = {str(x): str(y) for x, y in (aliases or {}).items()}
        E = np.load(d / 'embeddings.npy', mmap_mode='r')
        assert E.shape == (len(keys), meta['dim']), (name, E.shape)

        def find(sym):
            if sym in keys: return keys[sym]
            if sym in aliases and aliases[sym] in keys: return keys[aliases[sym]]
            return None
        rows = np.full(len(panel), -1); method = np.full(len(panel), 'none', dtype=object)
        for j, s in enumerate(panel):
            if s in keys:
                rows[j], method[j] = keys[s], 'symbol'
            elif s in aliases and aliases[s] in keys:
                rows[j], method[j] = keys[aliases[s]], 'alias'
            else:
                gid = panel_id[j] or s2i.get(s, '')
                for other in sorted(i2s.get(gid, ())):
                    r = find(other)
                    if r is not None:
                        rows[j], method[j] = r, 'ensembl'; break
        ok = rows >= 0
        X = np.zeros((len(panel), meta['dim']), np.float16)
        X[ok] = np.asarray(E[rows[ok]], dtype=np.float32).astype(np.float16)
        np.save(OUT / f'{name}.npy', X)
        mask[:, k] = ok; how[name] = method
        tmask = np.array([s in targets for s in panel])
        info[name] = dict(dim=int(meta['dim']), source_dtype=meta['dtype'], source_keys=len(keys), index=meta['index'],
                          panel_covered=int(ok.sum()), panel_genes=len(panel),
                          by_method={m: int((method == m).sum()) for m in ('symbol', 'alias', 'ensembl', 'none')},
                          targets_covered=int((ok & tmask).sum()), targets_on_panel=int(tmask.sum()),
                          targets_missing=sorted(panel[tmask & ~ok].tolist()),
                          finite=bool(np.isfinite(X[ok].astype(np.float32)).all()))
        print(f"{name:20s} d={meta['dim']:5d}  panel {ok.sum():,}/{len(panel):,}  {info[name]['by_method']}  "
              f"targets {info[name]['targets_covered']}/{info[name]['targets_on_panel']}", flush=True)
    np.save(OUT / 'coverage_mask.npy', mask)
    pd.DataFrame({n: how[n] for n in SOURCES}, index=pd.Index(panel, name='gene')).to_csv(OUT / 'coverage.csv')
    off_panel = sorted(targets - set(panel))
    (OUT / 'sources.json').write_text(json.dumps(dict(repo='arcinstitute/PIE_sources', commit=COMMIT, columns=SOURCES,
                                                      id_map=['lib/qc_cells + lib/effects_v4 gene tables'] + (['HGNC complete set'] if a.hgnc else []),
                                                      dtype='float16', rows='data/gene_names.csv order', sources=info,
                                                      targets_off_panel=off_panel), indent=1) + '\n')
    print(f'targets not on the panel: {off_panel}')


if __name__ == '__main__':
    main()
