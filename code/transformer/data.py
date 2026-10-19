"""Data for the vcc_2026_4 transformer: effect tables, evidence, basal features, labels and masks, gathered on the GPU.

Rows are (knockdown k, target line t). Targets: K562, RPE1, HepG2, Jurkat, HCT116, HEK293T in validation mode (+ H1 in
final mode); CD4 T (Rest, Stim8hr, Stim48hr) are evidence only. In validation mode H1 is not loaded at all.
Evidence for row (k, t): one slot per source s != t ([lfc, log se, is_measured, lfc_type one-hot, se_type one-hot,
p_de = 0, is_pie = 0]; the whole slot is 0 when s did not measure k) plus a PIE slot ([lfc_pred, 0, is_measured, 0..,
p_de, is_pie = 1]), always 0 when t is HCT116, HEK293T or H1 (PIE saw them).
Knockdown features: k's basal features in t, efficiency = median of k's own lfc over the evidence sources (never t;
+ N(0, eff_noise) in training only), efficiency_known. Targets: lfc_t(k, .), se_t, DE = BH 5% in t over measured genes
(two-sided normal p from lfc / se, k's own gene excluded). Every NaN becomes 0 with is_measured / loss_mask 0;
loss_mask is also 0 for k's own gene and for genes no loaded source measured.
"""
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import torch
from anndata.io import read_elem
from scipy.stats import norm

V = Path.home() / 'vcc'
SOURCES = [('replogle_k562_gw', 'K562'), ('replogle_rpe1', 'RPE1'), ('nadig_hepg2', 'HepG2'), ('nadig_jurkat', 'Jurkat'),
           ('orion_hct116', 'HCT116'), ('orion_hek293t', 'HEK293T'), ('marson_cd4_rest', 'CD4_Rest'),
           ('marson_cd4_stim8hr', 'CD4_Stim8hr'), ('marson_cd4_stim48hr', 'CD4_Stim48hr'), ('arc_h1', 'H1')]
SLOT = {line: i for i, (_, line) in enumerate(SOURCES)}
TARGETS = {'validation': ['K562', 'RPE1', 'HepG2', 'Jurkat', 'HCT116', 'HEK293T'],
           'final': ['K562', 'RPE1', 'HepG2', 'Jurkat', 'HCT116', 'HEK293T', 'H1']}
EVIDENCE_ONLY = ['CD4_Rest', 'CD4_Stim8hr', 'CD4_Stim48hr']
PIE_LEAKED = {'HCT116', 'HEK293T', 'H1'}
PIE_LINES = ['K562', 'RPE1', 'HepG2', 'Jurkat', 'A', 'B', 'C']
LINES = [line for _, line in SOURCES] + ['A', 'B', 'C']
LINE = {line: i for i, line in enumerate(LINES)}
BASAL_STEM = {**{line: src for src, line in SOURCES}, 'A': 'challenge_A', 'B': 'challenge_B', 'C': 'challenge_C'}
LFC_TYPES = ['mean_cpm_ratio', 'sum_ratio']
SE_TYPES = ['delta_cells', 'orion_cell_moments_delta', 'donor_moderated']
N_SLOTS = len(SOURCES) + 1                   # + PIE
PIE_SLOT = len(SOURCES)
SLOT_FEATURES = 10
BASAL_FEATURES = 5                            # detection_rate, detection_missing, pct_rank_S, log1p_cpm_S, is_measured
KD_FEATURES = BASAL_FEATURES + 2              # + efficiency, efficiency_known
FDR = 0.05
CHUNK = 2000


def kd_keys(obs):
    """Knockdown identity: the target's Ensembl ID, else 'symbol:<target_gene>'."""
    ids = obs['target_gene_id'].to_numpy(); sym = obs['target_gene'].astype(str).to_numpy()
    return np.array([str(t) if isinstance(t, str) and t.startswith('ENSG') else f'symbol:{s}' for t, s in zip(ids, sym)])


def panel(effects_dir):
    genes = pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str).to_numpy()
    with h5py.File(Path(effects_dir) / 'replogle_k562_gw.h5ad', 'r') as f:
        ids = read_elem(f['var'])['gene_id'].astype(str).to_numpy()
    return genes, ids


def bh_reject(p, alpha=FDR):
    n = len(p)
    if n == 0: return np.zeros(0, bool)
    o = np.argsort(p); q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.zeros(n, bool); out[o] = q < alpha
    return out


def de_full(effects_dir, source, cache_dir, kd_cols):
    """DE labels over the whole panel (rows x 18,533 bool), BH per row over measured genes, own gene excluded; cached."""
    cache = Path(cache_dir) / f'de_{source}.npz'
    path = Path(effects_dir) / f'{source}.h5ad'
    if cache.exists():
        z = np.load(cache)
        if int(z['mtime_ns']) == path.stat().st_mtime_ns:
            return np.unpackbits(z['bits'], axis=1, count=int(z['n_genes'])).astype(bool)
    with h5py.File(path, 'r') as f:
        n, G = f['layers/lfc'].shape; out = np.zeros((n, G), bool)
        for s in range(0, n, CHUNK):
            L = f['layers/lfc'][s:s + CHUNK].astype(np.float64); S = f['layers/se'][s:s + CHUNK].astype(np.float64)
            with np.errstate(invalid='ignore', divide='ignore'):
                P = 2 * norm.sf(np.abs(L / S))
            for i in range(len(L)):
                m = np.isfinite(P[i]); own = kd_cols[s + i]
                if own >= 0: m[own] = False
                cols = np.flatnonzero(m); out[s + i, cols] = bh_reject(P[i, cols])
    Path(cache_dir).mkdir(parents=True, exist_ok=True)
    np.savez(cache, bits=np.packbits(out, axis=1), n_genes=G, mtime_ns=np.int64(path.stat().st_mtime_ns))
    return out


def gene_ranking(effects_dir, cache_dir, lines):
    """Per panel gene: number of the given target lines measuring it, mean |lfc / se| over their rows; cached."""
    cache = Path(cache_dir) / f'gene_rank_{"_".join(lines)}.csv'
    if cache.exists():
        return pd.read_csv(cache, index_col=0)
    genes, _ = panel(effects_dir); n_lines = np.zeros(len(genes), int); zsum = np.zeros(len(genes)); zn = np.zeros(len(genes))
    for src, line in SOURCES:
        if line not in lines: continue
        with h5py.File(Path(effects_dir) / f'{src}.h5ad', 'r') as f:
            n = f['layers/lfc'].shape[0]; meas = np.zeros(len(genes), bool)
            for s in range(0, n, CHUNK):
                L = f['layers/lfc'][s:s + CHUNK]; S = f['layers/se'][s:s + CHUNK]
                ok = np.isfinite(L); meas |= ok.any(0)
                z = np.where(ok, np.abs(L / np.where(ok, S, 1)), 0); zsum += z.sum(0); zn += ok.sum(0)
        n_lines += meas
    rank = pd.DataFrame(dict(n_lines=n_lines, mean_abs_z=np.where(zn > 0, zsum / np.maximum(zn, 1), 0.0)),
                        index=pd.Index(genes, name='gene'))
    Path(cache_dir).mkdir(parents=True, exist_ok=True); rank.to_csv(cache)
    return rank


def select_genes(cfg, genes, targets_300):
    """Panel columns of the run's gene set: full panel, or our 300 targets + top genes up to n_genes."""
    if cfg['gene_set'] == 'full':
        return np.arange(len(genes))
    rank = gene_ranking(cfg['paths']['effects'], cfg['paths']['cache'], TARGETS['validation'])
    order = rank.sort_values(['n_lines', 'mean_abs_z'], ascending=False).index
    pos = {g: i for i, g in enumerate(genes)}
    chosen = [g for g in genes if g in targets_300]
    for g in order:
        if len(chosen) >= cfg['n_genes']: break
        if g not in targets_300: chosen.append(g)
    return np.sort(np.array([pos[g] for g in chosen]))


def basal_table(effects_dir, line, genes):
    """(18,533 x 5) basal features of a line from its baseline_v2 table, or the fallback from baseline_panel."""
    stem = BASAL_STEM[line]; v2 = Path(effects_dir) / f'{stem}_baseline_v2.csv'
    if v2.exists():
        b = pd.read_csv(v2, index_col=0).reindex(genes)
        det = b['detection_rate'].to_numpy(float); miss = (~np.isfinite(det)) | (b['detection_missing'].fillna(1).to_numpy() > 0)
        feats = np.stack([np.nan_to_num(det), miss.astype(float), np.nan_to_num(b['pct_rank_S'].to_numpy(float)),
                          np.nan_to_num(b['log1p_cpm_S'].to_numpy(float)), b['measured'].fillna(0).to_numpy(float)], 1)
        return feats.astype(np.float32), 'v2'
    b = pd.read_csv(Path(effects_dir) / f'{stem}_baseline_panel.csv', index_col=0).reindex(genes)
    cpm = b['control_mean_cpm'].to_numpy(float); det = b['detection_rate'].to_numpy(float); meas = np.isfinite(cpm)
    in_S = np.ones(len(genes), bool)
    for f in Path(effects_dir).glob('*_baseline.csv'):
        in_S &= pd.read_csv(f, index_col=0)['control_mean_cpm'].reindex(genes).notna().to_numpy()
    ref = np.sort(cpm[in_S & meas]); pct = np.zeros(len(genes))
    pct[meas] = np.searchsorted(ref, cpm[meas], side='right') / max(len(ref), 1)
    feats = np.stack([np.nan_to_num(det), (~np.isfinite(det)).astype(float), pct, np.zeros(len(genes)), meas.astype(float)], 1)
    return feats.astype(np.float32), 'fallback (no log CPM)'


class Bank:
    """Everything a batch needs, on one device: per-source lfc / se on the run's genes (float16), DE labels of the target
    sources, PIE tables, basal features, efficiencies and the row index (k, t)."""

    def __init__(self, cfg, device='cuda', verbose=True):
        self.cfg = cfg; self.mode = cfg['mode']; self.device = torch.device(device)
        E = Path(cfg['paths']['effects']); cache = cfg['paths']['cache']
        self.genes, self.panel_ids = panel(E)
        self.targets_300 = set(pd.read_csv(V / 'data/pert_counts.csv')['target_gene'].astype(str))
        self.gcols = select_genes(cfg, self.genes, self.targets_300); G = len(self.gcols)
        id_col = {g: i for i, g in enumerate(self.panel_ids) if g.startswith('ENSG')}
        sym_col = {g: i for i, g in enumerate(self.genes)}
        self.target_lines = TARGETS[self.mode]
        self.loaded = [line for _, line in SOURCES if self.mode == 'final' or line != 'H1']
        if verbose: print(f'mode: {self.mode} | targets: {self.target_lines} | evidence-only: {EVIDENCE_ONLY} | '
                          f'genes: {G} ({cfg["gene_set"]})' + ('' if 'H1' in self.loaded else ' | H1 not loaded'), flush=True)
        # knockdown vocabulary
        obs = {}
        for src, line in SOURCES:
            if line not in self.loaded: continue
            with h5py.File(E / f'{src}.h5ad', 'r') as f:
                obs[line] = read_elem(f['obs'])
        pie_paths = {line: E / f'pie_{line.lower() if line not in "ABC" else line}.h5ad' for line in PIE_LINES}
        want_pie = cfg['pie'] == 'on' or (cfg['pie'] == 'auto' and all(p.exists() for l, p in pie_paths.items() if l not in 'ABC'))
        if cfg['pie'] == 'on':
            missing = [str(p) for l, p in pie_paths.items() if l not in 'ABC' and not p.exists()]
            if missing: raise FileNotFoundError(f'pie=on but missing {missing}')
        self.pie_on = bool(want_pie)
        pie_obs = {}
        if self.pie_on:
            for line, p in pie_paths.items():
                if p.exists():
                    with h5py.File(p, 'r') as f: pie_obs[line] = read_elem(f['obs'])
        target_keys = {(self.panel_ids[sym_col[g]] if self.panel_ids[sym_col[g]].startswith('ENSG') else f'symbol:{g}')
                       for g in self.targets_300 if g in sym_col}
        keys = set(target_keys)
        for o in list(obs.values()) + list(pie_obs.values()): keys |= set(kd_keys(o))
        self.kd_vocab = np.array(sorted(keys)); kd_id = {k: i for i, k in enumerate(self.kd_vocab)}; K = len(self.kd_vocab)
        self.kd_symbol = np.array([k[7:] if k.startswith('symbol:') else '' for k in self.kd_vocab], dtype=object)
        kd_col = np.array([id_col.get(k, sym_col.get(k[7:], -1)) if not k.startswith('symbol:') else sym_col.get(k[7:], -1)
                           for k in self.kd_vocab])
        for i, k in enumerate(self.kd_vocab):
            if kd_col[i] >= 0 and not self.kd_symbol[i]: self.kd_symbol[i] = self.genes[kd_col[i]]
        self.kd_col = torch.as_tensor(kd_col, device=self.device)
        self.target_keys = sorted(target_keys)
        # source tables on the run's genes
        n_src = len(SOURCES)
        self.lfc, self.se, self.de = [None] * n_src, [None] * n_src, [None] * n_src
        src_row = np.full((K, n_src), -1); eff = np.full((K, n_src), np.nan, np.float32)
        self.lfc_type = np.zeros(n_src, int); self.se_type = np.zeros(n_src, int); self.n_de_bh = {}
        has_data = np.zeros(G, bool)
        for s, (src, line) in enumerate(SOURCES):
            if line not in self.loaded: continue
            o = obs[line]; k = kd_keys(o); rows = np.array([kd_id[x] for x in k])
            src_row[rows, s] = np.arange(len(k)); eff[rows, s] = o['efficiency'].to_numpy(np.float32)
            self.lfc_type[s] = LFC_TYPES.index(str(o['lfc_type'].iloc[0])); self.se_type[s] = SE_TYPES.index(str(o['se_type'].iloc[0]))
            L = np.empty((len(k), G), np.float16); S = np.empty((len(k), G), np.float16)
            with h5py.File(E / f'{src}.h5ad', 'r') as f:
                for a in range(0, len(k), CHUNK):
                    L[a:a + CHUNK] = f['layers/lfc'][a:a + CHUNK][:, self.gcols]
                    S[a:a + CHUNK] = f['layers/se'][a:a + CHUNK][:, self.gcols]
            has_data |= np.isfinite(L).any(0)
            self.lfc[s] = torch.from_numpy(L).to(self.device); self.se[s] = torch.from_numpy(S).to(self.device)
            if line in self.target_lines:
                de = de_full(E, src, cache, kd_col[rows])
                self.n_de_bh[line] = (de.sum(1), o['n_de_bh'].to_numpy())
                self.de[s] = torch.from_numpy(de[:, self.gcols]).to(self.device)
            if verbose: print(f'  {line:13s} {len(k):6d} knockdowns', flush=True)
        self.src_row = torch.as_tensor(src_row, device=self.device); self.eff = torch.as_tensor(eff, device=self.device)
        self.has_data = torch.as_tensor(has_data, device=self.device)
        # PIE
        self.pie_lfc, self.pie_p = {}, {}; pie_row = np.full((K, len(PIE_LINES)), -1)
        for j, line in enumerate(PIE_LINES):
            if line not in pie_obs: continue
            k = kd_keys(pie_obs[line]); pie_row[[kd_id[x] for x in k], j] = np.arange(len(k))
            with h5py.File(pie_paths[line], 'r') as f:
                self.pie_lfc[j] = torch.from_numpy(f['layers/lfc_pred'][:][:, self.gcols].astype(np.float16)).to(self.device)
                self.pie_p[j] = torch.from_numpy(f['layers/p_de'][:][:, self.gcols].astype(np.float16)).to(self.device)
        self.pie_row = torch.as_tensor(pie_row, device=self.device)
        self.line_pie = torch.as_tensor([PIE_LINES.index(l) if (self.pie_on and l in PIE_LINES and l not in PIE_LEAKED
                                                                and PIE_LINES.index(l) in self.pie_lfc) else -1
                                         for l in LINES], device=self.device)
        self.line_slot = torch.as_tensor([SLOT.get(l, -1) for l in LINES], device=self.device)
        # basal features over the whole panel, per line
        basal, self.basal_kind = [], {}
        for line in LINES:
            try:
                b, kind = basal_table(E, line, self.genes)
            except FileNotFoundError:
                b, kind = np.zeros((len(self.genes), BASAL_FEATURES), np.float32), 'missing'
            basal.append(b); self.basal_kind[line] = kind
        self.basal = torch.as_tensor(np.stack(basal), device=self.device)
        self.gcols_t = torch.as_tensor(self.gcols, device=self.device)
        self.loaded_slots = torch.as_tensor([SLOT[l] for l in self.loaded], device=self.device)
        # rows (k, t) of the target lines, and the gene-level holdout
        rk, rl, rt = [], [], []
        for line in self.target_lines:
            s = SLOT[line]; r = np.flatnonzero(src_row[:, s] >= 0)
            rk.append(r); rl.append(np.full(len(r), LINE[line])); rt.append(src_row[r, s])
        self.row_kd = np.concatenate(rk); self.row_line = np.concatenate(rl); self.row_t = np.concatenate(rt)
        uniq = np.unique(self.kd_vocab[self.row_kd]); rng = np.random.default_rng(cfg['seed'])
        held = set(rng.permutation(uniq)[:int(round(cfg['holdout_frac'] * len(uniq)))])
        self.heldout_keys = held
        is_held = np.array([self.kd_vocab[k] in held for k in self.row_kd])
        self.train_rows = np.flatnonzero(~is_held); self.val_rows = np.flatnonzero(is_held)
        self._rows_t = (torch.as_tensor(self.row_kd, device=self.device), torch.as_tensor(self.row_line, device=self.device),
                        torch.as_tensor(self.row_t, device=self.device))
        if verbose:
            print(f'rows: {len(self.row_kd)} (train {len(self.train_rows)}, held-out genes {len(self.val_rows)}; '
                  f'{len(held)} of {len(uniq)} target genes held out) | PIE evidence: {"on" if self.pie_on else "off"} | '
                  f'basal: {sorted(set(self.basal_kind.values()))}', flush=True)

    # ------------------------------------------------------------------ batches
    def rows_for(self, line, keys):
        """Prediction rows for a line without its own table use (kd id, line, -1)."""
        kd_id = {k: i for i, k in enumerate(self.kd_vocab)}
        kd = torch.as_tensor([kd_id[k] for k in keys], device=self.device)
        return kd, torch.full_like(kd, LINE[line]), torch.full_like(kd, -1)

    def batch(self, rows=None, kd=None, line=None, trow=None, train=False, generator=None):
        """Tensors for a batch, from row ids of this bank (rows) or explicit (kd, line, trow) for prediction."""
        if rows is not None:
            r = torch.as_tensor(rows, device=self.device)
            kd, line, trow = self._rows_t[0][r], self._rows_t[1][r], self._rows_t[2][r]
        B, G = len(kd), len(self.gcols); dev = self.device
        t_slot = self.line_slot[line]
        feats = torch.zeros(B, G, N_SLOTS, SLOT_FEATURES, device=dev)
        wsum = torch.zeros(B, G, device=dev); wl = torch.zeros(B, G, device=dev)
        for s in self.loaded_slots.tolist():
            r = self.src_row[kd, s]; avail = (r >= 0) & (t_slot != s)
            idx = r.clamp(min=0)
            L = self.lfc[s][idx].float(); S = self.se[s][idx].float()
            m = avail[:, None] & torch.isfinite(L) & torch.isfinite(S)
            L = torch.where(m, L, 0.0); S = torch.where(m, S, 1.0); mf = m.float()
            f = feats[:, :, s]
            f[..., 0] = L; f[..., 1] = torch.log(S) * mf; f[..., 2] = mf
            f[..., 3 + self.lfc_type[s]] = mf; f[..., 5 + self.se_type[s]] = mf
            w = mf / S.square(); wsum += w; wl += w * L
        pie_line = self.line_pie[line]
        if self.pie_on and (pie_line >= 0).any():
            for j in torch.unique(pie_line[pie_line >= 0]).tolist():
                sel = pie_line == j; r = self.pie_row[kd[sel], j]; ok = r >= 0; idx = r.clamp(min=0)
                L = self.pie_lfc[j][idx].float(); P = self.pie_p[j][idx].float()
                m = ok[:, None] & torch.isfinite(L) & torch.isfinite(P); mf = m.float()
                f = torch.zeros(int(sel.sum()), G, SLOT_FEATURES, device=dev)
                f[..., 0] = torch.where(m, L, 0.0); f[..., 2] = mf; f[..., 8] = torch.where(m, P, 0.0); f[..., 9] = mf
                feats[sel, :, PIE_SLOT] = f
        pw_mean = torch.where(wsum > 0, wl / wsum.clamp(min=1e-12), 0.0)
        slot_mask = feats[..., 2] > 0
        # knockdown features
        e = self.eff[kd].clone()                                   # (B, sources): k's own lfc in each source
        own_line = torch.nonzero(t_slot >= 0).squeeze(1)
        e[own_line, t_slot[own_line]] = torch.nan                  # never t itself
        eff = torch.nanmedian(e, dim=1).values; known = torch.isfinite(eff)
        eff = torch.where(known, eff, 0.0)
        if train and self.cfg['train']['eff_noise'] > 0:
            eff = eff + known.float() * self.cfg['train']['eff_noise'] * torch.randn(B, device=dev, generator=generator)
        kcol = self.kd_col[kd]
        kd_basal = torch.where((kcol >= 0)[:, None], self.basal[line, kcol.clamp(min=0)], 0.0)
        kd_feats = torch.cat([kd_basal, eff[:, None], known.float()[:, None]], 1)
        basal = self.basal[line][:, self.gcols_t]
        out = dict(kd=kd, line=line, slot_feats=feats, slot_mask=slot_mask, pw_mean=pw_mean, basal=basal,
                   kd_feats=kd_feats, kd_col=kcol, gene_cols=self.gcols_t)
        # targets
        if (trow >= 0).any():
            lfc_t = torch.full((B, G), torch.nan, device=dev); se_t = torch.ones(B, G, device=dev)
            de = torch.zeros(B, G, dtype=torch.bool, device=dev)
            for s in torch.unique(t_slot[trow >= 0]).tolist():
                sel = (t_slot == s) & (trow >= 0)
                lfc_t[sel] = self.lfc[s][trow[sel]].float(); se_t[sel] = self.se[s][trow[sel]].float()
                de[sel] = self.de[s][trow[sel]]
            measured = torch.isfinite(lfc_t) & torch.isfinite(se_t)
            own = self.gcols_t[None, :] == kcol[:, None]
            loss_mask = measured & ~own & self.has_data[None, :]
            out.update(lfc_t=torch.where(measured, lfc_t, 0.0), se_t=torch.where(measured, se_t, 1.0),
                       de=de & loss_mask, loss_mask=loss_mask)
        return out
