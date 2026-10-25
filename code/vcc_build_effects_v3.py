"""Build the v3 inputs from our public data: effect tables, control profiles and contexts.json.

  orion       Orion per-knockdown sums (<line>_pseudobulk.h5ad)  -> effect folder (LFC + cell-level SE) + control profile
  h1          Arc H1 sums (h1_pseudobulk.h5ad, same layout)      -> same as orion (held-out validation context)
  tcell       T-cell effect table (TCELL_<cond>_effects.h5ad)   -> effect folder (LFC + donor-level SE) + control profile
  make-config effect folders                                     -> contexts.json for vcc_transfer_model_v3.py

Every effect folder also gets context.json (name / background / platform), so make-config needs no typing.

Colab
  !python vcc_build_effects_v3.py orion --pseudobulk /content/drive/MyDrive/vcc/orion/HCT116_pseudobulk.h5ad \
        --genes "$DATA/gene_names.csv" --name ORION_HCT116 --background HCT116 --out $LIB/effects/ORION_HCT116
  !python vcc_build_effects_v3.py tcell --effects $LIB/tcell/effects/TCELL_Stim48hr_effects.h5ad \
        --genes "$DATA/gene_names.csv" --name TCELL_Stim48hr --priority 1 --out $LIB/effects/TCELL_Stim48hr
  !python vcc_build_effects_v3.py make-config --genes "$DATA/gene_names.csv" --out $LIB/contexts.json \
        --dirs $LIB/effects/ORION_HCT116 $LIB/effects/ORION_HEK293T $LIB/effects/TCELL_Stim48hr $LIB/effects/TCELL_Rest
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import pandas as pd
import h5py
import numpy as np
import scipy.sparse as sp

from vcc_common import CONTROL, alignment, read_names, save_json, save_profile, write_effects, moments, validate_counts

try:
    from anndata.io import read_elem, sparse_dataset
except ImportError:
    from anndata.experimental import read_elem, sparse_dataset

PSEUDO = 1e-9  # evaluator epsilon (cell-eval2 vcc2026); not used for modelling labels
# Modelling labels add this many CPM to both sides of every ratio. With the evaluator's 1e-9,
# a low-count gene with zero counts in a ~100-cell group gets LFC near -30 with a small
# delta-method SE: it defeats the normal-prior shrinkage and passes the expression gate in
# lines where the gene is expressed. The evaluator only scores genes above 5 CPM in controls,
# where 0.5 CPM moves an LFC by at most ~0.13; below that, noise collapses toward zero while
# genuine induction from ~0 CPM survives (e.g. 0 -> 10 CPM is still +4.4).
LABEL_PSEUDOCOUNT = 0.5


def mean_variance(mean, var_of_mean, unit):
    """Variance of a group's mean CPM, floored at Poisson noise plus one pseudo-count.

    `unit` is the CPM one count adds to the group mean (1e6 / group's total counts). An
    all-zero group has zero empirical variance, which made its label look precise; the
    floor (mean + unit) * unit keeps it as uncertain as one unseen count. Overdispersed
    genes have empirical variance above Poisson, so the floor rarely binds for them.
    """
    return np.maximum(var_of_mean, (np.maximum(mean, 0) + unit) * unit)


def label_pseudocount(a):
    pc = float(getattr(a, "pseudocount", LABEL_PSEUDOCOUNT))
    if not pc > 0: raise ValueError("label pseudocount must be positive CPM")
    return pc


def reader(f, key):
    """rows(list of int) -> dense float64 block, for dense or sparse h5ad matrices."""
    node = f[key]
    if isinstance(node, h5py.Group):
        ds = sparse_dataset(node)
        return lambda r: np.asarray(ds[np.asarray(r)].toarray(), dtype=np.float64)
    return lambda r: np.asarray(node[np.asarray(r), :], dtype=np.float64)


def names_of(df, col="target_gene"):
    return np.asarray(df[col] if col in df.columns else df.index, dtype=str)


def context_file(out, name, background, platform, priority, kind):
    save_json(Path(out) / "context.json", dict(name=name, background=background, platform=platform,
                                             source_priority=priority, kind=kind))


# ---------------------------------------------------------------- Orion / H1 (per-knockdown sums)
def cmd_sums(a, kind):
    genes = read_names(a.genes, "gene_name"); out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    pc = label_pseudocount(a)
    with h5py.File(a.pseudobulk, "r") as f:
        obs = read_elem(f["obs"]); var = read_elem(f["var"])
        control = a.control
        if control is None:
            control = str(read_elem(f["uns/control_label"])) if "uns/control_label" in f else CONTROL
        names = names_of(obs); n = obs["n_cells"].values.astype(np.float64)
        src = alignment(var.index.astype(str), genes); ok = src >= 0
        if "measured" in var.columns:
            meas = np.zeros(len(genes), bool); meas[ok] = var["measured"].values.astype(bool)[src[ok]]
        else:
            meas = ok.copy()
        rc, r1, r2, re = (reader(f, f"layers/{k}") for k in ("sum_counts", "sum_cpm", "sum_cpm2", "n_expr"))
        controls = np.flatnonzero(names == control)
        if len(controls) != 1:
            raise ValueError("sums input needs one aggregated control row")
        ci = int(controls[0])

        def aligned(block):
            outb = np.full((block.shape[0], len(genes)), np.nan); outb[:, ok] = block[:, src[ok]]; return outb
        c_cnt, c_cpm, c_cpm2, c_det = (aligned(r([ci]))[0] for r in (rc, r1, r2, re))
        nc = n[ci]
        c_pool = c_cnt / np.nansum(np.where(meas, c_cnt, 0)) * 1e6
        c_mean = c_cpm / nc; c_var = np.maximum(c_cpm2 / nc - c_mean ** 2, 0)
        c_unit = 1e6 / max(float(np.nansum(np.where(meas, c_cnt, 0))), 1.0)
        c_serel2 = mean_variance(c_mean, c_var / nc, c_unit) / (c_mean + pc) ** 2
        keep = np.flatnonzero((names != control) & (n >= a.min_cells))
        # Check EACH guide/row before selecting among successful guides.
        keep, qc = filter_sum_guides(keep, names, rc, r1, aligned, n, genes,
                                    c_pool, c_mean, meas, a)
        pd.DataFrame(qc).to_csv(out / "guide_qc.csv", index=False)
        targets = names[keep]
        # duplicate names (e.g. several guides): keep the row with most cells
        order = np.lexsort((-n[keep], targets)); first = np.r_[True, targets[order][1:] != targets[order][:-1]]
        keep = np.sort(keep[order][first]); targets = names[keep]
        print(f"{kind}: {len(keep):,} knockdowns with >= {a.min_cells} cells; control '{control}' with {int(nc):,} cells; "
              f"{int(meas.sum()):,}/{len(genes):,} genes measured")

        def rows():
            for b in range(0, len(keep), a.chunk):
                r = keep[b:b + a.chunk]
                cnt, cpm, cpm2 = (aligned(x(r)) for x in (rc, r1, r2))
                nk = n[r][:, None]
                pool = cnt / np.maximum(np.nansum(np.where(meas, cnt, 0), 1, keepdims=True), 1) * 1e6
                mean = cpm / nk
                effect_base = pool if a.lfc_estimator == "pooled" else mean
                control_base = c_pool if a.lfc_estimator == "pooled" else c_mean
                lfc = np.log2((effect_base + pc) / (control_base + pc))
                mean = cpm / nk; var_ = np.maximum(cpm2 / nk - mean ** 2, 0)
                unit = 1e6 / np.maximum(np.nansum(np.where(meas, cnt, 0), 1, keepdims=True), 1)
                se = np.sqrt(mean_variance(mean, var_ / nk, unit) / (mean + pc) ** 2 + c_serel2) / np.log(2)
                lfc[:, ~meas] = np.nan; se[:, ~meas] = np.nan
                for k in range(len(r)):
                    yield b + k, lfc[k], se[k]
                print(f"  {min(b + a.chunk, len(keep)):,}/{len(keep):,}", end="\r", flush=True)
        if not len(targets): raise ValueError("no successful target groups after guide QC")
        write_effects(out, genes, targets, rows(), dict(
            kind=kind, source=str(Path(a.pseudobulk).resolve()), control_label=control, min_cells=a.min_cells,
            lfc=f"log2 ratio of {a.lfc_estimator} CPM, pseudocount {pc} CPM", label_pseudocount=pc,
            label_estimator=a.lfc_estimator, baseline_convention="ratio_of_sums_cpm",
            estimator_matches_official=False, guide_qc="drop measured on-target LFC >= threshold before choosing largest group",
            se_kind="cell_level (delta method from per-cell CPM spread, Poisson one-count floor; overconfident vs replicate SE)", legacy=False),
            len(keep))
    prof = dict(mean_cpm=np.where(meas, c_mean, np.nan).astype(np.float32),
                pooled_cpm=np.where(meas, c_pool, np.nan).astype(np.float32),
                detection=np.where(meas, c_det / nc, np.nan).astype(np.float32),
                var_cpm=np.where(meas, c_var, np.nan).astype(np.float32),
                se_mean=np.where(meas, np.sqrt(c_var / nc), np.nan).astype(np.float32),
                measured=meas, n_cells=np.int64(nc), legacy_baseline=np.bool_(False))
    prof["log_mean"] = np.log1p(prof["pooled_cpm"])
    save_profile(out / "control_profile.npz", genes, prof, dict(source=str(a.pseudobulk), kind=kind, baseline="controls row"))
    context_file(out, a.name, a.background, a.platform, a.priority, kind)
    print(f"\nwrote {out}  (lfc.npy, se.npy, control_profile.npz, context.json)")


# ---------------------------------------------------------------- T cells (already an effect table)
def cmd_tcell(a):
    genes = read_names(a.genes, "gene_name"); out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    with h5py.File(a.effects, "r") as f:
        obs = read_elem(f["obs"]); var = read_elem(f["var"])
        names = names_of(obs); src = alignment(var.index.astype(str), genes); ok = src >= 0
        meas = np.zeros(len(genes), bool)
        meas[ok] = var["measured"].values.astype(bool)[src[ok]] if "measured" in var.columns else True
        ctrl = np.full(len(genes), np.nan); ctrl[ok] = var["ctrl_cpm"].values.astype(float)[src[ok]]
        rx = reader(f, "X"); rs = reader(f, "layers/se") if "layers/se" in f else None
        ncell = obs["n_cells"].values if "n_cells" in obs.columns else np.ones(len(names))
        order = np.lexsort((-ncell, names)); first = np.r_[True, names[order][1:] != names[order][:-1]]
        candidates = np.flatnonzero(names != CONTROL)
        vpos = {str(g): i for i, g in enumerate(var.index)}
        own_parts = []
        for start in range(0, len(candidates), a.chunk):
            part = candidates[start:start+a.chunk]; values = rx(part)
            own_parts.extend([values[i, vpos[t]] if t in vpos else np.nan for i, t in enumerate(names[part])])
        own = np.asarray(own_parts)
        passed = ~np.isfinite(own) | (own < a.min_ontarget_lfc)
        pd.DataFrame(dict(row=candidates, target=names[candidates], ontarget_lfc=own,
                          passed=passed, measured=np.isfinite(own))).to_csv(out / "guide_qc.csv", index=False)
        candidates = candidates[passed]
        order = np.lexsort((-ncell[candidates], names[candidates]))
        ordered = candidates[order]
        first = np.r_[True, names[ordered][1:] != names[ordered][:-1]] if len(ordered) else np.array([], bool)
        keep = np.sort(ordered[first]); targets = names[keep]
        print(f"tcell: {len(keep):,} knockdowns, {int(meas.sum()):,}/{len(genes):,} genes measured, SE: {rs is not None}")

        def rows():
            for b in range(0, len(keep), a.chunk):
                r = keep[b:b + a.chunk]
                L = np.full((len(r), len(genes)), np.nan); S = np.full_like(L, np.nan)
                L[:, ok] = rx(r)[:, src[ok]]
                if rs is not None:
                    S[:, ok] = rs(r)[:, src[ok]]
                L[:, ~meas] = np.nan; S[:, ~meas] = np.nan
                for k in range(len(r)):
                    yield b + k, L[k], S[k]
                print(f"  {min(b + a.chunk, len(keep)):,}/{len(keep):,}", end="\r", flush=True)
        if not len(targets): raise ValueError("no successful target groups after guide QC")
        write_effects(out, genes, targets, rows(), dict(
            kind="tcell", source=str(Path(a.effects).resolve()),
            lfc="replicate-paired log2 ratio of pseudobulk CPM (tcell_build_effects.py)",
            se_kind="donor_level (limma-style moderated)", legacy=rs is None), len(keep))
    pooled = np.where(meas, ctrl, np.nan).astype(np.float32)
    prof = dict(mean_cpm=pooled, pooled_cpm=pooled, log_mean=np.log1p(pooled),
                detection=np.full(len(genes), np.nan, np.float32), var_cpm=np.full(len(genes), np.nan, np.float32),
                se_mean=np.full(len(genes), np.nan, np.float32), measured=meas, n_cells=np.int64(0),
                legacy_baseline=np.bool_(True))     # baseline from pseudobulk, not single cells
    save_profile(out / "control_profile.npz", genes, prof, dict(source=str(a.effects), kind="tcell",
                 baseline="var['ctrl_cpm'] (pseudobulk; mean_cpm == pooled_cpm)"))
    context_file(out, a.name, a.background, a.platform, a.priority, "tcell")
    print(f"\nwrote {out}")


def filter_sum_guides(keep, names, rc, r1, aligned, n, genes, c_pool, c_mean, meas, a):
    gp = {g: i for i, g in enumerate(genes)}; retained, qc = [], []; pc = label_pseudocount(a)
    for start in range(0, len(keep), a.chunk):
        r = keep[start:start+a.chunk]
        cnt, cpm = aligned(rc(r)), aligned(r1(r))
        pool = cnt/np.maximum(np.nansum(np.where(meas, cnt, 0), axis=1, keepdims=True), 1)*1e6
        mean = cpm/n[r, None]
        base, ctrl = (pool, c_pool) if a.lfc_estimator == "pooled" else (mean, c_mean)
        for i, row in enumerate(r):
            j = gp.get(names[row])
            own = np.log2((base[i, j]+pc)/(ctrl[j]+pc)) if j is not None and meas[j] else np.nan
            passed = not np.isfinite(own) or own < a.min_ontarget_lfc
            qc.append(dict(row=int(row), target=names[row], ontarget_lfc=float(own),
                           measured=bool(np.isfinite(own)), passed=bool(passed)))
            if passed: retained.append(row)
    return np.asarray(retained, int), qc


def cmd_cells(a):
    """Raw Replogle/H1/etc. cells: guide QC, official 2026 arithmetic CPM labels.

    Normalization is per cell to 1e6 on the aligned gene axis, matching
    cell-eval2 vcc2026's arithmetic CPM convention. Labels add LABEL_PSEUDOCOUNT CPM
    (not the evaluator's 1e-9) so unscored low-count genes cannot produce huge labels.
    Delta-method SE remains cell-level and is identified as such in metadata.
    """
    genes = read_names(a.genes, "gene_name"); out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    pc = label_pseudocount(a)
    h = ad.read_h5ad(a.cells, backed="r")
    try:
        source_genes = h.var[a.gene_col].astype(str).values if a.gene_col else h.var_names.astype(str)
        if len(set(source_genes)) != len(source_genes):
            raise ValueError("raw source gene symbols must be unique; aggregate duplicate symbols first")
        idx = alignment(source_genes, genes); ok = idx >= 0
        labels = h.obs[a.obs_col].astype(str).values
        if a.target_split:
            labels = np.array([x.split(a.target_split)[0] for x in labels])
        control_values = a.control.split(',')
        ci = np.flatnonzero(np.isin(labels, control_values))
        if not len(ci): raise ValueError("no controls; set --control to the source's exact control label")
        matrix = h.layers[a.layer] if a.layer else h.X
        libraries = []
        for start in range(0, h.n_obs, a.block):
            m = matrix[start:start+a.block]
            validate_counts(m, a.cells)
            libraries.append(np.asarray(m[:, idx[ok]].sum(1)).ravel())
        libraries = np.concatenate(libraries)
        positive = libraries[libraries > 0]
        if not len(positive): raise ValueError("no positive raw libraries")
        scale = 1e6

        def normalized_stats(rows):
            total = np.zeros(ok.sum()); squares = np.zeros(ok.sum())
            for start in range(0, len(rows), a.block):
                rr = rows[start:start+a.block]
                m = matrix[rr]
                m = (m.tocsr() if sp.issparse(m) else sp.csr_matrix(m))[:, idx[ok]].astype(float)
                m = sp.diags(scale/np.maximum(libraries[rr], 1)) @ m
                total += np.asarray(m.sum(0)).ravel(); squares += np.asarray(m.multiply(m).sum(0)).ravel()
            avg = total/len(rows); natural = avg
            var = np.maximum(squares/len(rows)-avg**2, 0)/len(rows)
            unit = float(np.mean(scale/np.maximum(libraries[rows], 1)))/len(rows)
            relvar = mean_variance(natural, var, unit)/(natural+pc)**2
            full = np.full(len(genes), np.nan); rv = full.copy()
            full[ok] = natural; rv[ok] = relvar
            return full, rv

        ctrl, ctrl_var = normalized_stats(ci)
        prof = moments(h, ci, idx, a.layer, a.block)
        save_profile(out / "control_profile.npz", genes, prof, dict(kind=a.kind, baseline_convention="ratio_of_sums_cpm"))
        gp = {g: i for i, g in enumerate(genes)}; qc = []; groups = {}
        guides = h.obs[a.guide_col].astype(str).values if a.guide_col else labels
        # A missing guide column means all cells of the target form one group.
        grouped = pd.DataFrame({'target': labels, 'guide': guides}).groupby(['target', 'guide'], sort=True).indices
        for (target, guide), rr in grouped.items():
            if target in control_values: continue
            if len(rr) < a.min_cells: continue
            mean, rv = normalized_stats(rr)
            lfc = np.log2((mean+pc)/(ctrl+pc))
            j = gp.get(target); own = lfc[j] if j is not None else np.nan
            passed = not np.isfinite(own) or own < a.min_ontarget_lfc
            qc.append(dict(target=target, guide=guide, n_cells=len(rr), ontarget_lfc=own,
                           measured=np.isfinite(own), passed=passed))
            if passed and (target not in groups or len(rr) > groups[target][0]):
                groups[target] = (len(rr), lfc, np.sqrt(rv+ctrl_var)/np.log(2))
        targets = sorted(groups)
        if not targets: raise ValueError("no successful target groups after guide QC")
        pd.DataFrame(qc).to_csv(out / "guide_qc.csv", index=False)
        write_effects(out, genes, targets, ((i, groups[t][1], groups[t][2]) for i, t in enumerate(targets)),
                      dict(kind=a.kind, source=str(Path(a.cells).resolve()), label_estimator="per_cell_cpm",
                           normalization_target=scale, label_pseudocount=pc, official_epsilon=PSEUDO,
                           estimator_matches_official=True,
                           baseline_convention="ratio_of_sums_cpm", se_kind="cell_level_delta_per_cell_cpm_poisson_floor",
                           legacy=False), len(targets))
        context_file(out, a.name, a.background, a.platform, a.priority, a.kind)
    finally:
        h.file.close()
    print(f"wrote {len(targets)} successful raw-cell targets to {out}")


def cmd_manifest(a):
    """Explicit raw-cell manifest avoids guessing source labels/layers or downloading huge inputs."""
    records = pd.read_csv(a.manifest).fillna('').to_dict('records')
    for record in records:
        required = ['cells', 'name', 'background', 'obs_col', 'control']
        if any(not record.get(k) for k in required): raise ValueError(f"manifest needs {required}")
        options = dict(cells=record['cells'], genes=a.genes, out=str(Path(a.out)/record['name']),
                       name=record['name'], background=record['background'], obs_col=record['obs_col'],
                       control=record['control'], platform=record.get('platform') or '3prime',
                       kind=record.get('kind') or 'replogle', guide_col=record.get('guide_col') or None,
                       gene_col=record.get('gene_col') or None, layer=record.get('layer') or None,
                       target_split=record.get('target_split') or None, priority=int(record.get('priority') or 0),
                       min_cells=20, min_ontarget_lfc=-0.1, block=512, pseudocount=LABEL_PSEUDOCOUNT)
        cmd_cells(argparse.Namespace(**options))


# ---------------------------------------------------------------- contexts.json
def cmd_config(a):
    ctx = []
    for d in a.dirs:
        d = Path(d); c = json.loads((d / "context.json").read_text())
        ctx.append(dict(name=c["name"], background=c["background"], platform=c["platform"],
                        source_priority=c.get("source_priority", 0), effects=str(d.resolve()),
                        profile=str((d / "control_profile.npz").resolve())))
    bgs = {}
    for c in ctx:
        bgs.setdefault(c["background"], []).append(c["name"])
    if len(bgs) < 2:
        raise SystemExit("need at least two biological backgrounds")
    cfg = dict(genes=str(Path(a.genes).resolve()), effect_cap=a.effect_cap,
               ridge_alphas=[10, 100, 1000, 10000], contexts=ctx, signal_fraction=0.0,
               line_features=True, baseline_mode="bulk", selection_backgrounds=["CD4_T", "HCT116", "HEK293T", "H1"],
               transfer_grid=[dict(temperature=t, gate_center=c, gate_scale=1.0, detection_center=0.05,
                                   detection_scale=0.05) for t in [0.1, 0.5, 2.0] for c in [0.0, 1.8]],
               neighbours=str(Path(a.neighbours).resolve()) if a.neighbours else None, neighbour_k=5)
    save_json(a.out, cfg)
    print(json.dumps(cfg, indent=1)); print("backgrounds (holdout names):", sorted(bgs))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for kind, bg, plat in (("orion", None, "orion"), ("h1", "H1", "flex")):
        s = sub.add_parser(kind); s.add_argument("--pseudobulk", required=True); s.add_argument("--genes", required=True)
        s.add_argument("--out", required=True); s.add_argument("--name", required=True)
        s.add_argument("--background", required=bg is None, default=bg); s.add_argument("--platform", default=plat)
        s.add_argument("--priority", type=int, default=0); s.add_argument("--control")
        s.add_argument("--min-ontarget-lfc", type=float, default=-0.1)
        s.add_argument("--lfc-estimator", choices=["per_cell", "pooled"], default="per_cell")
        s.add_argument("--min-cells", type=int, default=20); s.add_argument("--chunk", type=int, default=256)
        s.add_argument("--pseudocount", type=float, default=LABEL_PSEUDOCOUNT, help="CPM added to both sides of label ratios")
    t = sub.add_parser("tcell"); t.add_argument("--effects", required=True); t.add_argument("--genes", required=True)
    t.add_argument("--out", required=True); t.add_argument("--name", required=True)
    t.add_argument("--background", default="CD4_T"); t.add_argument("--platform", default="flex")
    t.add_argument("--priority", type=int, default=0); t.add_argument("--chunk", type=int, default=256)
    t.add_argument("--min-ontarget-lfc", type=float, default=-0.1)
    raw = sub.add_parser("cells")
    for key in ["cells", "genes", "out", "name", "background", "obs-col", "control"]:
        raw.add_argument("--"+key, required=True)
    for key in ["guide-col", "gene-col", "layer", "target-split"]: raw.add_argument("--"+key)
    raw.add_argument("--kind", default="replogle"); raw.add_argument("--platform", default="3prime")
    raw.add_argument("--priority", type=int, default=0); raw.add_argument("--min-cells", type=int, default=20)
    raw.add_argument("--min-ontarget-lfc", type=float, default=-0.1)
    raw.add_argument("--block", type=int, default=512)
    raw.add_argument("--pseudocount", type=float, default=LABEL_PSEUDOCOUNT, help="CPM added to both sides of label ratios")
    manifest = sub.add_parser("manifest")
    for key in ["manifest", "genes", "out"]: manifest.add_argument("--"+key, required=True)
    m = sub.add_parser("make-config"); m.add_argument("--dirs", nargs="+", required=True)
    m.add_argument("--genes", required=True); m.add_argument("--out", required=True)
    m.add_argument("--neighbours", help="CSV gene_a,gene_b,score from coessentiality/STRING/CORUM")
    m.add_argument("--effect-cap", type=float, default=8.0)
    a = ap.parse_args()
    {"orion": lambda: cmd_sums(a, "orion"), "h1": lambda: cmd_sums(a, "h1"),
     "cells": lambda: cmd_cells(a), "manifest": lambda: cmd_manifest(a),
     "tcell": lambda: cmd_tcell(a), "make-config": lambda: cmd_config(a)}[a.cmd]()
