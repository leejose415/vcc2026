"""Profile an unseen line using ALL available raw-count control cells.

Transfer context uses ratio-of-sums CPM (bulk_cpm), not mean per-cell CPM.
Both are retained: mean_cpm/var_cpm/se_mean keep their original per-cell meaning
for DE diagnostics and the existing cell generator. This does not convert LFCs.

Example: python vcc_context_profile_v3.py --counts context_A.h5ad --all-controls
         --genes gene_names.csv --out profiles/A.npz
"""
import argparse
from pathlib import Path

import anndata as ad
import numpy as np

from vcc_common import alignment, moments, program_levels, read_names, save_profile, select_controls


def profile(counts, genes_csv, output, layer=None, obs_col="target_gene",
            control="non-targeting", all_controls=False, block=512):
    genes = read_names(genes_csv, "gene_name")
    a = ad.read_h5ad(counts, backed="r")
    try:
        idx = alignment(a.var_names.astype(str), genes)
        rows = select_controls(a, obs_col, control, all_controls)
        stats = moments(a, rows, idx, layer, block)
        # Explicit transfer representation. Do not overwrite mean_cpm/log_mean:
        # the DE scorer and the cell generator still need per-cell mean CPM.
        stats["baseline_cpm"] = stats["pooled_cpm"].copy()
        stats["baseline_log"] = np.log1p(stats["pooled_cpm"])
        stats["baseline_convention"] = np.array("ratio_of_sums_cpm")
        stats["program_names"], stats["baseline_program_levels"] = program_levels(stats["pooled_cpm"], genes)
        finite = np.isfinite(stats["baseline_log"])
        centre = float(stats["baseline_log"][finite].mean())
        stats["baseline_log_centered"] = stats["baseline_log"] - centre
        # These identifiers permit auditing overlap when making synthetic predictions.
        stats["control_obs_names"] = np.asarray(a.obs_names[rows], dtype=str)
        save_profile(output, genes, stats, dict(counts=str(Path(counts).resolve()), layer=layer,
                     obs_col=obs_col, control=control, all_controls=all_controls,
                     transfer_baseline_convention="ratio_of_sums_cpm",
                     transfer_baseline_field="pooled_cpm",
                     per_cell_mean_field="mean_cpm",
                     baseline_log_center=centre,
                     library_denominator="sum over the provided gene_names.csv axis"))
    finally:
        a.file.close()
    print(f"Saved {output}: {len(rows):,} controls, {len(genes):,} aligned genes")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--counts", required=True)
    p.add_argument("--genes", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--layer")
    p.add_argument("--obs-col", default="target_gene")
    p.add_argument("--control", default="non-targeting")
    p.add_argument("--all-controls", action="store_true")
    p.add_argument("--block", type=int, default=512)
    a = p.parse_args()
    profile(a.counts, a.genes, a.out, a.layer, a.obs_col, a.control, a.all_controls, a.block)
