# Cell QC and perturbation QC
Every single-cell source passes the same three cell-QC steps, then perturbation QC. The rules were fixed before any modelling and applied identically to every source. They are recorded per source in `results/qc/<source>/qc_method.md` and `checks.json`.

## 0. Raw files (`scripts/build_raw_cells.py`, `scripts/convert_new_contexts.py`)
- **One file per cell line** (`lib/cells/<NAME>.h5ad`): integer UMI counts (CSR int32), the source's own genes, and a standard `obs`:
  - `cell_line`, `dataset`, `condition`;
  - `target_gene` (`non-targeting` for controls);
  - `n_targets`, `perturbation`, `guide`, `batch` (lane);
  - `total_umi`, `umi_source` (`library` = the source's per-cell total over all genes; `matrix` = the sum over stored genes), `n_genes_detected`.
- **No cell is dropped at this step.**
- **New contexts** (KOLF2.1J, Adamson, Feng, Jiang): CRISPRi or knockout single-target cells only.
  - Feng's genome-wide screens have too few controls (12–36 per line). There, the **rest-of-line cells** (all assigned cells of other knockdowns in the same donor line) stand in as the reference.

## 1. UMI cut-offs (`scripts/qc_replogle_umi.py`)
- On x = log10(UMI + 1) of **control cells**:
  - low = median − 3 MAD, high = median + 4 MAD;
  - MAD is scaled by 1.4826 and floored at 0.05;
  - an absolute floor of 200 UMIs.
- **Per lane** (gem group / batch). A lane with too few controls falls back to the whole-file controls. The lane minimum is 200 controls for K562, 100 for RPE1 and **50 for every later source**.
- **Pre-filtered sources:** if the lowest histogram bin holds ≥ 1% of the peak bin, the authors already filtered the file (a sharp left edge, as in HepG2, Jurkat and H1). The low cut is then skipped.
- **Knockdowns that lose many cells:** with ≥ 30 cells, a knockdown whose share removed by the per-lane cut exceeds max(2 × the control rate, control rate + 5 points) is flagged. Its low cut is then taken from its **own** cells, so that strong knockdowns that shrink the library are not removed for it.

## 2. Apply (`scripts/qc_replogle_apply.py`)
1. Every cell: removed below 200 UMIs or above its lane's high cut.
2. Controls and unflagged knockdowns: removed below the lane's low cut.
3. Flagged knockdowns: removed below their own low cut.

Outputs barcode lists only (`kept_cells.txt`, `removed_cells.csv`).

## 3. Write (`scripts/write_replogle_qc.py`)
- **Source checks:**
  - every count is a whole number;
  - where a UMI column exists, it is the raw depth: ≥ 99.9% of cells have row sum ≤ UMI and ≥ 0.9 × UMI. Released matrices omit lowly expressed genes.
- **Output** `lib/qc_cells/<source>_qc.h5ad`: kept cells, counts unchanged, standardised `obs` (`source`, `background`, `batch`, `target_gene`, `is_control`, `target_gene_id`, `guide_id`, `n_umi`, `qc_low_ref`). The full QC record is in `uns['qc']`.
- **After writing:**
  - cell count = barcode list;
  - 100 random cells identical to the source;
  - source file size and modification time unchanged.
- **H1:** its QC'd file is used for training only. The scoring harness keeps the unfiltered H1 cells.

## 4. Perturbation QC (`scripts/qc_perturbations.py`)
Parameters: MIN_CELLS 30, MIN_CTRL_CPM 5, MAX_ONTARGET_LFC −0.5, P_MAX 0.05. Per-cell CPM = counts / n_umi × 10⁶.
- **Per construct (guide):**
  - measurable = the target is in the panel and its control mean CPM ≥ 5;
  - on-target LFC = log2((construct mean CPM + 0.5) / (control mean CPM + 0.5));
  - p = one-sided Mann–Whitney U (construct < control);
  - **keep = not measurable, or (LFC ≤ −0.5 and p < 0.05).** A knockdown that can't be verified is kept, not discarded.
- **Per knockdown:** cells are pooled from kept constructs. The knockdown is kept if n_cells ≥ 30. *Efficiency* = the pooled on-target LFC.
- **Sanity check:** 50 random constructs from each of three groups (dropped; −1 < LFC ≤ −0.5; LFC ≤ −1), with their count of downstream genes moving |log2FC| > 0.5.
- **Table-only sources** (Orion pseudobulk, Marson CD4 T): `scripts/qc_table_perturbations.py` applies the same on-target rule. Their own two-sided test is converted to a one-sided p.

## Fixed decisions (apply to every source)
- On-target threshold LFC ≤ −0.5 for every source; MIN_CELLS = 30; lane minimum 50 controls.
- SE: delta method, floor 0.05, no rescaling. **Missing genes are NaN, never 0.**
- Pre-filtered sources skip the low UMI cut. UMI checks are file-level.
- **Dropped datasets:** Tian (iPSC, neurons), MCF10A and LUHMES (weak knockdowns). HFF kept as 0 h cells only.
- **Same cell line never votes for itself:** at prediction time, every table whose `obs.background` matches the target line is excluded, not only same-named files.
