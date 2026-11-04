# Data sources
None of the data is in this repository. Each source below is public. Raw files go under `~/vcc/lib/`, and the challenge's own data under `~/vcc/data/` (validation lines A–C) and `~/vcc/data_final/` (final lines D–F).

## Challenge data (Arc Institute)

| What | How to get it |
|---|---|
| H1 hESC training screen (Arc 2025): about 300 knockdowns, used for Eval A | Arc's challenge release (`vcc datasets`, or the challenge website). Our copy was `lib/h1/H1_full_2025.h5ad` |
| Control cells and target lists: A–C (validation) and D–F (final, 18,400 controls each, 300 targets) | `vcc datasets download …` with the `vcc` CLI ≥ 0.3.0; `data_final/extracted/{context_D,E,F}.h5ad`, `pert_counts.csv`, `gene_names.csv` |

## Perturbation screens (training sources)

| Source | Cell line(s) | Accession / URL | Our layer |
|---|---|---|---|
| Replogle 2022, genome-wide | K562 | figshare / plus.figshare article **20029387** (`K562_gwps_raw_singlecell_01.h5ad`) | `replogle_k562_gw` |
| Replogle 2022, essential | RPE1 | figshare article **20029387** (file 35775606) | `replogle_rpe1` |
| Nadig 2025 | HepG2, Jurkat | GEO **GSE264667** (`*_raw_singlecell_01.h5ad`) | `nadig_hepg2`, `nadig_jurkat` |
| X-Atlas Orion (Xaira) | HCT116, HEK293T | Hugging Face `Xaira-Therapeutics/X-Atlas-Orion` (`scripts/orion_stream.py` streams only the needed targets) | `orion_hct116`, `orion_hek293t` (+ pseudobulk tables) |
| Marson lab, CD4⁺ T cells | primary CD4 T (rest, stim 8 h, stim 48 h) | public replicate-level tables. *Accession to fill in* | `marson_cd4_*` |
| Feng (HipSci iPSC CRISPRi) | 19 donor lines (targeted) + genome-wide fitness / non-fitness screens | public count CSVs. *Accession to fill in* | `feng_targeted_*`, `feng_gw_*` |
| KOLF2.1J CRISPRi | KOLF2.1J hiPSC | public `KOLF_Strong_Perturbations.h5ad`. *Accession to fill in* | `kolf_kolf21j` |
| Adamson 2016 | K562 (UPR screens 10X001 / 005 / 010) | public (GEO). *Accession to fill in* | `adamson_k562_*` |
| Jiang (cytokine-response CRISPR screens) | A549, BxPC3, HAP1, HT29, K562, MCF7 × IFNβ, IFNγ, insulin, TGFβ, TNFα | public (GEO). *Accession to fill in*; `scripts/extract_jiang_counts.R` | `jiang_*` (knockout) |
| Tahoe-100M (Tahoe Therapeutics) | 50 cancer lines × 379 drugs | Hugging Face `tahoebio/Tahoe-100M` (3,388 parquet shards; `scripts/convert_new_contexts.py`) | `tahoe_*` (drug modality, low weight) |
| PerturBase | hESC (Genga 2019), HFF (Hein 2021) | perturbase.cn raw releases (PRJNA524331 and others) | explored, not in the final 117 |

## Gene-level resources

| Resource | Use | Source |
|---|---|---|
| HGNC complete set (2026-10-28) | symbol mapping (previous / alias) | genenames.org |
| DepMap 24Q4 `CRISPRGeneEffect.csv` | gene features (F5 term, Q2) | figshare article 27993248 |
| STRING v12.0 (human links + info) | gene features (tested) | stringdb-downloads.org |
| MSigDB Hallmark 2020 | pathway programs (Step 22) | Enrichr library via `gseapy.get_library` |
| Seurat `cc.genes` (Tirosh 2016) | cell-cycle phases (Step 22) | `satijalab/seurat` `data/cc.genes.rda` |
| Cellosaurus | TP53 status of source lines (Step 21) | api.cellosaurus.org |
| ESM2 / STRING / DepMap / text gene embeddings | transformer gene features | exported by `scripts/export_pie_features.py` (not in this release) |

## Notes
- **Two source symbols on one Ensembl ID:** keep the panel's symbol and drop the alias; never sum.
- **Blank Ensembl IDs** are filled from the panel by symbol unless the ID is already taken.
- Some scripts fall back to our (now deleted) `gs://vcc2026` bucket for raw files; use the URLs above instead.
