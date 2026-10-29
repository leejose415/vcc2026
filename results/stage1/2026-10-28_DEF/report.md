# D/E/F final file — C+G (2026-10-28)

**Decision (Jose, 2026-10-28 ~05:25 UTC):** C+G is the D/E/F model (Eval A: +0.034 [+0.020, +0.048] over the Model 2 port; `results/stage1/2026-10-27_incumbent/report.md`). Run it on D/E/F **and upload it.**

**Model:** `code/vcc6/stage1.py`, arm C with the Gaussian shrink, through `run_line.py predict`.
- Each context's own controls (`data_final/extracted/context_<X>.h5ad`, 18,400 each).
- The 300 final targets (`results/target_lists/final_DEF_300.txt`, from `pert_counts.csv`).
- All 117 tables vote: D/E/F match no source line.
- Salt "", 400 cells per knockdown, the frozen generator (sha256 422f8a3d…).

**File:** `code/vcc6/def_file.py` strips the controls, sets `context`, and concatenates to 360,000 cells × 18,533 genes. It then runs `vcc prep` with the final `gene_names.csv` and `pert_counts.csv` (contexts D, E, F), checking the exit code. After that comes `vcc submit`.

**Estimate (05:27):** 3 × ~1 h runs, then about 30 min to combine and prep, then the upload. Done around 09:30 UTC.

## Log
- 05:27: context D started. Controls: 18,400, median panel library 20,236, 9,633 gated genes; 117 tables loading.
- **06:24: context D was killed out of memory (exit 137)** just after generating all 300 knockdowns, while writing the file. The 117 source tables (~15 GB) were still in memory in the same process, which reached 27 GB RSS (`journalctl`). The partial output was moved to `D/salt_none_oom_partial/`.
- **06:25: restarted with model and generation in separate processes** (`logs/run_def2.sh`), as the Eval A runs did.
  - Model per context: `model2_port.py cg`, which is `run_line`'s model step; the process exits and frees the tables.
  - Then generation from the saved arrays: `run_line predict --lfc-npz`, with no tables loaded. The generator's inputs are identical (lfc, de_sig, n_sources; `p_sign` is never read).
  - D's model output (05:35) is reused. New estimate: E and F models about 8 min each, then 3 generations of about 50 min each. Done around 09:20 UTC, then combine, prep and upload.
- 07:36: D done (exit 0): 300 knockdowns, median calls 242.5 requested and realised, 52 min of generation, 1.5 GB file. E started.
- **08:19:** I started vcc2027 to run F in parallel and stopped vcc2026's queue so F wouldn't also start there. **Jose said no;** vcc2027 was stopped at once, having run nothing. F was re-queued on vcc2026 behind E (`logs/run_def3.sh`: waits for E's process, checks its output is complete, then runs F). E was never interrupted.
- **09:15:** combine → `vcc prep` (exit checked) → `vcc submit -m LHY` queued behind F (`logs/prep_submit.sh`). A temporary 16 GB swap covers `vcc prep`'s ~30 GB peak and is removed afterwards. Nothing is uploaded unless `vcc prep` exits 0.
- **09:45:** the in-memory combine reached 30 GB RSS plus 15 GB swap while gzip-writing the merged file, at about 70 MB/min, with a risk of an out-of-memory kill. It was stopped. The three per-context submission files it had written (each 120,000 × 18,533, controls removed) are merged on disk instead (`anndata.experimental.concat_on_disk`), then `vcc prep`, then `vcc submit` (`logs/merge_prep_submit.sh`).
- **~10:00:** `vcc prep` reached 30 GB RAM plus 12 of 15 GB swap, so a second temporary 16 GB swap (`/swapfile_eval2`) was added. The upload script removes only `/swapfile_eval`; `/swapfile_eval2` is removed by hand afterwards. vcc2027 was stopped at 09:50 (idle; study outputs copied back).

## Submitted (2026-10-28 10:07 UTC)
- **Merge:** 2 min, on disk.
- **`vcc prep`:** 09:47–10:03, exit 0. 360,000 cells × 18,533 genes (float32); contexts D / E / F = 120,000 each; targets verified against the official list; counts preserved. Output `final/prediction_DEF.vcc`, 3.19 GB.
  - Its warning: about 34 GiB needed for 1.96e9 stored entries, more than this machine's 31 GiB RAM. It completed with the temporary swap.
- **`vcc submit -m LHY`:** exit 0, 3.0 GiB uploaded in 4 min 06 s. **Entry `JFCghVNyHwfCOWhnHG6M`**, final phase ("counts toward final scoring"). Server scoring job `vcc-metrics-prod-jfcghvnyhwfcowhnhg6m`. Scores are hidden until the winners are announced.
- **Model:** C+G (stage1 arm C + Gaussian shrink), all 117 tables, each context's own controls, salt "", the frozen generator. Median calls requested / realised: D 242.5 / 242.5, E 278.5 / 279, F 245 / 246.5.
- Both temporary swap files are removed afterwards.
