# US-all compare tables (all settings: baseline + AR2 + MRTS)

This guide explains how to regenerate:

- `comparison_full_table.csv`
- `comparison_top2.csv`

with **all numeric settings** from `settings.csv` compared together, including:

- Morris baseline lane
- AR2 lane (`ar2 == yes`)
- MRTS lane (`mrts > 0`)

---

## 1) What script to run

Use this entrypoint:

- `us-all-results-proposed.R`

It now delegates to:

- `refactor_results/12_run_proposed_results.R`

The refactored script computes scores and writes output tables.

---

## 2) How settings are included

The script reads `settings.csv` and uses:

- `all_numeric_settings`: all rows with numeric `setting` IDs (e.g., `1..124`, `201..209`)
- `done_ar2`: settings where `ar2 == yes`
- `done_mrts`: settings where `mrts` is a positive integer
- `done_morris`: fixed historical baseline lane used as baseline/provenance flag

So AR2 and MRTS are both included in one unified comparison pass.

---

## 3) How result files are discovered

For each setting `s`, the script searches for:

- `us-all-<s>.RData`

across result directories in priority order:

1. directories in `US_ALL_RESULTS_DIRS` (comma-separated, if set)
2. `US_ALL_RESULTS_DIR` (if set)
3. `results`
4. `results_new`

This allows baseline/AR2 files in `results` and MRTS files in `results_new` to be combined automatically.

---

## 4) Run commands

From `code/analysis/ozone/US-all`:

```powershell
# (optional) explicit search order
$env:US_ALL_RESULTS_DIRS = "results,results_new"

# run comparison build
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" "us-all-results-proposed.R"
```

If you only use one directory:

```powershell
$env:US_ALL_RESULTS_DIR = "results"
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" "us-all-results-proposed.R"
```

---

## 5) Outputs and meaning

The script writes:

- `comparison_full_table.csv`
  - Long table by setting × metric × quantile
  - Includes metadata/flags:
    - `is_baseline`
    - `is_proposed`
    - `is_ar2`
    - `is_mrts`
    - `model_lane`
    - `mrts`

- `comparison_top2.csv`
  - Top-2 settings for Brier relative score at target quantiles (`0.95, 0.98, 0.99, 0.995`)
  - Also includes `is_ar2`, `is_mrts`, `model_lane`, and `mrts`

- `comparison_paired_same_basis.csv`
  - Paired baseline-vs-extension comparison (AR2 and MRTS settings both included as proposed lanes)

- `us-all-results-proposed.RData`
  - Saved objects for reproducible downstream analysis

---

## 6) Common checks

After running, confirm:

1. Console summary shows searched directories and missing settings.
2. `comparison_full_table.csv` contains rows with settings in AR2 (`101..124`) and MRTS (`201..209`) when those files exist.
3. `comparison_top2.csv` includes `model_lane` and `is_mrts`/`is_ar2` flags.

---

## 7) Notes

- Gaussian setting `1` must be available; otherwise relative scores cannot be computed.
- Missing files are skipped and reported (the pipeline is fault-tolerant).
- Non-numeric setting IDs (e.g., `5a`) are not part of the all-numeric namespace used for these two comparison tables.
