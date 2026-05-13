# Recap: `us-all-results-extend.R`

## What this script does

`us-all-results-extend.R` is a **thin wrapper** that delegates to the canonical implementation in `refactor_results/12_run_extend_results.R`.

At a high level it:

1. Sets the working directory to the script location.
2. Sources `refactor_results/12_run_extend_results.R`, which handles all scoring and output.

The extend lane covers **Morris baseline + AR2 + MRTS** settings from `settings.csv`, all scored against result files in `code/analysis/ozone/US-all/results/`.

---

## What `12_run_extend_results.R` does

1. Loads shared experiment objects from `us-all-setup.RData` and scoring helpers from `../../../R/auxfunctions.R`.
2. Reads `settings.csv` to determine which settings to score (IDs 1-209).
3. Loads `results/us-all-<setting>.RData` for each available setting.
4. Scores each setting using `QuantScore` and `BrierScore` (quantile scores, Brier scores, CRPS, coverage, PIT calibration).
5. Produces comparison outputs with the `_extend` suffix:
   - `output/us-all/results/us-all-results-extend.RData`
   - `output/us-all/tables/comparison_full_table_extend.csv`
   - `output/us-all/tables/comparison_top2_extend.csv`
   - `output/us-all/tables/comparison_paired_same_basis_extend.csv`
   - `output/us-all/tables/comparison_scalar_metrics_extend.csv`
   - `output/us-all/tables/comparison_top2_extend.xlsx`

---

## Why it is comparable (same 基準)

The extend lane is comparable to the Gaussian baseline because it preserves the exact evaluation basis:

- **Same data setup**: uses the same `us-all-setup.RData`.
- **Same CV splits**: uses the same `cv.lst` folds.
- **Same scoring definitions**: calls `QuantScore(...)` and `BrierScore(...)` from the same helper source.
- **Same quantile grid / thresholds**: identical `probs` and thresholds derived the same way.
- **Same reference model**: relative scores normalized by Gaussian setting `1` (`us-all-1.RData`).
- **No in-place overwrite of baseline pipeline**: baseline script `us-all-results.R` remains frozen for reproducibility.

---

## Naming convention

The `_extend` suffix distinguishes this lane from the other two:

| suffix     | source                                              | baseline                  |
| ---------- | --------------------------------------------------- | ------------------------- |
| `_extend`  | `code/analysis/ozone/US-all/results/`               | setting 1 (this lane)     |
| `_proposed`| `code/analysis/ozone_prop/US-all/results/`          | reuses extend's setting 1 |
| `_pool`    | merged extend ∪ proposed (no re-scoring)            | extend's setting 1        |

See `output/us-all/README.md` for the full suffix-to-source table and pool ID layout.

---

## How to run

```ps
cd D:\Github\spatial-skew-t\code\analysis\ozone\US-all
$env:US_ALL_SUMMARY_DRAWS='5000'   # optional; defaults to 400 if omitted
Rscript us-all-results-extend.R
```

Or call the refactored runner directly:

```ps
Rscript refactor_results/12_run_extend_results.R
```

---

## Notes

- `settings$setting` can include non-numeric IDs (like `5a`), so the runner uses `as.integer(setting)` where needed and skips non-integer rows safely.
- Setting `2` has special prediction behavior in the legacy workflow; the runner includes a contract/scoring guard path to avoid aborting the whole run.
- Classification diagnostics are **not** computed in this lane (those were added later in the proposed lane). Pool files mark extend rows with `NA` for classification columns.
- The historical name for this runner was `us-all-results-proposed.R` / `12_run_proposed_results.R` — renamed to `_extend` to disambiguate after the `ozone_prop` lane introduced a new meaning for "proposed".
