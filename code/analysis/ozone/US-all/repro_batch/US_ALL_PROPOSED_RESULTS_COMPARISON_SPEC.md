# US-all Proposed Results Comparison Spec

This document is the implementation-aligned spec for:

- `code/analysis/ozone/US-all/us-all-results-proposed.R`

Its purpose is to define how proposed models are compared against Morris baseline under the same benchmark basis.

---

## 1) Goals

1. Keep Morris baseline workflow reproducible and untouched.
2. Compare proposed settings with baseline under the same data/CV/scoring basis.
3. Generate reproducible comparison tables for reporting.

---

## 2) Canonical files

### Required inputs

- `us-all-setup.RData`
- `settings.csv`
- `results/us-all-<setting>.RData`
- `../../../R/auxfunctions.R`

### Aggregation script

- `us-all-results-proposed.R`

### Outputs

- `us-all-results-proposed.RData`
- `comparison_full_table.csv`
- `comparison_top2.csv`
- `comparison_paired_same_basis.csv`

---

## 3) Setting scope and namespaces

### Baseline (fixed)

`done_morris <- c(1:5, 7:9, 11:13, 15:17, 33:36, 38:41, 43:46, 51:74)`

### Proposed (dynamic from settings metadata)

Proposed settings are discovered from:

- `settings$ar2 == "yes"`

Then coerced via numeric setting id (`setting_num`) and filtered to non-NA integers.

### Final requested set

- `all_requested <- sort(unique(c(done_morris, done_proposed)))`

---

## 4) Same-basis comparability contract

Comparability is ensured because both baseline and proposed settings are scored with:

1. Same setup object (`us-all-setup.RData`)
2. Same CV folds (`cv.lst`)
3. Same quantile grid (`probs`)
4. Same thresholds (`quantile(Y, probs=...)`)
5. Same score functions (`QuantScore`, `BrierScore`)
6. Same Gaussian reference (setting `1`) for relative metrics

This is the required basis for fair model ranking comparison.

---

## 5) Result object contract for each setting

Each `results/us-all-<setting>.RData` must provide:

- object `fit`
- `fit` as list with one entry per fold (`length(cv.lst)` expected)
- fold-level posterior predictions:
  - `fit[[d]]$yp`

Expected scoring-compatible dimensions:

- `pred <- fit[[d]]$yp` with 3 dimensions
- interpreted as: `[iter, n_validation_sites, n_time]`
- must align with `validate <- Y[val.idx, ]`

---

## 6) Robustness behavior in `us-all-results-proposed.R`

The script is intentionally fault-tolerant:

- Missing result file -> tracked in `skipped_missing_file`
- Missing `fit` / bad fold contract -> tracked in `skipped_bad_contract`
- Scoring error (`QuantScore`/`BrierScore`) -> tracked in `skipped_scoring_error`

Processing continues for other settings; invalid settings are skipped.

Hard stop condition:

- If setting `1` (Gaussian reference) is unavailable, script stops.

---

## 7) Computed outputs and their meaning

### `comparison_full_table.csv`

Long-form table over:

- setting × quantile × metric (`brier`, `quantile`)

Includes:

- absolute means and SE
- relative-to-Gaussian value
- model metadata (`method`, `knots`, `thresh`, `CMAQ`, `TS`, `ar2`)
- baseline/proposed flags

### `comparison_top2.csv`

Top-2 rows for Brier relative score at target quantiles:

- `0.95`, `0.98`, `0.99`, `0.995`

### `comparison_paired_same_basis.csv`

Paired baseline-vs-proposed table by matched configuration:

- match keys: `method`, `knots`, `thresh`, `CMAQ`, `TS`

Includes delta metrics:

- `delta_rel_brier`
- `delta_rel_quant`

---

## 8) Differences from legacy integration assumptions

This spec supersedes earlier assumptions that are no longer accurate:

1. No separate `settings_proposed.csv` is required.
   - Current implementation reads only `settings.csv` and uses `ar2` marker.
2. Proposed settings are not hard-coded in the script.
   - They are inferred dynamically from metadata.
3. Output set includes paired same-basis table.
   - `comparison_paired_same_basis.csv` is part of standard deliverables.
4. Fault handling is explicit and recorded.
   - Skipped settings are tracked in summary vectors and saved in `.RData`.

---

## 9) Recommended workflow

1. Use `repro_batch/batch_auto_check.R` to ensure result files exist.
2. Run `us-all-results-proposed.R` from `code/analysis/ozone/US-all`.
3. Inspect summary lines for skipped settings.
4. Use the three CSV outputs for ranking and reporting.
5. Treat paired table as primary evidence for same-basis baseline vs proposed comparison.

---

## 10) Naming note

This file intentionally uses a descriptive name:

- `US_ALL_PROPOSED_RESULTS_COMPARISON_SPEC.md`

to emphasize that it is a **results-comparison spec aligned with current implementation**, not a generic future integration draft.