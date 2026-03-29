# Recap: `us-all-results-proposed.R`

## What this script does

`us-all-results-proposed.R` is a **comparison aggregator** that keeps the original Morris baseline intact and adds the proposed (AR2) models as a parallel lane.

At a high level it:

1. Loads shared experiment objects from `us-all-setup.RData` and shared scoring functions from `../../../R/auxfunctions.R`.
2. Defines fixed baseline settings (`done_morris`) exactly as in the original reproduction.
3. Detects proposed settings from `settings.csv` where `ar2 == "yes"`.
4. Loads `results/us-all-<setting>.RData` for baseline + proposed settings.
5. Scores each setting using the same `QuantScore` and `BrierScore` framework.
6. Produces comparison outputs:
   - `us-all-results-proposed.RData`
   - `comparison_full_table.csv`
   - `comparison_top2.csv`
   - `comparison_paired_same_basis.csv`

---

## Why it is comparable (same 基準)

The script is comparable to the baseline because it preserves the exact evaluation basis:

- **Same data setup**: uses the same `us-all-setup.RData`.
- **Same CV splits**: uses the same `cv.lst` folds.
- **Same scoring definitions**: still calls `QuantScore(...)` and `BrierScore(...)` from the same helper source.
- **Same quantile grid / thresholds**: uses identical `probs` and thresholds derived the same way.
- **Same reference model**: relative scores are still normalized by Gaussian setting `1`.
- **No in-place overwrite of baseline pipeline**: baseline script `us-all-results.R` remains frozen for reproducibility.

So comparisons are apples-to-apples, not apples-to-different-fruit-salad.

---

## Main differences vs `us-all-results.R`

## 1) Scope of settings

- `us-all-results.R`:
  - hard-coded for settings `1..74`
  - evaluates only a fixed baseline `done` subset
- `us-all-results-proposed.R`:
  - keeps baseline `done_morris`
  - adds proposed settings dynamically from `settings$ar2 == "yes"` (e.g., `101+`)

## 2) Array sizing

- `us-all-results.R`: fixed arrays with third dimension `74`
- `us-all-results-proposed.R`: arrays sized by `max(all_requested)` to hold both baseline and proposed IDs

## 3) Robustness / fault tolerance

- `us-all-results.R`: assumes expected files/contracts; failures can stop run
- `us-all-results-proposed.R`:
  - tracks `skipped_missing_file`
  - tracks `skipped_bad_contract`
  - tracks `skipped_scoring_error`
  - continues processing other settings when one fails

## 4) Output products

- `us-all-results.R`:
  - primary output is `us-all-results.RData`
  - contains broader legacy plotting/analysis blocks
- `us-all-results-proposed.R`:
  - focused comparison outputs for baseline vs proposed:
    - full table by setting × quantile × metric
    - top-2 rankings at target quantiles
    - **paired same-basis table** matched by `(method, knots, thresh, CMAQ, TS)`

## 5) Comparison design

- `us-all-results.R`: baseline-only summary and visualization
- `us-all-results-proposed.R`: explicit two-track comparison architecture (baseline preserved, proposed added)

---

## Practical interpretation

If your question is:

> “Can I claim proposed model performance is compared under the same benchmark?”

For this script design, the answer is **yes**—because all scoring ingredients are held constant and only model results (settings lane) differ.

---

## Notes

- `settings$setting` can include non-numeric IDs (like `5a`), so the proposed script safely uses `setting_num <- as.integer(setting)` where needed.
- Setting `2` can have special prediction behavior in legacy workflow; the proposed script includes a contract/scoring guard path to avoid aborting the whole run.
