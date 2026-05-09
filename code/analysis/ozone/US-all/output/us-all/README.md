# US-all outputs

All generated outputs for the US-all ozone workflow are collected here.

## Folder layout

- `results/` — serialized R objects and score caches (`.RData`)
- `tables/` — CSV summaries and comparison tables (+ Excel workbooks)
- `plots/` — figures and map exports
- `logs/` — run logs or batch output

## File-name suffix convention

The suffix on a filename tells you **which subset of model fits** the comparison covers and **which directory of `fit$yp` files** it was scored from.

| suffix on file | data scope                                                                                      | scored from                                                        | Gaussian baseline (setting 1)                              | written by                                                |
| -------------- | ----------------------------------------------------------------------------------------------- | ------------------------------------------------------------------ | ---------------------------------------------------------- | --------------------------------------------------------- |
| `-a`           | CMAQ-vs-no-CMAQ paired lane                                                                     | `code/analysis/ozone/US-all/results/`                              | setting 1                                                  | `11_run_a_results.R`                                      |
| `_extend`      | Morris baseline + AR2 + MRTS lane (the "extended" setting space from `settings.csv`, IDs 1-209) | `code/analysis/ozone/US-all/results/us-all-{i}.RData`              | from same lane (`us-all-1.RData`)                          | `12_run_extend_results.R`                               |
| `_proposed`    | 10 skew-t fits from the `ozone_prop` backend (orig prop IDs 1-10 → unified IDs 2-11)            | `code/analysis/ozone_prop/US-all/results/ozone-prop-{i}.RData`     | reuses `code/analysis/ozone/US-all/results/us-all-1.RData` | `14_run_prop_results.R`                                   |
| `_pool`        | Pooled comparison: extend ∪ proposed in a single namespace (extend IDs 1-209, prop IDs 301-310) | merged from existing `_extend` + `_proposed` RData (no re-scoring) | setting 1 (extend's Gaussian)                              | `15_pool_extend_prop.R`                                   |
| `_mrts_cov`    | MRTS-covariance focused subset (separate dedicated lane)                                        | `code/analysis/ozone/US-all/results/`                              | varies                                                     | `us-all-results-mrts-cov.R`                               |

The same suffixes apply to both `results/*.RData` and `tables/comparison_*.csv` / `.xlsx`.

> **Note on the historical no-suffix files.** Before this convention, the runner was named `12_run_proposed_results.R` (now `12_run_extend_results.R`) and wrote `us-all-results-proposed.RData` plus `comparison_*.csv` (no suffix), which conflicted with the new `_proposed` lane (the ozone_prop fits). The runner and its plot companion (now `13_plot_extend_results.R`) were renamed and patched to use the `_extend` suffix, and the existing tables on disk were renamed accordingly. Any leftover historical files without a suffix (e.g. `comparison_top2.updated.xlsx`) are out-of-band manual exports.

## Pool ID layout (`_pool` files)

| pool ID range | source                                                                            | n          |
| ------------- | --------------------------------------------------------------------------------- | ---------- |
| 1             | Gaussian baseline (extend's `us-all-1.RData`)                                     | 1          |
| 2-209         | extend lane: Morris baseline + AR2 + MRTS — original `settings.csv` IDs preserved | 106 scored |
| 301-310       | ozone_prop fits — unified ID = 300 + `setting_orig` from `settings_prop.csv`      | 10         |

`baseline_setting = 1` for all `*.mean.ref.gau` columns. `done_morris`, `done_ar2`, `done_mrts`, `done_prop` are saved as separate vectors so each lane can be filtered. Every comparison table carries `setting_orig`, `source`, `is_baseline`, `is_proposed`, `is_ar2`, `is_mrts`, `is_prop`, and `model_lane` columns added by `decorate_setting_table`.

The `_pool` merge does **not** re-score — it concatenates the existing extend and proposed score arrays along the setting axis, keeping extend's setting 1 as the single Gaussian denominator. Brier scores at the baseline are byte-identical between the two lanes; CRPS / coverage / PIT drift by < 0.3 % due to MCMC draw-subsampling RNG.

Classification metrics in `_pool` files are populated **only for the ozone_prop rows (301-310)** — extend was scored before classification diagnostics were added.

## What goes where

- Use the refactored runners in `refactor_results/`. They write into this folder automatically.
- Keep legacy root-level outputs only for historical reference.
