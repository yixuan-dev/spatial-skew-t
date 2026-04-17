# US-all results refactor

This folder splits the old monolithic result scripts into reusable modules and purpose-specific runners.

## Why this split

The original scripts mixed several responsibilities in one file:

- loading/setup/bootstrap
- score computation
- score summarization and relative metrics
- comparison table generation
- exploratory plotting/diagnostics

This refactor keeps the scoring/comparison path modular and reproducible.

## Module files

- `00_bootstrap.R`
  - script working-directory setup
  - shared setup/data load (`us-all-setup.RData`, `settings.csv`, `auxfunctions.R`)
  - default quantile-score grid: `0.00, 0.10, ..., 0.90, 0.95, 0.98, 0.99, 0.995`
  - default Brier threshold grid: `0.00, 0.10, ..., 0.90, 0.95, 0.98, 0.99, 0.995`
- `01_score_engine.R`
  - result-file scoring loop
  - prediction contract checks
  - summary/relative-score computation
  - optional CRPS / coverage / PIT calibration diagnostics for Excel-ready exports
- `02_comparison_tables.R`
  - full comparison table
  - top-2 ranking table
  - paired same-basis table
  - scalar metrics table
  - uncertainty summary table
  - calibration-bin table

## Runner files

- `10_run_morris_results.R`
  - Morris baseline scoring lane (replacement for scoring part of `us-all-results.R`)
- `11_run_a_results.R`
  - CMAQ vs no-CMAQ lane from `us-all-results-a.R`
- `12_run_proposed_results.R`
  - baseline + AR2 proposed lane from `us-all-results-proposed.R`

## Output compatibility notes

All runner outputs now live together under `output/us-all/` with simple subfolders:

- `results/` for `.RData` outputs
- `tables/` for comparison CSVs and summaries
- `plots/` for figures
- `logs/` for logs

The runners keep the same key filenames used in your current workflow, but write them inside that shared output tree:

- Morris runner: `output/us-all/results/us-all-results-0401.RData`, `output/us-all/results/us-all-results.RData`
- A runner: `output/us-all/results/us-all-results-combined.RData`, `output/us-all/results/us-all-results-a.RData`
- Proposed runner: `output/us-all/results/us-all-results-proposed.RData` + CSVs in `output/us-all/tables/`

Additional Excel-ready tables from the proposed runner:

- `comparison_scalar_metrics.csv`
  - one row per setting
  - includes CRPS summaries and placeholder columns for `LOO-ELPD` / `WAIC`
- `comparison_uncertainty_summary.csv`
  - one row per setting
  - includes coverage targets/gaps plus PIT-based calibration summaries
- `comparison_calibration_bins.csv`
  - one row per setting x PIT bin
  - supports Excel histogram/reliability-style views
- `comparison_brier_split.csv`
  - one row per setting x event quantile x band type
  - reports same-threshold Brier splits for `all`, `below_threshold`, and `above_threshold`

Current limitation:

- `LOO-ELPD` / `WAIC` stay `NA` unless future result files save pointwise log-likelihood values.
- Existing result files expose posterior predictive draws (`fit[[d]]$yp`) but not pointwise log-likelihood arrays.

Legacy exploratory map/diagnostic plotting blocks are intentionally not auto-executed in these runners; they can be moved into dedicated plotting scripts if needed.
