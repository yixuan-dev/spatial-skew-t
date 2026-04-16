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
  - default quantile grids
- `01_score_engine.R`
  - result-file scoring loop
  - prediction contract checks
  - summary/relative-score computation
- `02_comparison_tables.R`
  - full comparison table
  - top-2 ranking table
  - paired same-basis table

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

Legacy exploratory map/diagnostic plotting blocks are intentionally not auto-executed in these runners; they can be moved into dedicated plotting scripts if needed.
