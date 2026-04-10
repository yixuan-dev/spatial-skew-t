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

Runners write the same key output filenames used in your current workflow:

- Morris runner: `us-all-results-0401.RData`, `us-all-results.RData`
- A runner: `us-all-results-combined.RData`, `us-all-results-a.RData`
- Proposed runner: `us-all-results-proposed.RData` + comparison CSVs

Legacy exploratory map/diagnostic plotting blocks are intentionally not auto-executed in these runners; they can be moved into dedicated plotting scripts if needed.
