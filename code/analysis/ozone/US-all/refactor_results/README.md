# US-all results refactor

This folder splits the old monolithic result scripts into reusable modules and purpose-specific runners.

Metric formulas and field definitions are documented in `METRIC_CALCULATION.md`.

## math notation gaussian 95% confidence interval

$$\hat{X} \pm  z_{0.95} \cdot \frac{\sigma(\hat{X})}{\sqrt{n}}$$

, where $z_{0.95} \approx 1.96$.

## 執行所有 settings 的比較

```ps
cd D:\Github\spatial-skew-t\code\analysis\ozone\US-all
$env:US_ALL_SUMMARY_DRAWS='5000'   # 全部 draws；若要預設 400，就把這行拿掉
Rscript refactor_results/12_run_proposed_results.R
```

yp 是 5000 x 400 x 31，所以 5000 就是全量。

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
  - classification diagnostics (TP/TN/FP/FN, accuracy, precision, recall, specificity, F1)
- `02_comparison_tables.R`
  - full comparison table
  - top-2 ranking table for all comparable metrics
  - Excel workbook export for per-metric top-2 sheets
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

- `comparison_top2.csv`
  - long-format top-2 table across all comparable metrics
- `comparison_top2.xlsx`
  - `all_metrics` sheet plus one sheet per metric
- `comparison_scalar_metrics.csv`
  - one row per setting
  - includes CRPS summaries and placeholder columns for `LOO-ELPD` / `WAIC`
- `comparison_classification_metrics.csv`
  - one row per setting x event quantile
  - includes confusion-matrix counts plus `accuracy`, `precision`, `recall`, `specificity`, and `F1`
- `comparison_uncertainty_summary.csv`
  - one row per setting
  - includes coverage targets/gaps plus PIT-based calibration summaries
- `comparison_calibration_bins.csv`
  - one row per setting x PIT bin
  - supports Excel histogram/reliability-style views
- `comparison_brier_split.csv`
  - one row per setting x event quantile x band type
  - reports same-threshold Brier splits for `all`, `below_threshold`, and `above_threshold`

`comparison_top2` 欄位解讀：

| 情況                 | `ranking_basis`     | `rank_value`           | `score_value` | 代表指標                                                                                                      |
| -------------------- | ------------------- | ---------------------- | ------------- | ------------------------------------------------------------------------------------------------------------- |
| 相對 Gaussian 排名   | `rel_to_gaussian`   | 真正拿來排名的相對分數 | 原始分數本身  | `brier`, `quantile`, `crps`, `brier_split`                                                                    |
| 原始分數直接排名     | `raw_score`         | 真正拿來排名的原始分數 | 原始分數本身  | `accuracy`, `precision`, `recall`, `specificity`, `f1`, `pit_ks`, `pit_uniformity_mae`, `pit_uniformity_rmse` |
| 與 target 的距離排名 | `abs_gap_to_target` | 與理想值的絕對差       | 原始分數本身  | `coverage`, `pit_mean`, `pit_variance`                                                                        |

因此：

- `rank_value` 一律是 Top-2 排名依據
- `score_value` 一律是該模型的原始指標值
- 若 `ranking_basis = rel_to_gaussian`，就會出現 `rank_value` 和 `score_value` 不同
- 若 `ranking_basis = raw_score`，通常 `rank_value = score_value`
- 若 `ranking_basis = abs_gap_to_target`，通常 `rank_value = |score_value - target_value|`

Current limitation:

- `LOO-ELPD` / `WAIC` stay `NA` unless future result files save pointwise log-likelihood values.
- Existing result files expose posterior predictive draws (`fit[[d]]$yp`) but not pointwise log-likelihood arrays.

Legacy exploratory map/diagnostic plotting blocks are intentionally not auto-executed in these runners; they can be moved into dedicated plotting scripts if needed.
