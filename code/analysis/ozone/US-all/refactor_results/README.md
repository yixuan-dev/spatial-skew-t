# US-all results refactor

This folder splits the old monolithic result scripts into reusable modules and purpose-specific runners.

Metric formulas and field definitions are documented in `METRIC_CALCULATION.md`.

## 執行所有 settings 的比較

### 1. extend lane (`/ozone/US-all/results/`)

```ps
cd D:\Github\spatial-skew-t\code\analysis\ozone\US-all
$env:US_ALL_SUMMARY_DRAWS='5000'   # 全部 draws；若要預設 400，就把這行拿掉
Rscript refactor_results/12_run_extend_results.R
```

yp 是 5000 x 400 x 31，所以 5000 就是全量。輸出全部帶 `_extend` 後綴 (`us-all-results-extend.RData`、`comparison_*_extend.csv` / `.xlsx`)。

### 2. proposed lane (`/ozone_prop/US-all/results/`)

```ps
cd D:\Github\spatial-skew-t\code\analysis\ozone\US-all
Rscript refactor_results/14_run_prop_results.R
```

只跑 10 個 ozone_prop fits，再加上 `/ozone/US-all/results/us-all-1.RData` 當 Gaussian baseline；輸出全部帶 `_proposed` 後綴。

### 3. pool — extend + proposed 合併比較

```ps
cd D:\Github\spatial-skew-t\code\analysis\ozone\US-all
Rscript refactor_results/15_pool_extend_prop.R
```

不會重新計算 score；只是把 `us-all-results-extend.RData` 與 `us-all-results-proposed.RData` 的 score 陣列在 setting 軸上拼接，extend 保留原本的 ID（1-209），ozone_prop fits 移到 unified ID 301-310。輸出全部帶 `_pool` 後綴。如果要改 prop 的偏移，調整 `15_pool_extend_prop.R` 裡的 `prop_id_base` 即可。

## 產生 MRTS 重點圖

先跑完 `12_run_extend_results.R` 之後，再執行：

```ps
cd D:\Github\spatial-skew-t\code\analysis\ozone\US-all
$env:US_ALL_PLOT_SUMMARY_DRAWS='400'  # optional；mean diagnostics 預設沿用 400 draws
Rscript refactor_results/13_plot_extend_results.R
```

這個 runner 會額外寫出：

- `output/us-all/tables/mrts_mean_diagnostics.csv`
- `output/us-all/plots/mrts_extreme_delta_profiles.png`
- `output/us-all/plots/mrts_extreme_split_brier.png`
- `output/us-all/plots/mrts_mean_error_gap.png`
- `output/us-all/plots/mrts_mean_bias.png`
- `output/us-all/plots/mrts_tail_mean_tradeoff.png`

圖的重點：

- `mrts_extreme_delta_profiles`
  - 看 AR2 / MRTS 在高 quantile (`0.90` 到 `0.995`) 相對同 basis baseline 的 tail score 改變
- `mrts_extreme_split_brier`
  - 把極端事件拆成 `above_threshold` 與 `below_threshold`
  - 分開看 miss 與 false alarm 的代價
- `mrts_mean_error_gap`
  - 看 MRTS 在平均層級 prediction 上對 `RMSE`、`MAE`、`|Bias|` 的優勢或劣勢
- `mrts_mean_bias`
  - 看 MRTS 是否傾向整體高估或低估平均值
- `mrts_tail_mean_tradeoff`
  - 直接把平均層級 `RMSE` 變化和 tail Brier 變化放在同一張圖
  - 左下角代表平均值與極端值都改善

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
- `03_plot_results.R`
  - MRTS-focused plotting helpers
  - mean-level diagnostics for predictive means
  - tail-vs-mean trade-off plots

## Runner files

| runner                          | scope                                                           | reads                                                                            | writes (`output/us-all/`)                                                                |
|---------------------------------|-----------------------------------------------------------------|----------------------------------------------------------------------------------|------------------------------------------------------------------------------------------|
| `10_run_morris_results.R`       | Morris baseline only (legacy `us-all-results.R`)                | `code/analysis/ozone/US-all/results/us-all-{i}.RData`                            | `results/us-all-results.RData`, `results/us-all-results-0401.RData`                      |
| `11_run_a_results.R`            | CMAQ-vs-no-CMAQ paired lane                                     | `code/analysis/ozone/US-all/results/`                                            | `results/us-all-results-a.RData`, `results/us-all-results-combined.RData`                |
| `12_run_extend_results.R`     | extend lane: Morris + AR2 + MRTS                                 | `code/analysis/ozone/US-all/results/us-all-{i}.RData`                            | `results/us-all-results-extend.RData`, `tables/comparison_*_extend.csv` / `.xlsx`        |
| `13_plot_extend_results.R`    | MRTS-focused plots (extend lane)                                | `results/us-all-results-extend.RData`                                            | `plots/mrts_*.png`, `tables/mrts_mean_diagnostics.csv`                                   |
| `14_run_prop_results.R`         | proposed lane: 10 ozone_prop fits + Gaussian baseline reused     | `code/analysis/ozone_prop/US-all/results/ozone-prop-{i}.RData` + `us-all-1.RData` | `results/us-all-results-proposed.RData`, `tables/comparison_*_proposed.csv` / `.xlsx`     |
| `15_pool_extend_prop.R`         | pool: extend ∪ proposed in one comparison (no re-scoring)       | `results/us-all-results-extend.RData` + `results/us-all-results-proposed.RData`   | `results/us-all-results-pool.RData`, `tables/comparison_*_pool.csv` / `.xlsx`            |

The historical name *"proposed"* in runner `12` referred to the AR2/MRTS extensions vs the Morris baseline, all within `/ozone/US-all/results/`. After the new ozone_prop lane was added, runner `12` and reader `13` were both patched (and the existing on-disk tables renamed) to use the `_extend` suffix, so `*_proposed` now unambiguously means the ozone_prop lane. The top-level wrapper [us-all-results-proposed.R](../us-all-results-proposed.R) still sources `12_*` (so its filename is now misleading — feel free to rename to `us-all-results-extend.R`).

## Output compatibility notes

All runner outputs now live together under `output/us-all/` with simple subfolders:

- `results/` for `.RData` outputs
- `tables/` for comparison CSVs and summaries
- `plots/` for figures
- `logs/` for logs

The runners keep the same key filenames used in your current workflow, but write them inside that shared output tree:

- Morris runner: `output/us-all/results/us-all-results-0401.RData`, `output/us-all/results/us-all-results.RData`
- A runner: `output/us-all/results/us-all-results-combined.RData`, `output/us-all/results/us-all-results-a.RData`
- Extend runner (12): `output/us-all/results/us-all-results-extend.RData` + `tables/comparison_*_extend.csv` / `.xlsx`
- Proposed runner (14): `output/us-all/results/us-all-results-proposed.RData` + `tables/comparison_*_proposed.csv` / `.xlsx`
- Pool runner (15): `output/us-all/results/us-all-results-pool.RData` + `tables/comparison_*_pool.csv` / `.xlsx`

See `output/us-all/README.md` for the full suffix-to-source table.

Additional Excel-ready tables from the proposed runner:

- `comparison_top2.csv`
  - long-format top-2 table across all comparable metrics
- `comparison_top2.xlsx`
  - `all_metrics` sheet plus one sheet per metric
- `comparison_scalar_metrics.csv`
  - one row per setting
  - includes MSPE, MAPE, CRPS summaries and placeholder columns for `LOO-ELPD` / `WAIC`
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
| 相對 Gaussian 排名   | `rel_to_gaussian`   | 真正拿來排名的相對分數 | 原始分數本身  | `brier`, `quantile`, `crps`, `mspe`, `mape`, `brier_split`                                                    |
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
