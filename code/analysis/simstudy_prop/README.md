# simstudy_prop

這個目錄是 proposed `prop model` 在 simstudy 的獨立工作區。

這裡的核心目標是：

- 使用和 `../simstudy/` 完全相同的模擬資料
- 跑 methods `1~5` 的 `prop` 版本
- 再和 `../simstudy/results/` 中的 Morris baseline 做直接比較

因此，`simstudy_prop` 不是另一套新的模擬資料來源，而是針對同一份 `simdata.RData` 的平行分析工作區。

- `run-prop.R`：執行 prop backend 的 simstudy-style 模擬（直接以這支腳本作為 CLI 入口）
- `run-prop-batch.R`：將 `settings_prop.csv` 的資料列轉成真正可執行的 batch runs
- `scores-prop.R`：Stage 1，從 `fits<suffix>/` 計算 Brier / Quantile 分數，存到 `scores<setting>-prop<suffix>.RData`
- `tables-prop.R`：Stage 2，讀 `scores<setting>-prop<suffix>.RData`，產出 CSV 表格 + `simresults<setting>-prop<suffix>.RData` 彙整物件
- `plots-prop.R`：Stage 3，讀 `simresults<setting>-prop<suffix>.RData`，產出 `output/plots/` 下的 PDF 圖
- `prop_load.R`：載入 `../../R/prop` backend 與 simstudy 資料
- `prop_simstudy_helpers.R`：CLI / 檔名 / seed / method catalog helper
- `settings_prop.csv`：prop 專用的批次 manifest
- `run-prop.md`：`simstudy_prop` 的 canonical 使用說明

真正的 backend 原始碼維持在：

- `../../R/prop/`

資料檔預設不複製到這裡。腳本會優先找：

1. `./simdata.RData`
2. `../simstudy/simdata.RData`

因此如果沿用原本 simstudy 資料，直接在這個目錄執行即可。

要切換到其他資料集（例如 deformed-covariance 的 `simdata_def.RData`），可以加 `--data=<path>`，
腳本會先試 `./<path>`，若不存在則自動回退到 `../simstudy/<basename>`：

```
Rscript run-prop.R --data=simdata_def.RData --setting=1 1 1 1 1 5
```

目前 `prop` 版 method catalog 固定為：

- 1: Gaussian
- 2: skew-t, `K = 1`
- 3: t, `K = 1`, threshold `q(0.80)`
- 4: skew-t, `K = 5`
- 5: t, `K = 5`, threshold `q(0.80)`

預設輸出路徑如下：

- `fits/`：model fit `.RData` 與 run plan，例如 `fits/3-1-2-p20.RData`
- `output/results/`：載入 fits 後的 analysis `.RData`（含 `mrts_brier_analysis.RData`）
- `output/tables/`：比較表（含 `mrts_brier_long.csv`、`mrts_brier_summary.csv`）
- `output/plots/`：後續比較圖

## Post-fit 分析 pipeline (`scores-prop.R` → `tables-prop.R` → `plots-prop.R`)

擬合完成後，分數計算、表格產出、繪圖是三個獨立步驟，可分別重跑：

```
Rscript scores-prop.R --setting=4                       # Stage 1：算分數，存 .RData
Rscript tables-prop.R --setting=4                       # Stage 2：讀 .RData，輸出表格
Rscript plots-prop.R  --setting=4                       # Stage 3：讀彙整物件，輸出 PDF

Rscript scores-prop.R --setting=1 --data=simdata_def.RData
Rscript tables-prop.R --setting=1 --data=simdata_def.RData
Rscript plots-prop.R  --setting=1 --data=simdata_def.RData
```

Stage 1（`scores-prop.R`）讀 `fits<suffix>/<setting>-<method>-<dataset>-p<K>.RData` 中的 `fit.1`，使用 `../../R/ar2/auxfunctions.R` 的 `BrierScore` / `QuantScore`，輸出單一 `.RData` 快取：

- `output/results/scores<setting>-prop<suffix>.RData`：4 維分數陣列 `[probs, dataset, method, prop_k]` + 參數區間 + `elapsed_sec`

Stage 2（`tables-prop.R`）從 Stage 1 快取彙整：

- `output/tables/score_long<setting>-prop<suffix>.csv`（per-dataset 長表）
- `output/tables/score_mean<setting>-prop<suffix>.csv`（dataset 平均 + 中位數）
- `output/tables/score_rel_gauss<setting>-prop<suffix>.csv`（相對 Gaussian）
- `output/tables/best_method_per_K<setting>-prop<suffix>.csv`
- `output/tables/lambda_coverage<setting>-prop<suffix>.csv`
- `output/results/simresults<setting>-prop<suffix>.RData`（彙整物件，供 Stage 3 使用；內含 `lambda` 原始陣列）

Stage 3（`plots-prop.R`）讀 Stage 2 的 `simresults<setting>-prop<suffix>.RData`，
輸出 `output/plots/` 下的 PDF 圖（relative-score-by-quantile、mean-vs-K、
lambda 95% CI）。

更完整的選項說明請看 [run-prop.md](run-prop.md) 的「Post-fit pipeline」節。
