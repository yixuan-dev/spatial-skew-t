# simstudy_prop

這個目錄是 proposed `prop model` 在 simstudy 的獨立工作區。

這裡的核心目標是：

- 使用和 `../simstudy/` 完全相同的模擬資料
- 跑 methods `1~5` 的 `prop` 版本
- 再和 `../simstudy/results/` 中的 Morris baseline 做直接比較

因此，`simstudy_prop` 不是另一套新的模擬資料來源，而是針對同一份 `simdata.RData` 的平行分析工作區。

- `run-prop.R`：執行 prop backend 的 simstudy-style 模擬（直接以這支腳本作為 CLI 入口）
- `run-prop-batch.R`：將 `settings_prop.csv` 的資料列轉成真正可執行的 batch runs
- `results-prop.R`：將 proposed outputs 和 baseline `../simstudy/results/` 做對照
- `analyze_mrts.R`：掃描 `fits/<setting>-<method>-<dataset>-p<K>.RData`，計算 Brier / Quantile score，分析 score 隨 MRTS rank `K` 的變化
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

## MRTS rank 分析（`analyze_mrts.R`）

針對同一個 `(setting, method, dataset)`、不同 MRTS rank `K` 的 fits 做 Brier / Quantile score 比較：

```
Rscript analyze_mrts.R                                # fits/ 中所有檔案
Rscript analyze_mrts.R --setting=4 --method=2 --dataset=1   # 只挑某個 block
```

輸入：`fits/<setting>-<method>-<dataset>-p<K>.RData`（讀取 `fit.1$yp`，並使用 `../../R/prop/auxfunctions.R` 中的 `BrierScore` / `QuantScore`）。

輸出：

- `output/tables/mrts_brier_long.csv`：每個 `(K, quantile)` 一列，含 threshold、Brier、quantile score、elapsed time
- `output/tables/mrts_brier_summary.csv`：每個 `K` 的 score 平均，並依 band 拆分（`bulk` = q ∈ [0.90, 0.95]、`tail` = q ≥ 0.98、`all` = 全部 quantile）
- `output/results/mrts_brier_analysis.RData`：`brier_long`、`brier_summary`、`brier_wide`（K × quantile 矩陣）
