# simstudy

這個目錄是 Morris baseline 模擬研究的主要工作區。涵蓋三個階段：

1. **資料生成**：`setup.R` / `setup_def.R` 產生 `simdata.RData` / `simdata_def.RData`
2. **模型擬合**：`run-settings.R` 用各種 method × MRTS K 組合對每個 dataset 跑 MCMC，把 fit 寫到 `results/` 或 `results_def/`
3. **後處理 pipeline**：`scores.R` (Stage 1) 從 fits 算 Brier / Quantile 分數，`tables.R` (Stage 2) 從快取產出 CSV 表格 + PDF

## 主要腳本

- `setup.R`、`setup_def.R`：產生 `simdata.RData`、`simdata_def.RData`
- `run-settings.R`：擬合 driver；method 1–8，可選 MRTS basis K（細節見 [run-settings.md](run-settings.md)）
- `scores.R`：**Stage 1**——從 `results<suffix>/<setting>-<method>-<dataset>[-K<mrts_k>].RData` 計算分數，輸出 `output/results/scores<setting><suffix>.RData`
- `tables.R`：**Stage 2**——讀 Stage 1 快取，輸出 `output/tables/` 下的 CSV 表格 + 預設 `output/plots/` 下的 PDF
- `mrts_cov_helpers.R`：CLI / 檔名 / seed / method catalog helper（兩個階段共用）
- `package_load.R`、`ar2_load.R`：Morris backend / AR(2) backend 套件載入
- `results-mrts-cov.R`：MRTS 共變量覆蓋率分析（與主 pipeline 獨立）
- `timeseriesplots.R`：時間序列預測診斷圖

shell launchers (`launch-*.sh`) 是 `run-settings.R` 的舊式批次入口；長期會被
manifest-style 替代。

## 資料來源

`scores.R`、`tables.R`、`run-settings.R` 預設找 `./simdata.RData`，
找不到時自動回退到 `../simstudy/simdata.RData`（在 simstudy 自身的目錄內這條
fallback 沒有作用，但保留是為了和 `simstudy_prop/` 的腳本對稱）。

要切換到其他資料集（例如 deformed-covariance 的 `simdata_def.RData`），加
`--data=<path>`：

```
Rscript scores.R --setting=1 --data=simdata_def.RData
Rscript tables.R --setting=1 --data=simdata_def.RData
```

腳本會先試 `./<path>`，若不存在則自動回退到 `../simstudy/<basename>`。

## Setting catalog（資料生成設定，對應 `--setting=<id>`）

`--setting` 是 `simdata.RData` 第 4 維的索引（`dim(y)[4]`）。Setup 對應如下
（見 [setup.R](setup.R)）：

| setting | 資料生成過程 |
|---|---|
| 1 | Gaussian |
| 2 | t, K = 1 |
| 3 | t, K = 5 |
| 4 | **Skew-t, K = 1, λ = 3** |
| 5 | Skew-t, K = 5, λ = 3 |
| 6 | Max-stable, ξ = 0.2 |
| 7 | Skew-t (setting 4) 但 threshold 以下做 exp 變換 |
| 8 | Brown-Resnick |
| 9 | Skew-t, K = 1, fixed AR(2), stronger serial dependence: φ=(0.80, -0.35) |
| 10 | Skew-t, K = 1, fixed AR(2), weaker serial dependence: φ=(0.12, -0.05) |

setting 4 / 5 是文件最常用的 skew-t 目標。`simdata_def.RData`（deformed
covariance）只有 setting 1–3，傳 `--setting=4` 會直接報「setting must be in 1..3」。

## Method catalog（分析方法，對應 `--methods=<spec>`）

Morris baseline 的 method 1–8：

- 1: Gaussian
- 2: Skew-t, K = 1
- 3: t, K = 1, threshold q(0.80)
- 4: Skew-t, K = 5
- 5: t, K = 5, threshold q(0.80)
- 6: Max-stable, threshold q(0.80)
- 7: Skew-t, K = 1 + AR(2) temporal (τ, z, knots：`temporaltau/z/w=TRUE`, `ar2_tau/z/w=TRUE`)
- 8: Skew-t, K = 5 + AR(2) temporal (同上三組 φ)

`scores.R` 預設 `--methods=1:5`（method 6 max-stable 因為 fit 物件結構不同，
通常另外處理；`--methods=1:8` 仍可加進來，缺少的參數欄位會被填成 NA）。

注意 setting 與 method 是兩條獨立的軸：例如 `--setting=4 --methods=1:5` 表示
「在 skew-t-K1 資料上跑 Gaussian / Skew-t / t / Skew-t / t 五個分析方法」。

## 目錄與輸出規則

```
results/                fits 來源（baseline + MRTS K）
results_def/            fits 來源（deformed-covariance 資料的 fits）

output/results/         Stage 1 / Stage 2 cache
output/tables/          Stage 2 CSV
output/plots/           Stage 2 PDF
```

依 `--data` 推導：

| 輸入 | fits 目錄 | 輸出 suffix |
|---|---|---|
| `simdata.RData`     | `results/`     | （無 suffix） |
| `simdata_def.RData` | `results_def/` | `_def` |

也就是說 `--data=simdata_def.RData` 不會和原本 `simdata.RData` 的輸出檔混在一起。

### 檔名規則（fits）

- baseline：`<results_dir>/<setting>-<method>-<dataset>.RData`
- MRTS：`<results_dir>/<setting>-<method>-<dataset>-K<mrts_k>.RData`

`scores.R` 把 mrts_k = 0 視為 baseline（無 `-K` 後綴），mrts_k > 0 對應到
`-K<mrts_k>` 變體。

### 檔名規則（output）

| 階段 | 檔名 |
|---|---|
| Stage 1 cache | `output/results/scores<setting><suffix>.RData` |
| Stage 2 long table | `output/tables/score_long<setting><suffix>.csv` |
| Stage 2 mean / median | `output/tables/score_mean<setting><suffix>.csv` |
| Stage 2 rel-vs-Gaussian | `output/tables/score_rel_gauss<setting><suffix>.csv` |
| Stage 2 best (method, K) | `output/tables/best_method_per_K<setting><suffix>.csv` |
| Stage 2 lambda 95% 覆蓋率 | `output/tables/lambda_coverage<setting><suffix>.csv` |
| Stage 2 aggregated objects | `output/results/simresults<setting><suffix>.RData` |
| Stage 2 plots | `output/plots/{bs,qs}_*-set<setting><suffix>-K<k>.pdf`、`lambda_ci_vs_dataset-set<setting><suffix>-method<m>-K<k>.pdf` |

## Post-fit pipeline (`scores.R` + `tables.R`)

擬合完成後，分數計算與表格產出是兩個獨立步驟。

### Stage 1: `scores.R`

```
Rscript scores.R --setting=<id>
                 [--data=<path>]
                 [--methods=<spec>]   default 1:5
                 [--datasets=<spec>]  default 1..nsets
                 [--mrts_k=<spec>]    default = auto-detect from results_dir
```

- `--setting` 必填，會以 `dim(y)[4]` 驗證
- `--mrts_k` 預設自動掃描 `results<suffix>/`：把 baseline (`<setting>-*-*.RData`) 視為 mrts_k = 0，把 `<setting>-*-*-K<K>.RData` 對應 mrts_k = K，最後取聯集

輸出：`output/results/scores<setting><suffix>.RData`，內容是一份 4 維陣列 +
參數區間：

```
quant.score, brier.score        [length(probs), n_datasets, n_methods, n_mrts_k]
beta.0/1/2, tau.alpha, tau.beta,
rho, nu, gamma, lambda           [length(intervals), n_datasets, n_methods, n_mrts_k]
elapsed_sec                      [n_datasets, n_methods, n_mrts_k]

probs       <- c(0.9, 0.91, ..., 0.995)        # 11 個
intervals   <- c(0.01, 0.025, ..., 0.99)        # 8 個
mrts_ks, datasets, methods, setting
data_path, data_suffix, fits_dir   # provenance
```

每 10 個 dataset checkpoint 一次。`fits_dir` 寫進 cache 是 provenance 用的，
不會影響 Stage 2 的輸出路徑。

### Stage 2: `tables.R`

```
Rscript tables.R --setting=<id> [--data=<path>] [--no-plots]
```

讀 Stage 1 cache，輸出：

- 5 份 CSV（見上面表格）
- 1 份 `simresults<setting><suffix>.RData`（彙整物件）
- 預設另出 `output/plots/` 下的 PDF（`--no-plots` 可關閉）

當 `--mrts_k=0`（只有 baseline）這種 `n_ks = 1` 的切片時，Stage 2 會自動跳過
mean-vs-K 系列圖（沒有意義），CSV 仍正常產出。

### 端到端範例

```powershell
# Setting 4 全部（auto-detect baseline + 所有 MRTS K）
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" .\scores.R --setting=4
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" .\tables.R --setting=4

# AR(2) 固定係數新設定（setting 9 / 10）+ AR(2) 分析法 7／8 smoke test
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" .\run-settings.R --setting=9 "1" 1 2 "(7,8)"
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" .\run-settings.R --setting=10 "1" 1 2 "(7,8)"
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" .\scores.R --setting=9 --methods="(7,8)" --datasets="1" --mrts_k=0
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" .\scores.R --setting=10 --methods="(7,8)" --datasets="1" --mrts_k=0
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" .\tables.R --setting=9 --no-plots
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" .\tables.R --setting=10 --no-plots

# Deformed-covariance dataset，3 個 setting
foreach ($s in 1..3) {
  & "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" .\scores.R --setting=$s --data=simdata_def.RData
  & "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" .\tables.R --setting=$s --data=simdata_def.RData
}

# 只看 baseline（不掃 MRTS）的 setting 4 切片
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" .\scores.R --setting=4 --mrts_k=0
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" .\tables.R --setting=4 --no-plots
```

## 與 `simstudy_prop/` 的對應

simstudy 與 simstudy_prop 共用同一份 `simdata.RData` / `simdata_def.RData`。
post-fit pipeline 結構相同，只差在第 4 維軸的名稱與 fits 目錄：

| | simstudy（這裡） | simstudy_prop |
|---|---|---|
| K 軸名稱 | `mrts_k` | `prop_k` |
| fits 目錄 | `results/` / `results_def/` | `fits/` / `fits_def/` |
| 檔名 K 後綴 | `-K<K>.RData` | `-p<K>.RData` |
| Stage 1 cache | `scores<setting><suffix>.RData` | `scores<setting>-prop<suffix>.RData` |
| 預設 method | 1:5（method 6 max-stable 另論） | 1:5 |

兩邊的 CSV 欄位結構完全一致（除了 `mrts_k` ↔ `prop_k`），跨資料夾比較時可以
直接 `rbind`。

## 例外：MRTS 共變量覆蓋率分析

`results-mrts-cov.R` 是針對「MRTS basis 增加時，beta 區間的覆蓋率」做專門分析，
輸入是 `comparison_mrts/`（不走 `output/`），輸出仍寫回 `comparison_mrts/`。
這條路徑刻意保留獨立，沒有併進 `scores.R` / `tables.R`。

如果之後想把覆蓋率比較也納入主 pipeline，那是 `tables.R` 加一份 CSV 的事，
不是另起腳本。
