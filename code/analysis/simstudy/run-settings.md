# run-settings.R 使用說明

此腳本可執行任意 data setting 的複現實驗，支援：

- `datasets` 向量輸入，例如 `1:5`、`20:28`
- `methods` 向量輸入，例如 `1:8`、`(2,3,8)`
- methods 1–5 的平行執行（PSOCK）
- method 6（max-stable）以 `ms_threads` 控制 C++ threads
- 以額外參數 `mrts_k` 將 methods 1–5 切換為指定 MRTS basis 個數的版本

## 命令列格式

`Rscript run-settings.R [--data=<path>|--data <path>] [--setting=<id>|--setting <id>] [datasets] [workers] [ms_threads] [methods] [mrts_k]`

- `--data`、`--setting` 僅接受放在 positional args 之前（兩者順序可互換）。
- `--data` 與 `--setting` 若有提供，需放在所有 positional args 之前（兩者順序可互換）。
- 若未提供 `--setting`，預設使用 `setting = 5`。
- 若未提供 `--data`，預設讀取 `./simdata.RData`。
- `mrts_k` 為選填；若提供，原本選到的 method 1–5 會改為只執行對應的 MRTS 版本（不再同時跑 baseline）。

## 參數說明

- `data`：資料檔路徑，對應 `--data`（預設 `./simdata.RData`）
- `setting`：單一整數 setting（有效範圍依 `simdata.RData` 的 setting 維度）
- `datasets`：dataset 向量表達式（範圍 1..50）
- `workers`：methods 1–5 的平行 worker 數
- `ms_threads`：method 6 的執行緒數
- `methods`：method 向量表達式（範圍 1..8）
- `mrts_k`：MRTS basis 個數，可為單一整數或向量；只會套用到 method 1–5

### `--data` 輸出目錄規則

- `--data` 預設：`./simdata.RData`，輸出目錄為 `results/`
- 若檔名為 `simdata_def.RData`，輸出目錄為 `results_def/`
- 規則為：`results` + 去掉副檔名後檔名中 `simdata` 之後的 suffix

## 支援的向量語法

- `1:5`
- `20:28`
- `c(1,2,6)`
- `(1,2,6)`，會自動轉為 `c(1,2,6)`

PowerShell 建議用字串傳入，避免 shell 先行解讀。

## 方法編號

- 1: Gaussian
- 2: Skew-t, K=1
- 3: t, K=1, threshold q(0.80)
- 4: Skew-t, K=5
- 5: t, K=5, threshold q(0.80)
- 6: Max-stable, threshold q(0.80)
- 7: Skew-t, K=1 + temporal AR(2) (`temporaltau/z/w=TRUE`, `ar2_tau/z/w=TRUE`)
- 8: Skew-t, K=5 + temporal AR(2) (同上三組 φ)

## MRTS 擴增規則

若 `mrts_k` 非空，runner 會把 `methods` 內選到的 method 1–5 改為：

- `1+mrts{K}`：method 1 加上 `K` 個 MRTS basis
- `2+mrts{K}`：method 2 加上 `K` 個 MRTS basis
- `3+mrts{K}`：method 3 加上 `K` 個 MRTS basis
- `4+mrts{K}`：method 4 加上 `K` 個 MRTS basis
- `5+mrts{K}`：method 5 加上 `K` 個 MRTS basis

例如：

- `methods="(1,4,6)"` 且 `mrts_k="15"`

實際執行的 method key 會是：

- `6`
- `1+mrts15`
- `4+mrts15`

若 `methods` 沒有包含 1–5，提供 `mrts_k` 不會新增任何 MRTS task，非 1–5 的方法（例如 6）仍照常執行。

## Seed 規則

- methods 1–5（含 MRTS 版本）：`method_id * 1000 + setting * 100 + dataset_id`
- method 6：`setting * 100 + dataset_id`

這可確保：

1. 同一個 `(dataset_id, method_id, mrts_k)` 重跑時可重現
2. 平行或序列執行順序不應改變該組合的亂數路徑
3. 不同 `mrts_k` 版本若 `dataset_id` 與 `method_id` 相同，會共用 seed（目前設計）

## 輸出檔名

baseline 方法輸出：

- `results/<setting>-<method_id>-<dataset_id>.RData`

MRTS 擴增版本輸出：

- `results/<setting>-<method_id>-<dataset_id>-K{K}.RData`

例如：

- `results/5-2-1.RData`
- `results/5-3-1-K15.RData`
- `results/5-2-28-K5.RData`

每次執行也會另外輸出：

- `results/run-plan-setting-<setting>.csv`

這份檔案會列出本次所有實際排入的 method key 與其參數。

## Windows PowerShell 範例

```powershell
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" ".\run-settings.R" --setting=5 "1:5" 4 2 "1:8"
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" ".\run-settings.R" --setting=5 "1:5" 4 2 "(1,4,6)" "15"
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" ".\run-settings.R" --setting=5 "1:5" 4 2 "(1,3,5)" "c(5,10,15)"
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" ".\run-settings.R" --data=.\simdata_def.RData --setting=5 "1:5" 4 2 "(7,8)"

# Fixed-phi AR(2) settings + methods 7／8 AR(2) (pilot)
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" ".\run-settings.R" --setting=9 "1" 1 2 "(7,8)"
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" ".\run-settings.R" --setting=10 "1" 1 2 "(7,8)"
```

## 平行策略說明

- methods 1–5（無 `mrts_k`）或 methods 1–5 的 MRTS 版本（有 `mrts_k`）：建立 `dataset × method_key` 任務網格，派發至 PSOCK workers
- method 6：依 dataset 逐一執行，內部用 `ms_threads`

若只想先 smoke test，建議：

- `datasets="1"`
- `methods="(1,4)"`
- `mrts_k="15"`
- `workers=1`

## MRTS vs NON-MRTS 彙整（results-mrts-cov.R）

當 `results/` 已有 baseline 與 MRTS 輸出檔時，可執行 `results-mrts-cov.R` 彙整比較。

常用環境變數：

- `SIMSTUDY_MRTS_RESULTS_DIR`：結果目錄（預設 `results`）
- `SIMSTUDY_MRTS_OUTPUT_DIR`：彙整輸出目錄（預設 `comparison_mrts`）
- `SIMSTUDY_MRTS_SETTINGS`：要彙整的 settings（例如 `4`、`c(4,5)`）
- `SIMSTUDY_MRTS_METHODS`：要比較的 baseline method id（預設 `1:5`）
- `SIMSTUDY_MRTS_DATASETS`：dataset 範圍（例如 `1:50`）
- `SIMSTUDY_MRTS_K`：要納入比較的 MRTS `k` 向量（例如 `c(10,15,20,25)`）

新增（分區比較與選 `k`）：

- `SIMSTUDY_MRTS_BULK_RANGE`：bulk 分位數區間，兩個機率值（預設 `c(0.90,0.95)`）
- `SIMSTUDY_MRTS_TAIL_MIN`：tail 起始分位數（預設 `0.98`）
- `SIMSTUDY_MRTS_FOCUS_QUANTILES`：重點分位數（預設 `c(0.95,0.98,0.99)`）
- `SIMSTUDY_MRTS_OBJECTIVE`：選 `k` 目標，`balanced` / `extreme-first` / `bulk-first`（預設 `balanced`）

`results-mrts-cov.R` 會將以下檔案寫到 `SIMSTUDY_MRTS_OUTPUT_DIR`（預設 `comparison_mrts/`）：

- `comparison_mrts_cov_paired.csv`：逐 setting × quantile 的 paired 比較
- `comparison_mrts_cov_summary.csv`：`family × k × quantile` 彙整
- `comparison_mrts_cov_summary_by_band.csv`：`family × k × (bulk/tail/other)` 彙整
- `comparison_mrts_cov_focus_quantiles.csv`：focus quantiles 的摘要
- `comparison_mrts_cov_best_k.csv`：依 objective 權重產生的最佳 `k` 排序
- `comparison_mrts_cov_report_config.csv`：本次彙整使用的設定快照
