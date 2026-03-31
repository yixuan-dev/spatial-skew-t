# run-setting5.R 使用說明

此腳本用於 data setting = 5（Skew-t, K=5, lambda=3）的複現實驗，支援：

- `datasets` 向量輸入（例如 `1:5`、`20:28`）
- `methods` 向量輸入（例如 `1:5`、`(1,2,6)`）
- methods 1–5 的平行執行（PSOCK）
- method 6（max-stable）以 `ms_threads` 控制 C++ threads

## 命令列格式

`Rscript run-setting5.R [datasets] [workers] [ms_threads] [methods]`

- `datasets`：dataset 向量表達式（範圍 1..50）
- `workers`：methods 1–5 的平行 worker 數
- `ms_threads`：method 6 的執行緒數
- `methods`：method 向量表達式（範圍 1..6）

### 支援的向量語法

- `1:5`
- `20:28`
- `c(1,2,6)`
- `(1,2,6)`（會自動轉為 `c(1,2,6)`）

> 建議在 PowerShell 以字串傳入（加雙引號），避免 shell 先行解讀。

## 方法編號

- 1: Gaussian
- 2: Skew-t, K=1
- 3: t, K=1, threshold q(0.80)
- 4: Skew-t, K=5
- 5: t, K=5, threshold q(0.80)
- 6: Max-stable, threshold q(0.80)

## Seed 規則（重點）

當 `datasets` 與 `methods` 改為向量時，每個 `(dataset, method)` 任務都在函數內獨立設定 seed：

- methods 1–5：`set.seed(method_id * 1000 + setting * 100 + dataset_id)`
- method 6：`set.seed(setting * 100 + dataset_id)`

這可確保：

1. 同一個 `(dataset_id, method_id)` 在重跑時可重現
2. 平行或序列執行順序不應改變該組合的亂數路徑
3. 只跑 methods 子集合時，不影響已執行方法本身的 seed 定義

## 輸出檔名

每次任務輸出為：

`results/5-<method_id>-<dataset_id>.RData`

例如：

- `results/5-2-1.RData`
- `results/5-6-28.RData`

## Windows PowerShell 範例

```powershell
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" ".\run-setting5.R" "1:5" 4 2 "1:6"
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" ".\run-setting5.R" "20:28" 3 2 "1:5"
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" ".\run-setting5.R" "1:5" 2 2 "(1,2,6)"
```

## 平行策略說明

- methods 1–5：建立 `dataset × method` 任務網格後，派發至 PSOCK workers。
- method 6：依 dataset 逐一執行（內部用 `ms_threads`）。

如果你只想先 smoke test，可先用：

- `datasets="1"`
- `methods="(2)"` 或 `methods="(6)"`
- `workers=1`
