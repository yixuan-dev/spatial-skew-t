# run-settings.R 使用說明

此腳本可執行任意 data setting 的複現實驗，支援：

- `datasets` 向量輸入，例如 `1:5`、`20:28`
- `methods` 向量輸入，例如 `1:6`、`(2,3,6)`
- methods 1–5 的平行執行（PSOCK）
- method 6（max-stable）以 `ms_threads` 控制 C++ threads
- 以額外參數 `mrts_k` 將 methods 1–5 切換為指定 MRTS basis 個數的版本

## 命令列格式

`Rscript run-settings.R [--setting=<id>|--setting <id>] [datasets] [workers] [ms_threads] [methods] [mrts_k]`

- 若提供 `--setting`，必須放在腳本名稱後的第一個參數位置。
- 若未提供 `--setting`，預設使用 `setting = 5`。
- `mrts_k` 為選填；若提供，原本選到的 method 1–5 會改為只執行對應的 MRTS 版本（不再同時跑 baseline）。

## 參數說明

- `setting`：單一整數 setting（有效範圍依 `simdata.RData` 的 setting 維度）
- `datasets`：dataset 向量表達式（範圍 1..50）
- `workers`：methods 1–5 的平行 worker 數
- `ms_threads`：method 6 的執行緒數
- `methods`：method 向量表達式（範圍 1..6）
- `mrts_k`：MRTS basis 個數，可為單一整數或向量；只會套用到 method 1–5

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
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" ".\run-settings.R" --setting=5 "1:5" 4 2 "1:6"
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" ".\run-settings.R" --setting=5 "1:5" 4 2 "(1,4,6)" "15"
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" ".\run-settings.R" --setting=5 "1:5" 4 2 "(1,3,5)" "c(5,10,15)"
```

## 平行策略說明

- methods 1–5（無 `mrts_k`）或 methods 1–5 的 MRTS 版本（有 `mrts_k`）：建立 `dataset × method_key` 任務網格，派發至 PSOCK workers
- method 6：依 dataset 逐一執行，內部用 `ms_threads`

若只想先 smoke test，建議：

- `datasets="1"`
- `methods="(1,4)"`
- `mrts_k="15"`
- `workers=1`
