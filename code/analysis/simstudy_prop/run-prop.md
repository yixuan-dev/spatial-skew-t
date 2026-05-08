# simstudy_prop `run-prop.R` / `run-prop-batch.R` 使用說明

這個目錄的工作目標，是在 **和 `code/analysis/simstudy` 完全相同的模擬資料** 上，
執行 methods `1~5` 的 `prop` 版本，並和 Morris baseline 做比較。

## 資料來源

`run-prop.R` 與 `results-prop.R` 預設會優先找：

1. `./simdata.RData`
2. `../simstudy/simdata.RData`

若未另外提供本地副本，預設直接使用 `../simstudy/simdata.RData`，
這樣可以確保 `simstudy_prop` 與 `simstudy` 的比較是建立在同一份模擬資料上。

要切換到其他資料集（例如 deformed-covariance 的 `simdata_def.RData`），
可以使用 `--data=<path>`：

- 先試 `./<path>`
- 不存在的話自動回退到 `../simstudy/<basename>`

因此 `--data=simdata_def.RData` 會直接讀到 `../simstudy/simdata_def.RData`，
不需要把檔案複製到本地。

## 命令列格式

```
Rscript run-prop.R [--data=<path>|--data <path>] \
                   [--setting=<id>|--setting <id>] \
                   [datasets] [workers] [ms_threads] [methods] [mrts_k]
```

`--data` 與 `--setting` 必須出現在 positional arguments 之前，順序不限。

這裡保留了和 `../simstudy/run-settings.R` 相同的 CLI 外型，
但語義有兩個重要差異：

- `methods` 只接受 `1~5`
- 最後一個位置參數 `mrts_k` 在這裡代表 `prop` 的 low-rank dimension / basis rank，而不是 MRTS covariates

## 方法編號

- 1: Gaussian
- 2: skew-t, `K = 1`
- 3: t, `K = 1`, threshold `q(0.80)`
- 4: skew-t, `K = 5`
- 5: t, `K = 5`, threshold `q(0.80)`

## 參數說明

- `--data`：可選，dataset 的檔案路徑。預設使用 `../simstudy/simdata.RData`
- `--setting`：單一整數 setting，範圍由載入的 dataset (`dim(y)[4]`) 決定
- `datasets`：dataset 向量表達式，例如 `1:5`、`(1,3,5)`
- `workers`：methods `1~5` 的 PSOCK worker 數
- `ms_threads`：保留同一 CLI 位置，目前 `prop` workflow 不另外使用
- `methods`：method 向量表達式，範圍 `1..5`
- `mrts_k`：`prop` basis rank，可為單一整數或向量

## 支援的向量語法

- `1:5`
- `c(1,2,5)`
- `(1,2,5)`，會自動轉成 `c(1,2,5)`

PowerShell 建議將向量規格放在引號中，避免 shell 先行解讀。

## 輸出位置 (`fits_dir`) 規則

`fits_dir` 的決定優先序：

1. 環境變數 `SIMSTUDY_PROP_FITS_DIR`（若有設定一律生效）
2. 從 `--data` 的 basename 推導：
   - `simdata.RData` → `fits/`
   - `simdata_def.RData` → `fits_def/`
   - 其他 `simdata<suffix>.RData` → `fits<suffix>/`
3. 預設 `fits/`

也就是說，當切換到 `simdata_def.RData` 時，輸出會自動分流到 `fits_def/`，
不會和原本 `fits/` 中的 baseline 結果混在一起。

## 目前 workflow 的資料與輸出分工

- `fits/`（或 `fits_def/`、`fits<suffix>/`）
  - model fit `.RData`，命名為 `<fits_dir>/<setting>-<method>-<dataset>-p<k>.RData`
  - `run-plan-setting-<setting>.csv`
- `output/results/`
  - 載入 fits 後產生的 analysis `.RData`
- `output/tables/`
  - comparison CSV tables
- `output/plots/`
  - 預留給比較圖

## Windows PowerShell 範例

```powershell
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" ".\run-prop.R" --setting=5 "1:5" 4 1 "1:5" "20"
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" ".\run-prop.R" --setting=3 "2" 2 1 "1" "20"
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" ".\run-prop.R" --setting=5 "1:3" 2 1 "(3,5)" "c(20,30)"

# Deformed-covariance dataset (simdata_def.RData), three settings
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" ".\run-prop.R" --data=simdata_def.RData --setting=1 "1:5" 4 1 "1:5" "20"
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" ".\run-prop.R" --data=simdata_def.RData --setting=2 "1:5" 4 1 "1:5" "20"
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" ".\run-prop.R" --data=simdata_def.RData --setting=3 "1:5" 4 1 "1:5" "20"

& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" ".\run-prop-batch.R" 1
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" ".\run-prop-batch.R" --run_id=phase3-02
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" ".\results-prop.R"
```

例如：

- `--setting=3 "2" 2 1 "1" "20"` 會寫出 `fits/3-1-2-p20.RData`
- `--setting=5 "1:3" 2 1 "(3,5)" "c(20,30)"` 會寫出 `fits/5-3-1-p20.RData`, `fits/5-5-1-p20.RData`, `fits/5-3-1-p30.RData`, ...
- `--data=simdata_def.RData --setting=1 "1" 1 1 "1" "20"` 會寫出 `fits_def/1-1-1-p20.RData`

## `settings_prop.csv` 批次執行

`run-prop-batch.R` 會把 `settings_prop.csv` 的資料列轉成真正的 batch runs。

- 列號是以資料列為準的 1-based 編號，不含 header
- 若不給 selector，會依 `enabled=yes` 與 `priority` 順序依序執行
- `runner_script` 欄位目前預期為 `run-prop.R`

批次腳本會另外輸出：

- `output/results/settings_prop_batch_selection.csv`
- `output/results/settings_prop_batch_status.csv`

## 比較目的

`results-prop.R` 的定位不是單純彙整 `prop` 自己的表現，而是回答：

- 在相同模擬資料下，`prop model` 是否可行？
- `prop` 相對於 Morris baseline 的 `QuantScore` / `BrierScore` 表現如何？
- 哪些 setting、哪種 method family、哪個 `prop_k` 比較值得往下追？

baseline 結果預設從 `../simstudy/results/` 讀取，
而 `prop` fits 預設從 `fits/`（或對應的 `fits<suffix>/`）讀取。
`results-prop.R` 會優先找新的 `-p<k>.RData` 命名，必要時也能回讀舊的 `-P<k>.RData` 檔案。
