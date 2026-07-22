# US-all runner guide (`us-all-run.R`)

本文件聚焦在 **實驗設定（settings）如何被控制**，以及 **哪些 setting 之間可直接比較**。

---

## 1) Quick start

先切到目錄：

```powershell
Set-Location "d:\Github\spatial-skew-t\code\analysis\ozone\US-all"
```

建議在 Windows 使用明確的 Rscript 路徑：

```powershell
$cand = Get-ChildItem 'C:\Program Files\R' -Directory | Sort-Object Name -Descending | Select-Object -First 1
$rs = Join-Path $cand.FullName 'bin\Rscript.exe'
& $rs us-all-run.R 114
```

---

## 2) 核心觀念：一個 runner，全部 setting 走同一流程

`us-all-run.R` 會讀取 `settings.csv`，用同一套流程跑不同設定。
因此你不需要再維護 `us-all-xxx.R` wrapper scripts。

```powershell
# 單一 setting
& $rs us-all-run.R 114

# 多個 setting
& $rs us-all-run.R 101 102 103

# 範圍
& $rs us-all-run.R 1:124

# MRTS settings
& $rs us-all-run.R 201:209
```

---

## 3) 實驗設定（`settings.csv`）與可比較性

### 3.1 欄位說明表

| 欄位      | 型別/值域                                  | 在 runner 的用途                                       |
| --------- | ------------------------------------------ | ------------------------------------------------------ |
| `setting` | 字串（可含數字+字母，如 `5a`）             | 設定 ID；seed 取前綴數字（`5a -> 5`）                  |
| `method`  | `gaussian` / `t` / `skew-t` / `max-stable` | 映射模型 lane；`skew-t` 會映射為 `method=t, skew=TRUE` |
| `knots`   | 正整數                                     | 對應 `nknots`                                          |
| `thresh`  | 數值                                       | 對應 `thresh.all`                                      |
| `CMAQ`    | `yes` / `no`                               | `yes` 用完整 `X`；`no` 只用 `X[,,1]`（intercept）      |
| `TS`      | `yes` / `no`                               | `yes` 啟用 `temporaltau/temporalw/temporalz`           |
| `ar2`     | `yes` / 空白                               | 主要為設定標記；實際 backend 由環境變數控制            |
| `mrts`    | 空白或正整數                               | 正整數 `k` 時啟用 MRTS covariates                      |

### 3.2 常用設定區段

| 區段      | 主要特性                                       | 常見用途            |
| --------- | ---------------------------------------------- | ------------------- |
| `1:74`    | 基本 US-all 組合（多種 `method/knots/thresh`） | baseline 與主比較   |
| `101:124` | `ar2=yes` 標記群                               | AR2 對照實驗        |
| `201:209` | `mrts = 5/10/15`                               | MRTS covariate 實驗 |

### 3.3 實驗「控制變因」表（跨 setting 共用）

> 下列條件在 runner 內固定或規則固定，構成可比較實驗的共同基礎。

| 控制項         | 固定內容（`us-all-run.R`）                                        | 可比較性意義                          |
| -------------- | ----------------------------------------------------------------- | ------------------------------------- |
| 資料來源       | 同一組 `Y/X/S/cv.lst`（由 `package_load.R` 或 `ar2_load.R` 載入） | 不同 setting 使用同一資料母體         |
| CV 切分        | 固定 2 folds，索引由 `cv.lst` 提供                                | 在同一切分下比較                      |
| Seed 規則      | `set.seed(setting_prefix * 100 + fold)`                           | 隨機機制一致（種子值隨 setting 變化） |
| MCMC 排程規則  | 由 `RUN_MODE` 映射（dev/prod）                                    | 比較時需固定 `RUN_MODE`               |
| 參數初始化骨架 | `gamma/rho/nu` 初始化與上下界、`min.s/max.s` 固定策略             | 降低非目標超參數干擾                  |
| 輸出結構       | `results/us-all-<setting>.RData`，每設定含兩個 fold 結果          | 後處理可一致比較                      |

### 3.4 一眼看懂：所有 setting 都相同的內容

| 類別                      | 所有 setting 都相同嗎？ | 固定內容                                                                 |
| ------------------------- | ----------------------- | ------------------------------------------------------------------------ |
| Cross validation fold 數  | 是                      | 固定 2 folds（`val = 1, 2`）                                             |
| Cross validation 切分來源 | 是                      | 都使用同一個 `cv.lst`                                                    |
| Train/Validation 切分邏輯 | 是                      | 每個 setting 都用相同 fold index 去切 `Y/X/S`                            |
| 每個 setting 的輸出結構   | 是                      | 都存成 `fit[[1]], fit[[2]]` 到 `us-all-<setting>.RData`                  |
| MCMC 呼叫骨架             | 是                      | 都用同一個 `call_common` 架構呼叫 `mcmc()` / `mcmc_ar2()`                |
| MCMC 固定參數（見 5.5）   | 大多是                  | 例如 `keep.knots=FALSE`, `iterplot=FALSE`, `gamma/rho/nu` 初始化與上下界 |

### 3.5 比較前檢查表（A/B 設定）

| 檢查項目              | 建議是否相同                     | 目的                          |
| --------------------- | -------------------------------- | ----------------------------- |
| `US_ALL_RUN_MODE`     | 必須相同                         | 避免 `iters/burn/update` 差異 |
| `US_ALL_MCMC_BACKEND` | 必須相同（除非正在比較 backend） | 避免混入 backend 影響         |
| `method`              | 通常相同                         | 要比較方法效果時才改          |
| `knots`               | 通常相同                         | 隔離 knot 數效果              |
| `thresh`              | 通常相同                         | 隔離 threshold 效果           |
| `CMAQ`                | 建議相同                         | 避免 covariate 維度變化       |
| `TS`                  | 比較 TS 效果時才改               | 隔離 temporal 效果            |
| `ar2`/backend         | 比較 AR2 效果時才改              | 隔離 AR2 效果                 |
| `mrts`                | 比較 MRTS 效果時才改             | 隔離 MRTS 效果                |

### 3.6 可直接對照的 setting 比較表（推薦）

#### A) 比較 TS（Temporal）效果：`TS=no -> TS=yes`

| knots | `skew-t, thresh=0` | `skew-t, thresh=50` | `t, thresh=75` |
| ----: | ------------------ | ------------------- | -------------- |
|     1 | `3 -> 51`          | `4 -> 52`           | `5 -> 53`      |
|     5 | `7 -> 54`          | `8 -> 55`           | `9 -> 56`      |
|     6 | `33 -> 57`         | `38 -> 58`          | `43 -> 59`     |
|     7 | `34 -> 60`         | `39 -> 61`          | `44 -> 62`     |
|     8 | `35 -> 63`         | `40 -> 64`          | `45 -> 65`     |
|     9 | `36 -> 66`         | `41 -> 67`          | `46 -> 68`     |
|    10 | `11 -> 69`         | `12 -> 70`          | `13 -> 71`     |
|    15 | `15 -> 72`         | `16 -> 73`          | `17 -> 74`     |

#### B) 比較 AR2 backend 效果：`TS baseline -> AR2`

| 基準設定（`TS=yes`, 非 AR2） | AR2 對照設定 | 映射規則                        | 控制條件                           |
| ---------------------------- | ------------ | ------------------------------- | ---------------------------------- |
| `51:74`                      | `101:124`    | `setting_ar2 = setting_ts + 50` | `method/knots/thresh/CMAQ/TS` 對齊 |

例：`51 -> 101`, `52 -> 102`, ..., `74 -> 124`。

#### C) 比較 MRTS `k` 效果：`no MRTS -> k=5/10/15`

| 比較主題                     | no MRTS（基準） | MRTS `k=5` | MRTS `k=10` | MRTS `k=15` | 控制條件                           |
| ---------------------------- | --------------- | ---------- | ----------- | ----------- | ---------------------------------- |
| `skew-t, knots=1, thresh=0`  | `51`            | `201`      | `204`       | `207`       | `CMAQ=yes`, `TS=yes`, backend 相同 |
| `skew-t, knots=1, thresh=50` | `52`            | `202`      | `205`       | `208`       | 同上                               |
| `t, knots=1, thresh=75`      | `53`            | `203`      | `206`       | `209`       | 同上                               |

### 3.7 `mrts` 欄位規則

- 空白/缺值：不加 MRTS covariates。
- 正整數 `k`：每個 fold 建立 MRTS basis 並追加到 `x` covariates。

目前 runner 會：

1. 用 `autoFRK::mrts(S_train, k, x = S_pred)` 建 train/pred basis（含 fallback）。
2. 移除 near-constant basis 欄位（`sd <= 1e-10`）。
3. 把保留 basis append 到 `X.o / X.p`。

---

## 4) 環境變數

### `US_ALL_MCMC_BACKEND`

- `legacy`（預設，使用 `mcmc()`）
- `ar2`（使用 `mcmc_ar2()`）

```powershell
$env:US_ALL_MCMC_BACKEND = "legacy"
```

### `US_ALL_RUN_MODE`

- `dev`：`iters=2000`, `burn=1000`, `update=200`
- `prod`：`iters=30000`, `burn=25000`, `update=500`（預設）

```powershell
$env:US_ALL_RUN_MODE = "dev"
```

### `US_ALL_RESULTS_DIR`

結果輸出路徑（預設 `results_new`）。不存在會自動建立。

```powershell
$env:US_ALL_RESULTS_DIR = "results_mrts_cov_dev"
```

### `US_ALL_SETTINGS`

可用環境變數指定 settings（優先於命令列參數）：

```powershell
$env:US_ALL_SETTINGS = "201:209"
& $rs us-all-run.R
```

---

## 5) MCMC 輸入參數表（`us-all-run.R` 實際呼叫）

### 5.1 `RUN_MODE` 對應迭代控制

| `US_ALL_RUN_MODE` | `iters` | `burn` | `update` | 用途             |
| ----------------- | ------: | -----: | -------: | ---------------- |
| `dev`             |    2000 |   1000 |      200 | 快速 smoke test  |
| `prod`            |   30000 |  25000 |      500 | 正式實驗（預設） |

### 5.2 `run_mcmc()` 核心輸入（`call_common`）

| 參數                      | 來源                        | 備註                                          |
| ------------------------- | --------------------------- | --------------------------------------------- |
| `y`, `s`, `x`             | fold 訓練資料               | `y.o`, `S.o`, `X.o`                           |
| `x.pred`, `s.pred`        | fold 驗證資料               | `X.p`, `S.p`                                  |
| `method`                  | `settings.csv::method` 映射 | `skew-t -> method=t`                          |
| `skew`                    | `settings.csv::method` 映射 | `skew-t -> TRUE`                              |
| `keep.knots`              | 固定值                      | `FALSE`                                       |
| `thresh.all`              | `settings.csv::thresh`      | 數值閾值                                      |
| `thresh.quant`            | 固定值                      | `FALSE`                                       |
| `nknots`                  | `settings.csv::knots`       | knot 數                                       |
| `iters`, `burn`, `update` | `US_ALL_RUN_MODE`           | 見 5.1                                        |
| `iterplot`                | 固定值                      | `FALSE`                                       |
| `beta.init`               | 載入物件                    | 來自 `package_load.R` / `ar2_load.R`          |
| `tau.init`                | 方法相依                    | `gaussian` 用 `tau.init_default`；否則 `0.05` |
| `gamma.init`              | 固定值                      | `0.5`                                         |
| `rho.init`, `rho.upper`   | 固定值                      | `1`, `5`                                      |
| `nu.init`, `nu.upper`     | 固定值                      | `0.5`, `10`                                   |
| `min.s`, `max.s`          | 固定值                      | `c(-2.25,-1.55)` / `c(2.35,1.30)`             |

### 5.3 條件式參數（依 `TS` 與 backend）

| 條件                      | 追加參數                                               | 說明               |
| ------------------------- | ------------------------------------------------------ | ------------------ |
| `TS=yes`                  | `temporaltau=TRUE`, `temporalw=TRUE`, `temporalz=TRUE` | 啟用時間向度更新   |
| backend=`ar2` 且 `TS=yes` | `ar2_tau=TRUE`, `ar2_w=TRUE`, `ar2_z=TRUE`             | AR2 參數與 TS 同步 |
| backend=`ar2` 且 `TS=no`  | `ar2_tau=FALSE`, `ar2_w=FALSE`, `ar2_z=FALSE`          | 關閉 AR2 時序項    |

### 5.4 `max-stable` lane 主要輸入

| 參數                      | 值/來源                                       |
| ------------------------- | --------------------------------------------- |
| `y`                       | `t(y.o)`                                      |
| `x`, `xp`                 | `t(X.o[,,2])`, `t(X.p[,,2])`（需 `CMAQ=yes`） |
| `s`, `sp`                 | `S.o`, `S.p`                                  |
| `thresh`                  | `settings.csv::thresh`                        |
| `knots`                   | `S.o`                                         |
| `iters`, `burn`, `update` | 同 `RUN_MODE` 映射                            |
| `thin`                    | 固定 `1`                                      |

### 5.5 `mcmc()` 參數：哪些在所有 setting 都一樣？

> 下面把 `run_mcmc` 參數分成「跨 setting 固定」與「隨 setting 改變」。

| 參數                                    | 跨 setting 是否固定 | 固定/變動規則                                  |
| --------------------------------------- | ------------------- | ---------------------------------------------- |
| `keep.knots`                            | 固定                | 永遠 `FALSE`                                   |
| `thresh.quant`                          | 固定                | 永遠 `FALSE`                                   |
| `iterplot`                              | 固定                | 永遠 `FALSE`                                   |
| `gamma.init`                            | 固定                | 永遠 `0.5`                                     |
| `rho.init`, `rho.upper`                 | 固定                | 永遠 `1`, `5`                                  |
| `nu.init`, `nu.upper`                   | 固定                | 永遠 `0.5`, `10`                               |
| `min.s`, `max.s`                        | 固定                | 永遠 `c(-2.25,-1.55)` / `c(2.35,1.30)`         |
| `iters`, `burn`, `update`               | 條件固定            | 在同一 `RUN_MODE` 下固定；`dev` 與 `prod` 不同 |
| `method`, `skew`                        | 變動                | 由 `settings.csv::method` 映射                 |
| `nknots`                                | 變動                | 由 `settings.csv::knots` 決定                  |
| `thresh.all`                            | 變動                | 由 `settings.csv::thresh` 決定                 |
| `x`, `x.pred`                           | 變動                | 受 `CMAQ` 與 `mrts` 影響（covariate 維度可變） |
| `temporaltau`, `temporalw`, `temporalz` | 條件變動            | `TS=yes` 時為 `TRUE`                           |
| `ar2_tau`, `ar2_w`, `ar2_z`             | 條件變動            | backend=`ar2` 時依 `TS` 切換                   |
| `tau.init`                              | 條件變動            | `gaussian` 用 `tau.init_default`；其餘 `0.05`  |

---

## 6) 推薦執行範例

### A. 快速測試 MRTS 設定（dev）

```powershell
Set-Location "D:\Github\spatial-skew-t\code\analysis\ozone\US-all"
$cand = Get-ChildItem 'C:\Program Files\R' -Directory | Sort-Object Name -Descending | Select-Object -First 1
$rs = Join-Path $cand.FullName 'bin\Rscript.exe'

$env:US_ALL_RUN_MODE = "dev"
$env:US_ALL_MCMC_BACKEND = "legacy"
$env:US_ALL_RESULTS_DIR = "results_mrts_cov_dev"

& $rs us-all-run.R 201:209
```

### B. 只跑單一 MRTS 設定

```powershell
$env:US_ALL_RUN_MODE = "dev"
$env:US_ALL_RESULTS_DIR = "results_mrts_cov_dev_201"
& $rs us-all-run.R 201
```

### C. 清空環境變數回預設

```powershell
Remove-Item Env:\US_ALL_SETTINGS -ErrorAction SilentlyContinue
Remove-Item Env:\US_ALL_MCMC_BACKEND -ErrorAction SilentlyContinue
Remove-Item Env:\US_ALL_RUN_MODE -ErrorAction SilentlyContinue
Remove-Item Env:\US_ALL_RESULTS_DIR -ErrorAction SilentlyContinue
```

---

## 7) 平行執行（PowerShell 7）

可用 `run-parallel-ps7.ps1`。目前檔案中預設設定是：

- `$settings = (1..74) | Where-Object { $_ -ne 2 }`
- `$backend = "ar2"`
- `$runMode = "prod"`

若你要跑 MRTS 201~209，請把腳本中的 `$settings` 改成：

```powershell
$settings = 201..209
```

然後再執行：

```powershell
.\run-parallel-ps7.ps1
```

---

## 8) 常見問題

### Q1. 一定要建立 `us-all-201.R` 這類檔案嗎？

不用。直接 `us-all-run.R 201` 或 `us-all-run.R 201:209` 就可以。

### Q2. 為什麼看到 `MRTS source: autoFRK::mrts` 與欄位 drop 訊息？

這是正常訊息。runner 會在每個 fold 建 basis 並移除 near-constant 欄位，避免無訊號欄位造成數值問題。

### Q3. max-stable 可以和 MRTS 同時用嗎？

目前不支援；設定若同時是 `max-stable` 且 `mrts` 有值，runner 會停止並報錯。

---

## 9) 最短命令速查

```powershell
# dev + MRTS sweep
$env:US_ALL_RUN_MODE = "dev"
$env:US_ALL_RESULTS_DIR = "results_mrts_cov_dev"
& $rs us-all-run.R 201:209
```

---

## 10) Cincinnati 地圖管線（setting 204，2026-07-22 新增）

與 CV runner 分離的「全資料配適 → 網格預測 → 地圖」三段式，仿 `us-all-full-71.R` 系列：

| 階段 | 腳本 | 輸出 |
| --- | --- | --- |
| 配適 | `us-all-full-204.R` | `results/us-all-full-204.RData`（含 `fit/S.o/X.o/S.p/X.p/mrts_meta/cincy`） |
| 健康檢查 | `check-full-204.R` | 斷言失敗即 exit 1，**通過前不得進預測** |
| 預測 | `predict-cincy-204.R` | `us-all-pred-cincy-204.RData` |
| 地圖 | `make-map-cincy-204.R` | `plots/cincy-204-exceed75.pdf` |

- 三支皆吃 `dev` 參數/`US_ALL_RUN_MODE=dev`，dev 檔案帶 `-dev` 後綴，不會覆寫正式檔。
- 預測窗口：Cincinnati (1.0681, −0.0257) ±0.20（±200 km，34×34 網格）。
- MRTS basis 由 `mrts_basis.R` 提供（自 `us-all-run.R` 抽出共用）；配適時一併建好
  train/pred 兩側並存入 fit 檔，預測端**不重建**。

### ⚠️ z.init 修正的可比性斷點（2026-07-22）

`code/R/mcmc_cont_lambda.R`（legacy backend）的 `z.init` 預設值已從 `0` 修正為
`NULL → 0.6745/sqrt(tau.init)`（同 `ar2`/`prop`；詳見 `tex/z_init_bug`）。
在此之前所有 legacy 的 **skew-t + TS** 配適（settings 51、52、54…、201/202/204/205/207/208）
z 全程凍結為 0、λ 隨先驗遊走，實質為對稱 t。**修正前後的 skew-t + TS 結果不可直接比較**；
非 skew 或非 TS 的設定不受影響。另 `US_ALL_RESULTS_DIR` 過去未被 runner 讀取（寫死
`results/`），現已生效——舊紀錄若聲稱寫到別的目錄，實際都在 `results/`。

---
