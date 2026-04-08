# us-all-run.R PowerShell 使用指南

## 📌 快速開始

### 1️⃣ 定位至工作目錄

```powershell
Set-Location "D:\Github\spatial-skew-t\code\analysis\ozone\US-all"
```

**提醒：** 所有命令必須在此目錄執行。腳本會自動設置工作目錄，但建議先手動切換確保穩定性。

---

## 🎯 基本用法

### 執行單一或多個設定

```powershell
# 執行設定 114
Rscript us-all-run.R 114

# 執行多個設定
Rscript us-all-run.R 1 2 3 5 54

# 執行範圍（1 到 124）
Rscript us-all-run.R 1:124

# 執行多個範圍
Rscript us-all-run.R 1:20 50:60 100:114
```

if Rscript is in PATH, you can simply run:

```powershell
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" us-all-run.R 114
```

### 使用 PowerShell 7 平行腳本（`run-parallel-ps7.ps1`）

如果你要同時跑多個 settings，建議直接使用已儲存的腳本：`run-parallel-ps7.ps1`。

```powershell
Set-Location "D:\Github\spatial-skew-t\code\analysis\ozone\US-all"
.\run-parallel-ps7.ps1
```

腳本內容如下（與目前檔案一致）：

```powershell
Set-Location "D:\Github\spatial-skew-t\code\analysis\ozone\US-all"

$settings = (1..74) | Where-Object { $_ -ne 2 }
$throttle = 4               # 同時最多幾個進程（建議先 2~4）
$backend = "ar2"        # legacy / ar2
$runMode = "prod"          # dev / prod
$results = "results_new"   # 預設輸出資料夾
$wd = "D:\Github\spatial-skew-t\code\analysis\ozone\US-all"
$rscript = (Get-Command Rscript -ErrorAction Stop).Source

$settings | ForEach-Object -Parallel {
    $s = $_
    Set-Location $using:wd
    $env:US_ALL_MCMC_BACKEND = $using:backend
    $env:US_ALL_RUN_MODE = $using:runMode
    $env:US_ALL_RESULTS_DIR = $using:results

    & $using:rscript "us-all-run.R" $s
} -ThrottleLimit $throttle
```

> 若終端機不是 PowerShell 7，`ForEach-Object -Parallel` 會不可用。請先確認版本，或改用 `Start-Job` 版本。

---

## 🔧 環境變數配置

四個核心環境變數控制腳本行為：

### 1. `US_ALL_MCMC_BACKEND` — 選擇 MCMC 函數

| 值       | 函數         | 載入檔案         | 說明               |
| -------- | ------------ | ---------------- | ------------------ |
| `legacy` | `mcmc()`     | `package_load.R` | **預設**，原始實驗 |
| `ar2`    | `mcmc_ar2()` | `ar2_load.R`     | 新版 AR(2) 模型    |

**範例：**
```powershell
# 使用 Legacy（預設）
Rscript us-all-run.R 114

# 使用 AR2 後端
$env:US_ALL_MCMC_BACKEND = "ar2"
Rscript us-all-run.R 114
```

**重要提醒：** 若要新增第三個後端（例如 `mcmc_ar3`），需要：
- 在 `settings.csv` 中新增對應行
- 建立載入檔案 `ar3_load.R`（定義 `mcmc_ar3()` 及依賴函數）
- 在 `us-all-run.R` 第 46–49 行修改 backend 邏輯：
  ```r
  if (backend == "ar2") {
    source("./ar2_load.R", chdir = TRUE)
  } else if (backend == "ar3") {
    source("./ar3_load.R", chdir = TRUE)
  } else {
    source("./package_load.R", chdir = TRUE)
  }
  ```

---

### 2. `US_ALL_RUN_MODE` — 控制迭代次數

| 值     | iters | burn  | update | 用途               |
| ------ | ----- | ----- | ------ | ------------------ |
| `dev`  | 2000  | 1000  | 200    | 快速檢查、診斷     |
| `prod` | 30000 | 25000 | 500    | **預設**，完整實驗 |

**範例：**
```powershell
# Dev 模式（快速，~1-2 分鐘/fold）
$env:US_ALL_RUN_MODE = "dev"
Rscript us-all-run.R 114

# Prod 模式（完整，~30+ 分鐘/fold）
$env:US_ALL_RUN_MODE = "prod"
Rscript us-all-run.R 1:124
```

**輸出範例：**
```
RUN_MODE: dev | iters: 2000 | burn: 1000 | update: 200
```

---

### 3. `US_ALL_RESULTS_DIR` — 指定結果存放目錄

**預設值：** `results_new`（相對於當前目錄）

| 設定         | 存放位置                                                          |
| ------------ | ----------------------------------------------------------------- |
| 不設定       | `D:\Github\spatial-skew-t\code\analysis\ozone\US-all\results_new` |
| 設定絕對路徑 | 自訂位置                                                          |

**範例：**
```powershell
# 使用預設目錄
Rscript us-all-run.R 114
# → 結果存在：results_new/us-all-114.RData

# 自訂目錄
$env:US_ALL_RESULTS_DIR = "D:\scratch\exp_2026_03\run1"
Rscript us-all-run.R 1:124
# → 結果存在：D:\scratch\exp_2026_03\run1\us-all-1.RData 等

# 相對路径
$env:US_ALL_RESULTS_DIR = "..\results_archive"
Rscript us-all-run.R 50:60
```

**特性：**
- ✅ 目錄不存在時自動建立（含上層目錄）
- ✅ 會列印：`Created results directory: ...`
- ✅ 跨平台相容（Windows/Linux/Mac）

---

### 4. `US_ALL_SETTINGS` — 環境變數指定設定

**推薦場景：** 批次投稿或 HPC 作業系統

```powershell
# 用環境變數指定（逗號分隔）
$env:US_ALL_SETTINGS = "1,2,3,5,54"
Rscript us-all-run.R

# 或用範圍表達式
$env:US_ALL_SETTINGS = "1:20,50:60,100:124"
Rscript us-all-run.R
```

**優先順序：**
1. `US_ALL_SETTINGS` 環境變數（若已設定）
2. 命令列引數（若無環境變數）
3. 報錯：需要提供設定

---

## 📋 完整使用範例

### 場景1：快速測試 AR2 後端

```powershell
Set-Location "D:\Github\spatial-skew-t\code\analysis\ozone\US-all"

$env:US_ALL_MCMC_BACKEND = "ar2"
$env:US_ALL_RUN_MODE = "dev"
$env:US_ALL_RESULTS_DIR = "temp_test"

Rscript us-all-run.R 114 115 116
```

**預期輸出：**
```
Created results directory: temp_test

==============================
Setting: 114 | Backend: ar2
RUN_MODE: dev | iters: 2000 | burn: 1000 | update: 200
...
CV 1 finished. Fold sec: 45.23 | Avg sec/dataset: 45.23
CV 2 finished. Fold sec: 47.89 | Avg sec/dataset: 46.56
Saved: temp_test/us-all-114.RData
```

### 場景2：完整投稿實驗

```powershell
Set-Location "D:\Github\spatial-skew-t\code\analysis\ozone\US-all"

# 清除舊環境變數
Remove-Item Env:\US_ALL_MCMC_BACKEND -ErrorAction SilentlyContinue
Remove-Item Env:\US_ALL_RUN_MODE -ErrorAction SilentlyContinue

# 使用預設：legacy backend, prod mode, results_new 目錄
Rscript us-all-run.R 1:124

# 結果保存在 results_new/ 子資料夾中
```

### 場景3：比較兩個後端

```powershell
Set-Location "D:\Github\spatial-skew-t\code\analysis\ozone\US-all"

# Legacy 執行
$env:US_ALL_MCMC_BACKEND = "legacy"
$env:US_ALL_RESULTS_DIR = "results_legacy"
Rscript us-all-run.R 101 102 103

# AR2 執行
$env:US_ALL_MCMC_BACKEND = "ar2"
$env:US_ALL_RESULTS_DIR = "results_ar2"
Rscript us-all-run.R 101 102 103

# 比較結果
# legacy 版本：results_legacy/us-all-101.RData
# ar2 版本：results_ar2/us-all-101.RData
```

---

## 📝 環境變數管理

### 查看當前設定

```powershell
# 查看特定變數
$env:US_ALL_MCMC_BACKEND
$env:US_ALL_RUN_MODE
$env:US_ALL_RESULTS_DIR

# 或一次查看所有 US_ALL_* 變數
(Get-ChildItem Env:US_ALL_*).foreach({'$($_.Name) = $($_.Value)'} | Format-Table -AutoSize
```

### 清除變數（恢復預設）

```powershell
# 清除單個變數
Remove-Item Env:\US_ALL_MCMC_BACKEND -ErrorAction SilentlyContinue

# 清除所有 US_ALL_* 變數
(Get-ChildItem Env:US_ALL_*).ForEach({ Remove-Item Env:\$($_.Name) })
```

### 在 PowerShell 設定檔中永久設定

編輯 `$PROFILE`（或新建）：

```powershell
# 列印當前 profile 路徑
Write-Host $PROFILE

# 編輯（用 Notepad 或你慣用的編輯器）
notepad $PROFILE
```

新增至 profile：

```powershell
# US-all 專案的預設設定
$env:US_ALL_MCMC_BACKEND = "legacy"
$env:US_ALL_RUN_MODE = "prod"
$env:US_ALL_RESULTS_DIR = "results_new"
```

---

## 🚀 進階擴展

### 新增 MCMC 後端的步驟

假設要新增 `mcmc_ar3()` 後端：

#### Step 1：建立載入檔案 `ar3_load.R`

在同目錄下建立 `ar3_load.R`：

```r
# ar3_load.R
# 載入 AR3 相關函數及依賴

library(your_package)

# 定義或匯入 mcmc_ar3 函數
mcmc_ar3 <- function(...) {
  # AR3 邏輯
}

# 其他必要的依賴函數
source("./other_dependencies.R", chdir = TRUE)
```

#### Step 2：修改 `us-all-run.R` 的後端判斷

在第 46–49 行修改：

```r
if (backend == "ar2") {
  source("./ar2_load.R", chdir = TRUE)
} else if (backend == "ar3") {  # 新增
  source("./ar3_load.R", chdir = TRUE)
} else {
  source("./package_load.R", chdir = TRUE)
}
```

同時更新驗證邏輯（第 44–45 行）：

```r
backend <- tolower(Sys.getenv("US_ALL_MCMC_BACKEND", unset = "legacy"))
if (!backend %in% c("legacy", "ar2", "ar3")) {  # 加入 "ar3"
  stop("US_ALL_MCMC_BACKEND must be one of: legacy, ar2, ar3", call. = FALSE)
}
```

#### Step 3：更新 `settings.csv`

確保新設定行有 `ar3` 欄位設置。

#### Step 4：測試

```powershell
$env:US_ALL_MCMC_BACKEND = "ar3"
$env:US_ALL_RUN_MODE = "dev"
Rscript us-all-run.R 101
```

---

## ⚠️ 常見問題

### Q1：如何查看實際執行的參數？

每個 setting 開始時會列印完整配置：

```
==============================
Setting: 114 | Backend: ar2
RUN_MODE: dev | iters: 2000 | burn: 1000 | update: 200
Method: skew-t | K: 20 | Threshold: 0.02
CMAQ: yes | TS: yes | AR2 row: yes
```

### Q2：結果文件去哪了？

檢查 `US_ALL_RESULTS_DIR` 設定：

```powershell
# 預設位置
Get-ChildItem results_new/ | Where-Object Name -like "us-all-*.RData"

# 或自訂位置
Get-ChildItem $env:US_ALL_RESULTS_DIR
```

### Q3：如何在現有 AR2 基礎上新增 AR3？

見上文 **進階擴展** 章節。簡言之：
1. 建立 `ar3_load.R`
2. 修改 us-all-run.R 的 backend 判斷
3. 測試

### Q4：範圍表達式有限制嗎？

- ✅ `1:124` 支持
- ✅ `1:30,50:60` 同時多個範圍
- ✅ `1` 單一數字混用
- ❌ 不支持倒序（`124:1` 作為 124 到 1）、浮點數

---

## 📊 計時和監控

每個 fold 完成時會列印：

```
CV 1 finished. Fold sec: 123.45 | Avg sec/dataset: 123.45
CV 2 finished. Fold sec: 145.67 | Avg sec/dataset: 134.56
```

- **Fold sec**：該 fold 花費的秒數
- **Avg sec/dataset**：到目前為止的平均（用於估算剩餘時間）

---

## 📁 檔案結構

```
D:\Github\spatial-skew-t\code\analysis\ozone\US-all\
├── us-all-run.R              ← 統一執行入口（本討論的主角）
├── run-parallel-ps7.ps1      ← PowerShell 7 平行執行腳本
├── package_load.R            ← Legacy 後端依賴
├── ar2_load.R                ← AR2 後端依賴
├── ar3_load.R                ← （未來）AR3 後端依賴
├── settings.csv              ← 124 個設定的定義
├── us-all-setup.RData        ← 預先計算的 Y, X, S, cv.lst 等
├── results_new/              ← 預設輸出目錄
│   ├── us-all-1.RData
│   ├── us-all-2.RData
│   └── ...
└── US-ALL-RUN-GUIDE.md       ← 本檔案
```

---

## 📞 總結

| 任務                 | 命令                                                             |
| -------------------- | ---------------------------------------------------------------- |
| 快速測試 setting 114 | `$env:US_ALL_RUN_MODE="dev"; Rscript us-all-run.R 114`           |
| 完整實驗 1–124       | `Rscript us-all-run.R 1:124`                                     |
| 用 AR2 執行          | `$env:US_ALL_MCMC_BACKEND="ar2"; Rscript us-all-run.R 114`       |
| 自訂結果位置         | `$env:US_ALL_RESULTS_DIR="C:\my\path"; Rscript us-all-run.R 114` |
| 新增後端             | 新建 `newbackend_load.R` + 修改 us-all-run.R 的 backend 判斷邏輯 |
