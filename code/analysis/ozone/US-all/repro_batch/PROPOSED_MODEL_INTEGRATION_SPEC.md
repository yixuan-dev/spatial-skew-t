# Proposed model 插入規格（US-all）

本規格目標：
1. 不破壞 Morris 原始重現流程（`done` 與 `us-all-results.R` 可原封不動重跑）。
2. 讓你能把 proposed model 以同一資料、同一 CV、同一評分函數接進比較。

---

## A. 檔案新增清單（建議最小集合）

### 1) 模型執行腳本（每個設定一支）
- 命名：`us-all-<setting>.R`
- 位置：`code/analysis/ozone/US-all/`
- 建議：先做 3 個 proposed 設定（例如不同 K 或不同閾值）做 smoke test，再擴增。

### 2) proposed 設定表
- 新增：`settings_proposed.csv`
- 欄位沿用 `settings.csv`：
  - `setting, method, knots, thresh, CMAQ, TS, rerun, running`
  - 可再加：`model_tag, note`

### 3) proposed 匯總腳本（不要直接覆寫原版）
- 新增：`us-all-results-proposed.R`
- 功能：
  - 讀取 `us-all-setup.RData`
  - 讀取 baseline + proposed 的 `results/us-all-*.RData`
  - 計算 `QuantScore/BrierScore`
  - 產生與 `us-all-results.R` 同格式指標（至少 `brier.score.mean`、相對 Gaussian）
  - 輸出 `us-all-results-proposed.RData` 與比較表（CSV）

### 4)（可選）共用函數檔
- 新增：`proposed_model_api.R`
- 放 proposed 模型的主要 fitting/predict 介面，避免每支 `us-all-<setting>.R` 重複貼程式。

---

## B. setting 編碼規則（避免撞號）

目前 Morris 已使用 `1..74`。

### 建議規則
- **保留 Morris 不動**：`1..74`
- **proposed 從 101 起編號**（或 201 起）：
  - 例：`101, 102, 103, ...`
- 規則：
  - 不和 `1..74` 重疊
  - 每個 setting 對應唯一腳本 `us-all-<setting>.R`
  - `settings_proposed.csv` 必有同一 setting 的 metadata

### 範例（僅示意）
- 101: proposed, K=5, T=50, no TS
- 102: proposed, K=7, T=50, no TS
- 103: proposed, K=7, T=75, TS

---

## C. `us-all-<setting>.R` 契約（必須符合）

為了能被 `us-all-results*.R` 讀取，你的腳本需輸出：
- 物件名稱：`fit`
- 型別：`list` 長度 2（對應 2-fold CV）
- 每個 fold 至少要有：
  - `fit[[d]]$yp`（MCMC 後驗預測樣本）

`$yp` 維度需相容既有評分函數：
- `pred.d <- fit.d$yp[, , ]`
- 可被 `QuantScore(pred.d, probs, validate)`
- 可被 `BrierScore(pred.d, thresholds, validate)`

> 重點：只要 `fit[[d]]$yp` 契約一致，`us-all-results` 的比較邏輯就能沿用。

---

## D. `us-all-results.R` 接入策略（建議）

### 不建議
- 直接大改原 `us-all-results.R`（會破壞純 Morris 重現）。

### 建議
- 新增 `us-all-results-proposed.R`，分兩層：
  1. **baseline 層**：讀 Morris full-done（1..74）。
  2. **proposed 層**：讀你新增 setting（>=101）。

### 實作重點
1. 把陣列第三維從固定 74 改成動態 `max_setting` 或字典式儲存。
2. `done_all <- c(done_morris, done_proposed)`。
3. 相對分數仍用 Gaussian baseline（setting 1）：
   - `rel_bs[i, ] <- brier.score.mean[i, ] / brier.score.mean[1, ]`
4. 最終輸出：
   - `comparison_top2.csv`（每個 quantile 的 top1/top2，含 setting、model_tag、K/T/TS、rel BS）

---

## E. 公平比較規範（強烈建議固定）

1. 同一 `us-all-setup.RData`（同資料、同站點）
2. 同一 `cv.lst`（同 fold）
3. 同一 `probs/thresholds`
4. 同一評分函數 `QuantScore/BrierScore`
5. 同級 MCMC 預算（或清楚註記差異）

---

## F. 最小 smoke test 流程

1. 先新增 1 個 proposed setting（如 101）
2. 跑 `us-all-101.R` 產生 `results/us-all-101.RData`
3. 用小版 `done_all = c(1, 3, 38, 101)` 測 `us-all-results-proposed.R`
4. 確認 `comparison_top2.csv` 產生且欄位齊全
5. 再擴增 proposed settings

---

## G. 交付輸出（最終）

- `us-all-results-proposed.RData`
- `comparison_top2.csv`
- `comparison_full_table.csv`（所有模型 × quantile 的 rel BS）
- 寫入論文表格時，主張「排名與結論一致」，並補充數值可能有平台微差。
