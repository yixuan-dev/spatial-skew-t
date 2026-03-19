# F1 超精簡 Runbook（full-done reproduction）

目標：完成 F1 批次並通過檢查，再進入 F2。

## F1 範圍
- settings: `1,2,3,4,5,7,8,9,11,12,13,15,16,17`

## 執行前（一次）
- 工作目錄：`d:\Github\spatial-skew-t\code\analysis\ozone\US-all`
- 確認：`us-all-setup.RData` 已存在

## 逐項執行（順序）
1. 跑 `us-all-1.R`
2. 檢查 `results/us-all-1.RData` 存在
3. 跑 `us-all-2.R`
4. 檢查 `results/us-all-2.RData` 存在
5. 依序重複上述模式直到 `us-all-17.R`（跳過 6,10,14）

## F1 完成條件
以下檔案全部存在：
- `results/us-all-1.RData`
- `results/us-all-2.RData`
- `results/us-all-3.RData`
- `results/us-all-4.RData`
- `results/us-all-5.RData`
- `results/us-all-7.RData`
- `results/us-all-8.RData`
- `results/us-all-9.RData`
- `results/us-all-11.RData`
- `results/us-all-12.RData`
- `results/us-all-13.RData`
- `results/us-all-15.RData`
- `results/us-all-16.RData`
- `results/us-all-17.RData`

## 失敗處理（最短路徑）
- 若某一 setting 失敗：
  - 先記錄該 setting 編號
  - 先繼續下一個 setting
  - F1 跑完後只補跑失敗那幾個

## 進入下一批
- F1 完成條件滿足後，進入 F2（`33:36, 38:41, 43:46`）。
