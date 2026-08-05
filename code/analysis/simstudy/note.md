# z.init bug 修復前實驗掃描結果（提醒：日後是否需要重跑）

掃描日期：2026-07-17。
Bug 全文見 `tex/z_init_bug/z_init_bug.tex`（time-block 研究的 `z.init = 0.01` 事故）
與 `tex/root_cause_fix/z_init_fix.tex`（backend 預設值的根因修復）。

## 結論（TL;DR）

**本目錄（`code/analysis/simstudy`）目前磁碟上沒有任何需要因 z.init bug 重跑的實驗。**

- 所有受影響方法（m7–m10，時間性 z update）的現存 fit 檔案 mtime 均 ≥ 修復日。
- 唯一落在曖昧時間窗的 set9 m7/m8（20 個檔，2026-05-28 凌晨）已逐檔驗證健康（見下）。
- 其餘結果（`results_def`、`results_nonsta`、scores3/4、`old/`、max-stab、lm）
  全部只用非時間性方法，結構上不可能踩到這個 bug。

## Bug 影響面（哪些東西才可能中招）

z.init 初始化 bug 只影響 **copula 空間的時間性 z update**：

| 進入點 | 檔案 | 對應方法 |
|---|---|---|
| `updateZTS`（AR(1)） | `code/R/ar2/update_params.R:41`（`hn.cop`） | method 9、10 |
| `updateZTS_AR2` | `code/R/ar2/update_params_ar2.R:352`（`hn.cop`） | method 7、8 |

非時間性的 z update 不走 `hn.cop` 起始路徑，因此 **methods 1–5、11、12（含所有
MRTS 變體）、max-stable（m6）、lm baseline 一律不受影響**，無論何時跑的。

## 修復時間線（判定基準）

| 時間 | commit | `z.init` 預設 | 後果 |
|---|---|---|---|
| ~2026-05-27 之前 | （735a566 起） | `0` | **鏈凍結**：`hn.cop(0,σ) = −∞` → MH ratio NaN → φ_z ≡ 0、z 全程不動 |
| 2026-05-27 17:43 | aabe3f4 | `1` | 過渡修復：z* ≈ 0.47（τ=1 時），可用但偏離中心 |
| 2026-05-28 09:57 | f7770f6 | `NULL` → `0.6745/√τ`（HN 中位數） | 正式修復：z* = 0，對任何 τ 都在 copula 尺度正中心 |

本目錄的 `run-settings.R` **從不顯式傳 `z.init`**（time-block 研究才有 `z.init = 0.01`
的災難），所以判定只看「fit 是在哪個預設值下跑的」＝看檔案時間。

## 掃描明細

### `results/`（主研究，settings 9–19）

| 批次 | 檔案日期 | 判定 |
|---|---|---|
| set9 m7/m8（20 檔） | 2026-05-28 00:17–03:36（修復 commit 前的凌晨） | **曖昧窗，已驗證健康**（見下節） |
| set10–12 m7/m8 | 2026-05-28 11:57 之後 | 修復後 ✓ |
| set13–15 m7/m8 | 2026-05-29 | 修復後 ✓ |
| set9–15 m9/m10 | 2026-05-31 | 修復後 ✓ |
| set16 m7–m10 | 2026-07-03 | 修復後 ✓ |
| set17–19 m7/m9 | 2026-07-12 | 修復後 ✓ |
| 全部 m1–m5 及 MRTS | （任意日期） | 非時間性，無關 ✓ |

凍結 bug 時代（z.init=0）的 m7/m8 原始 fit 已全數被上述重跑覆蓋，磁碟上不存在。

### set9 m7/m8 曖昧窗的驗證（2026-07-17 執行）

set9 的 20 個檔案介於過渡修復（z.init=1）與正式修復之間，無法從時間戳判定用哪個
初始值——但兩者都不是災難值，只要鏈有收斂即無需重跑。逐檔檢查
`tex/z_init_bug` 提出的兩條自我一致性判準，全數通過：

1. **φ_1z 未凍結**：posterior SD 介於 0.06–0.37（凍結時恆為 0）。
2. **z 自我一致**：z̄ ÷ (√(2/π)·mean(1/√τ)) 介於 0.81–1.26（壞掉的 run 約 0.14，差 7 倍）。
3. **λ 量級正常**：多數 ≈ +3；壞掉的 run 是 −50 ~ −100。

與確定修復後的對照檔（13-7-1、13-8-1、9-9-1、16-7-1）行為一致。**不需重跑。**

⚠ 附帶觀察（與本 bug 無關）：9-7-1 / 9-7-8 / 9-7-9 的 λ ≈ −4.1 ~ −4.4（z̄ 健康），
屬已知的 (β₀, λ, z) ridge 不可辨識問題（`z_init_bug.tex` §Residual caveat），另案處理。

### 其他目錄

| 位置 | 內容 | 判定 |
|---|---|---|
| `results_def/` + `output/results/scores*_def.RData`（2026-05-09/10，修復前） | 只有 m1–5 + MRTS | 非時間性，無關 ✓ |
| `results_nonsta/` + `scores*_nonsta*` | run plan 只含 m1/m2/m4/m5/m11/m12 + MRTS | 非時間性，無關 ✓ |
| `output/results/scores3、scores4、posterior3、posterior4`（2026-05-28） | setting 3/4 run plan 只有 m1–5 + MRTS | 非時間性，無關 ✓ |
| `old/`（scores1–8.RData，2025-10） | 舊世代研究，6 個方法，score 物件無任何 φ_z 欄位＝無 temporal 成分 | 無關（且已封存淘汰）✓ |
| `max-stab/`、lm baseline | 無 MCMC z update | 無關 ✓ |

### 範圍外備註：`code/analysis/time_block_forecast/simstudy`

該研究才是 bug 現場（顯式 `z.init = 0.01`）。已於修復時處理：壞結果封存於
`results_BUGGY_backup/`，`run-settings.R` 已移除該參數並加註警告，前後對照由
`compare_brier_prepost.R` 保留 —— 該腳本與其產出已於 2026-08-05 移至
`code/analysis/time_block_forecast/block1_positive_control/`（simstudy 只留正式
實驗，探索與診斷腳本歸 positive_control）。**該目錄修復前的所有 score（CRPS/Brier/
energy/variogram）一律作廢，不得引用。**

## 日後判準（若挖出舊備份或懷疑某個 temporal fit）

任何 `temporalz`/`ar2_z` 的 fit 滿足以下任一條即丟棄重跑：

- 檔案 mtime < 2026-05-28；或
- `sd(fit$phi.z[,1]) == 0`（凍結）；或
- `mean(fit$z)` 遠小於 `sqrt(2/pi) * mean(1/sqrt(fit$tau))`（比值 ≪ 1，如 ~0.14）。

快速檢查片段：

```r
load("results/9-7-1.RData")   # -> fit.1
c(sd_phi1z = sd(fit.1$phi.z[, 1]),
  ratio    = mean(fit.1$z) / (sqrt(2/pi) * mean(1/sqrt(fit.1$tau))),
  lambda   = mean(fit.1$lambda))
# 健康：sd > 0、ratio ≈ 1、|lambda| ~ O(1)
```
