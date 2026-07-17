# 非平穩模擬：MCMC 結果與未解問題

> 記錄日期：2026-07（pilot 階段）
> 設計文件：[nonsta_settings_design.md](nonsta_settings_design.md)（講「怎麼建資料」）
> 本文件講「跑出什麼、什麼還沒解」。
>
> ⚠ **所有 MCMC 結果都是 pilot，datasets 1–3（n = 3）**，除了 setting 1 是既有快取的
> n = 10。異常在 n = 3 下已經很大且跨 method/dataset 一致，但**寫進論文前需補到 n = 10**
> （每個 setting 135 fits，約 3.3 h）。

---

## 已跑的 setting

| setting | surf_type | n (datasets) | fits | 狀態 |
|---------|-----------|--------------|------|------|
| 1 | invariant | 10（既有快取） | — | scores1_nonsta.RData |
| 6 | covdep_cov | 3 | 135 | scores6_nonsta.RData，fits 保留在 results_nonsta/ |
| 11 | invariant_strong | 3 | 135 | scores11_nonsta.RData，fits 保留 |

網格一律：methods 1–5 × K ∈ {0,3,5,8,10,12,15,20,25}。
MRTS 基底來源全部驗證為 `autoFRK::mrts`（訓練站點），**無 combined fallback 洩漏**。

---

## 結果 1（穩固）：recovery 的 headline 對照，setting 1 vs 6

同一個 f1(s)，setting 1 放**均值**、setting 6 放**相關 range** rho(s)=rho0·exp(0.7·f1)。

主力 skew-t（method 2, K=1）的 relative recovery RMSE（相對於各 method 自己 K=0）：

| K_MRTS | 5 | 10 | 15 | 20 | 25 |
|--------|-----|-----|-----|-----|-----|
| **Setting 1**（f1 在均值） | 0.18 | 0.15 | 0.072 | 0.057 | **0.054** |
| **Setting 6**（同一 f1 在 range） | 1.02 | 1.02 | 1.05 | 1.02 | **1.02** |

**結論成立**：同一個空間 pattern，換動差層級，MRTS 從恢復 94.6% 變成完全不動。
排除「資訊不足」的替代解釋。圖：`output/plots/headline_nonsta_1v6.pdf`。

> ⚠ recovery.rmse 在「真均值是平的」setting（4/6/8）上會**雙峰**——skew-t 的截距與 λ 共線
> （μ = X'β + λz_t，z_t 半常態），詳見 [nonsta_settings_design.md](nonsta_settings_design.md) 該節。
> setting 1 因為空間訊號強（sd 3.68）不受影響；上面 0.054 可信。

---

## 結果 2（穩固，論文級）：oracle Brier 天花板——9/11 setting 理論上不可能動 Brier

腳本：[oracle_brier_ceiling.R](oracle_brier_ceiling.R) → `output/oracle_brier_ceiling.csv`

因為 MRTS 的 β 跨時間 pooled，均值通道極限 = 時間平均均值面 μ̄。定義下限

  R* = BS(μ̄) / BS(空間常數)，  在真 skew-t 預測律下算。

R* 是**下限**，無 MCMC 可算。結果：

| id | surf_type | 動差 | sd(μ̄) | **R*** | 判定 |
|----|-----------|------|--------|--------|------|
| 1 | invariant | mean | 3.68 | 0.856 | 勉強（上限僅 14%） |
| 2 | varying | mean | 0.42 | 1.000 | 死路 |
| 3 | ns_dependence | dependence | 0.41 | 0.999 | 死路 |
| 4 | ps_cov | dependence | 0.00 | 1.000 | 死路 |
| 5 | ps_add | dependence | 0.39 | 0.999 | 死路 |
| 6 | covdep_cov | dependence | 0.00 | 1.000 | 死路 |
| 7 | covdep_add | dependence | 0.41 | 1.000 | 死路 |
| 8 | fuentes_cov | dependence | 0.00 | 1.000 | 死路 |
| 9 | gh_marginal | tails | 0.42 | 1.000 | 死路 |
| 10 | lambda_varying | skewness | 2.80 | **1.120** | **反效果** |
| 11 | invariant_strong | mean | 11.04 | **0.375** | positive control |

三個推論（見設計文件詳述）：
1. **setting 1 從設計上就當不了 Brier positive control**（天花板僅 14%，實測 1.5%）。
2. **settings 2–9 的 Brier 虛無不具資訊**——是設計性質，非關於 MRTS 的證據。
   特別是 replace_C（4/6/8）：邊際已驗證不變（KS p 0.075–0.864），而 Brier 是逐格邊際分數。
3. **setting 10 的 R* > 1**：餵真均值反而讓 Brier 變差（λ(s) 的均值位移只是偏態病灶的症狀）。

**Brier 不動 ≠ MRTS 沒用**：同一批 fit，method 3 的 energy score 改善 42%（0.630→0.579），
Brier 只動 0.9%。多變量分數（energy/vario）才有鑑別力。

---

## ✅ 已定案（2026-07-17）：兩個判決

### 判決 1：threshold vs symmetry —— **門檻是驅動因素**

Methods 11/12（Sym-t **無門檻**，與 3/5 唯一差異是門檻）在 setting 1（n=10）擬合完成
（fits 保留在 `results_nonsta/1-11-*` / `1-12-*`；分數在 `scores1_nonsta_m1112.RData`，
原 `scores1_nonsta.RData` 快取未被覆寫）。

relative energy（vs own K=0）：m3（有門檻）0.58–0.63 大降；**m11（無門檻）0.99 完全持平**。
m5 vs m12 同型。→ 增益是門檻專屬。

**但絕對值反轉了解讀**：K=0 的絕對 energy —— m2/m11 = 8.5、**m3 = 21.8、m5 = 137**。
POT censoring 在沒有均值基底時把預測場毀掉（baseline 差 2.6×–16×）；MRTS 是**部分救援**
（21.8×0.579≈12.6，仍輸給從不需要救的 8.5），不是綜效。正確結論：
「只有門檻 method 有被救的空間；絕對值上無門檻 method 在所有 K 都贏」。

### 判決 2：setting 11 的 Brier 反向異常 —— **在大 K 解決**

Setting 12（6-bump, sd 11）跑到 K=50（n=3）。Brier 相對比：K≤25 平-到-差（1.00–1.05），
**K=30 起轉好、K=50 改善 7–12%**（Gau 0.930、StK1 0.925、SyK1-T 0.880、StK5 0.913；
m5 被 fallback 爆值毀掉一個 dataset 不計）。與外推機制的預測完全一致：投影誤差
0.13@K25 → 0.057@K50，均值殘差小到讓銳利預測分布得利才翻正。

recovery 也實現了投影曲線的延伸手肘：StK1 0.93@K5 → 0.49@K10 → 0.14@K25 → **0.073@K50**。
K50 的 7–12% 仍遠高於 oracle 下限 R\*=0.495——大部分理論空間被 sharpness-extrapolation
trade-off 吃掉。「Brier 對均值恢復是壞指標」的結論保持，只是機制現在有實證支持。

報告 v2：`tex/nonsta_simstudy/nonsta_simstudy.tex`（7 頁，含 disambiguation 表與 setting 12 表）。

---

## ~~⚠⚠ 未解問題~~（已解決，見上）：setting 11 的 positive control 反向失敗

**這是目前最重要、還沒解釋完的異常。下次接手從這裡開始。**

### 現象

Setting 11（= setting 1 的均值面 ×3，sd(μ̄)=11.04）。oracle 說 Brier 最多可改善 62.5%（R*=0.375）。
**實測卻是均值結構越強、MRTS 讓 Brier 越差。** relative Brier（vs own K=0，n=3）：

| method | K=5 | K=10 | K=12 | K=25 |
|--------|-----|------|------|------|
| 1 Gauss | 1.059* | 1.101* | 1.114* | 1.051* |
| 2 SkewT-K1 | 1.052* | 1.103* | **1.163*** | 1.086* |
| 4 SkewT-K5 | 1.088* | 1.123* | **1.184*** | 1.112 |
| 5 SymT-K5 | 1.078 | 1.114* | 1.163* | 1.133 |
| 3 SymT-K1 | 1.016 | 1.039* | 1.008 | 0.959* |

（`*` = |R−1| > 2 SE。K=3 時 ≈ 1，一加 TPS 特徵函數就開始惡化。）

### 為什麼這是矛盾

同一批 fit：
- **recovery.rmse 大幅改善**：method 2 從 14.010（K=0）→ **0.659**（K=25），降 95%。
- **所有參數在 K=25 恢復到真值**（method 2, dataset 1；真值 λ=3, E[σ]=2.257, ρ=1, γ=0.9）：

  | K | λ̂ | E[σ̂] | ρ̂ | γ̂ |
  |---|-----|-------|-----|-----|
  | 0 | 0.15 | 11.16 | 2.04 | 0.99 |
  | 10 | 2.04 | 3.01 | 1.21 | 0.91 |
  | 25 | **3.18** | **1.94** | **0.96** | **0.94** |

  K=0 是嚴重誤設定：配不出 sd=11 的均值面，把空間變異塞進膨脹的 σ（真值的 5 倍）、λ 壓到 0.15。
  K=25 把每個參數都修對。

**所以：一個在每個參數上都正確設定、均值恢復 95% 的模型，Brier 竟然輸給嚴重誤設定的 K=0 模型 8–16%。**

### 候選解釋（都還沒驗證）

1. **樣本外預測分布的偏誤放大**：K=0 的膨脹 σ 意外地給了保留站點「較寬、較保險」的預測區間；
   K=25 的窄 σ 雖然對訓練站點正確，但 MRTS 基底在**測試站點**外推的均值有誤差，窄 σ 把這個
   均值誤差暴露成尾端機率的大偏差 → Brier 罰。即「正確的 σ + 有偏的外推均值」比「膨脹的 σ +
   常數均值」在**樣本外 Brier** 上更糟。
2. **Brier 在高分位對均值-尺度權衡的非線性**：需要看 relative Brier **逐 exceedance 分位**
   （p=0.90…0.995）在哪裡崩掉。← **這個診斷指令上次被分類器中斷，是下一步第一件事。**
3. **n=3 的偶然**：可能性低（跨 4 個 method、跨 K 單調、SE 小），但補到 n=10 才能排除。

### 下一步（接手指令）

```r
# 1. relative Brier by exceedance probability, setting 11 method 2 —— 找崩在哪個分位
#    (上次 /tmp/why5.R 被中斷，重跑)
# 2. 分解 in-sample vs out-of-sample：MRTS 均值在 100 訓練站點 vs 44 測試站點的 RMSE 差多少
# 3. 若機制 1 成立 —— 這其實是另一個論文級發現：
#    「mean-only augmentation 即使在有均值結構、且能完美恢復參數時，
#     樣本外機率預測仍可能變差，因為它移除了誤設定 σ 的無意保險。」
# 4. 確認機制後，補 setting 1 + 11 到 n=10（各 135 fits, ~3.3h）
```

### 已驗證不是原因

- 不是 MRTS 洩漏：basis source 全部 autoFRK::mrts（訓練站點），無 fallback。
- 不是資料生成錯：setting 11 = setting 1 ×3，cor(g1,g11)=1.0000，settings 1–10 位元不變。
- 不是均值沒恢復：recovery 降 95%，參數全對。

---

## 檔案索引

| 檔案 | 內容 |
|------|------|
| `output/results/scores{1,6,11}_nonsta.RData` | scored 結果 |
| `results_nonsta/{6,11}-*.RData` | 保留的 raw fits（各 135，約 37 GB/setting） |
| `output/oracle_brier_ceiling.csv` | oracle R* 表 |
| `output/plots/headline_nonsta_1v6.pdf` | recovery headline 圖 |
| `output/nonsta_diagnostics.{pdf,csv}` | DGP 非平穩性診斷 |
| [oracle_brier_ceiling.R](oracle_brier_ceiling.R) | oracle 計算（論文級） |
| [plot_headline_1v6.R](plot_headline_1v6.R) | headline 對照圖 |
