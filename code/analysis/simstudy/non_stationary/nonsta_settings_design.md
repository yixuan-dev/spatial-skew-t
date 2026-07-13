# `simdata_nonsta.RData` 實驗設計文件

> 對應實作：[setup_nonsta.R](../setup_nonsta.R)、[ns_cov_builders.R](ns_cov_builders.R)
> 診斷腳本：[diagnose_nonsta.R](diagnose_nonsta.R) → `output/nonsta_diagnostics.pdf` / `.csv`
> 姊妹資料集：[def_settings_design.md](def_settings_design.md)（變形路線）

---

## 設計目的

論文要主張的是「**mean-only 的 MRTS covariate 擴充有結構性的能力上限**」。舊版只有 3 個
setting，唯一的反例（setting 3）是低秩、振盪型的相依——一個特例。審稿人合理的反問是：
「你只證明了 MRTS 抓不到*這一種*非平穩，換成別種呢？」

這份資料集把「非平穩」展開成 **機制 × 動差層級** 的網格，讓結論從單點變成一條線：

| 非平穩所在 | Settings | MRTS 能不能撿到 |
|-----------|----------|----------------|
| **第一動差**（均值） | 1, 2 | **能** — 正控制組 |
| **第二動差**（相依） | 3–8 | **不能** — 而且換四種不同機制都一樣 |
| **高階動差**（偏態／尾重） | 9, 10 | 看有沒有誘導出均值 |

### 核心對照：Setting 1 vs Setting 6

兩者由**同一個** $f_1(s) = \cos(\pi\|s/10 - c_1\|)$ 驅動：

- **Setting 1** 把 $f_1$ 放進**均值**：$g(s) = 5f_1 + 3f_2$ → MRTS 完全恢復
- **Setting 6** 把 $f_1$ 放進**相關 range**：$\rho(s) = \rho_0 e^{0.7 f_1(s)}$ → MRTS 不論 $K$ 多大都無效

$f_1$ 是位置的平滑函數，MRTS 基底張得出來（setting 1 已經證明了這件事），模型也拿得到座標。
所以 setting 6 的失敗**不能**用「資訊不足」解釋——唯一的差別是真值把 $f_1$ 放在第二動差，
而 MRTS 只能把它放進均值。**動差層級錯了，不是資訊不夠。**

---

## 共同參數（全 10 個 settings 共用）

| 參數 | 值 | 說明 |
|------|-----|------|
| `beta.t` | `c(10, 0, 0)` | intercept-only；座標**不帶**線性趨勢 |
| `nu.t` | `0.5` | 指數族（Matérn 退化） |
| `gamma.t` | `0.9` | 邊際相關縮放（= 1 − nugget） |
| `rho.t` | `1` | 基準 range |
| `tau.alpha.t` / `tau.beta.t` | `3` / `8` | $\tau \sim \text{Gamma}(3/2, 4)$ → df 3 |
| `lambda` | `3` | 偏態參數（setting 10 讓它隨空間變動） |
| `dist` / `nknots` | `"t"` / `1` | skew-t-1 |
| `phi.z = phi.w = phi.tau` | `0` | 時間 i.i.d. |
| **站點** | 144 點，`runif([0,10]²)`，`set.seed(20)` | 與 `simdata.RData` 同網格 |
| `nt` / `nsets` / `ntest` | `50` / `50` / `44` | |
| `sigma.add` | `3` | 加法式隨機效應的振幅（settings 5/7/9） |

Cosine bump 基底（Tzeng & Huang 2018, Scenario 1；於 $s_{01} = s/10$ 上求值）：

$$f_1(s) = \cos\big(\pi \|s_{01} - c_1\|\big), \qquad f_2(s) = \cos\big(2\pi \|s_{01} - c_2\|\big)$$

$c_1 = (0, 1)$、$c_2 = (0.75, 0.25)$。這組基底**貫穿整份設計**：settings 1/2 用它當均值面、
setting 6 用 $f_1$ 當 range、setting 9 用 $f_1$ 當偏態、$f_2$ 當尾重、setting 10 用 $f_2$ 當 $\lambda$。

---

## Setting 總表

| ID | `surf_type` | 機制 | 注入路徑 | 動差 | 出處 |
|----|-------------|------|---------|------|------|
| 1 | `invariant` | 固定 cosine-bump 均值面 | 均值 | 1st | Tzeng & Huang (2018) |
| 2 | `varying` | 時變 cosine-bump 均值面 | 均值 | 1st | Tzeng & Huang (2018) |
| 3 | `ns_dependence` | 低秩 cosine 隨機效應 | 加法 | 2nd | Cressie & Johannesson (2008) |
| 4 | `ps_cov` | Paciorek–Schervish 局部各向異性 | **取代 C** | 2nd | Paciorek & Schervish (2006) |
| 5 | `ps_add` | 同 4，改當加法隨機效應 | 加法 | 2nd + 邊際 | 同上 |
| 6 | `covdep_cov` | $\rho(s) = \rho_0 e^{\beta f_1(s)}$ ★ | **取代 C** | 2nd | Schmidt et al. (2011) |
| 7 | `covdep_add` | 同 6，改當加法隨機效應 | 加法 | 2nd + 邊際 | 同上 |
| 8 | `fuentes_cov` | 4 區域加權混合 GP | **取代 C** | 2nd | Fuentes (2001, 2002) |
| 9 | `gh_marginal` | Tukey g-and-h 場 | 加法 | 3rd + 4th | Xu & Genton (2017) |
| 10 | `lambda_varying` | $\lambda(s) = 3 + 3f_2(s)$ ★ | 改 skew 項 | 3rd | Morris et al. (2017) 假設破壞 |

★ = 論文的兩個 headline setting。

### 三條注入路徑的差別

- **`replace_C`**（4/6/8）：把非平穩的相關矩陣直接交給 `rpotspatTS(cov.type = "precomputed")`。
  因為建構子保證**單位對角**，skew-t 的**邊際完全不變**，只有相依變成非平穩 → 分數差異可以
  乾淨地歸因於相依結構。這三個 setting 的**真均值恆為 10**（完全沒有均值結構）。
- **`additive`**（3/5/7/9）：在 skew-t 抽樣上加一個 mean-zero 的隨機場。邊際會被改變
  （skew-t 與該場的卷積），相依與邊際的效果混在一起。
- **`skew_term`**（10）：直接改 skew 項的係數。

### A/B 配對

Settings **4/5** 與 **6/7** 是配對對照：同一個非平穩結構，一個走 `replace_C`、一個走 `additive`，
用來拆開「邊際被改變」與「相依被改變」各自的影響。配對的兩個 setting **共用 base seed**
（`pair.seed = c(1,2,3,4,4,6,6,8,9,10)`），所以底層 skew-t 抽樣（`tau`/`z`/`knots` 與標準常態
innovation）**完全相同**，唯一的差別就是注入路徑。已驗證：`identical(tau.t[[4]], tau.t[[5]])` 為 `TRUE`。

---

## 各 setting 的數學定義

### Settings 1–3（沿用舊版，位元級不變）

$$\text{1: } g(s) = 5f_1 + 3f_2 \qquad \text{2: } g_t(s) = w_1(t)f_1 + w_2(t)f_2,\ (w_1,w_2)\sim N(0, \mathrm{diag}(25,9))$$

$$\text{3: } u_t(s) = \sum_{j=1}^{2}\xi_{tj}\cos(\kappa_j\|s - c_j\|),\quad \xi_{tj}\sim N(0,\tau_j^2),\ (\tau_1,\tau_2)=(5,3)$$

Setting 3 的 $\kappa = (\pi/10,\ \pi/5)$、中心 $(0,10)$ 與 $(7.5, 2.5)$，寫在原生的 $[0,10]^2$ 座標上。
積掉 $\xi$ 後得到低秩非平穩共變異數 $\sum_j \tau_j^2\phi_j(s)\phi_j(s')$，逐點變異數
$\sum_j \tau_j^2\phi_j(s)^2$ 隨位置變化。

**這是唯一一個相關會在長距離轉負的 setting**——低秩 cosine 基底做得到，單調的 Matérn 做不到。

### Settings 4/5：Paciorek–Schervish 局部各向異性

$$\Sigma(s) = R(\theta(s))\,\mathrm{diag}\big(a_1(s)^2, a_2(s)^2\big)\,R(\theta(s))^\top$$

$$a_1(s) = \rho_0 e^{\kappa_s u_2 + \kappa_a u_1},\quad a_2(s) = \rho_0 e^{\kappa_s u_2 - \kappa_a u_1},\quad \theta(s) = \frac{\pi}{2}\cdot\frac{s_2}{10}$$

其中 $u_k = s_k/10 - 0.5$，$\kappa_s = 0.8$（scale）、$\kappa_a = 1.0$（aspect）。三件事各自獨立變化：

| 性質 | 公式 | 沿哪個軸變 | 域內變化幅度 |
|------|------|-----------|-------------|
| **大小** | $\sqrt{\lvert\Sigma\rvert} = \rho_0^2 e^{2\kappa_s u_2}$ | $s_2$ | 2.2× |
| **長短軸比** | $a_1/a_2 = e^{2\kappa_a u_1}$ | $s_1$ | 7.4× |
| **方向** | $\theta(s)$ | $s_2$ | 0 → π/2 |

PS 共變異數（$\nu = 0.5$）：

$$C(s,s') = \lvert\Sigma(s)\rvert^{1/4}\lvert\Sigma(s')\rvert^{1/4}\left\lvert\tfrac{\Sigma(s)+\Sigma(s')}{2}\right\rvert^{-1/2}\exp\!\big(-\sqrt{Q}\big),\quad Q = (s-s')^\top\left[\tfrac{\Sigma(s)+\Sigma(s')}{2}\right]^{-1}(s-s')$$

正規化常數保證 $C(s,s) = 1$，正定性由 Paciorek & Schervish (2006) 證明。

> **⚠ `kappa_scale` 不能是 0。** 只讓 aspect 與 rotation 變化的話，$\lvert\Sigma(s)\rvert$ **恆為常數**
> ——橢圓只改形狀、面積不變。全向性的 variogram（以及**等向的擬合模型**）會把方向平均掉，
> 看到的是一個平穩過程。初版就是這個 bug：四象限 variogram 的分離度只有 0.025，比參考組的
> 噪音底線 0.031 還低。加上 scale 梯度後跳到 0.175（噪音底線的 5.7 倍）。

### Settings 6/7：共變數驅動的相依 ★

$$\rho(s) = \rho_0\exp\big(\beta_\rho f_1(s)\big),\qquad \rho_0 = 1,\ \beta_\rho = 0.7 \ \Rightarrow\ \rho(s)\in[0.50,\ 2.01]$$

等向 PS 核 $\Sigma(s) = \rho(s)^2 I_2$，代入後閉式化簡為：

$$C(s,s') = \frac{2\,\rho(s)\rho(s')}{\rho(s)^2 + \rho(s')^2}\exp\left(-\|s - s'\|\sqrt{\frac{2}{\rho(s)^2+\rho(s')^2}}\right)$$

單位對角自動成立，正定性由 PS 直接繼承（等向核是其特例）。`ps_cov_iso()` 已與通式
`ps_cov()` 交叉驗證，最大差 1.7e-16。

### Setting 8：Fuentes 區域混合

$$Y(s) = \sum_{k=1}^{4}w_k(s)Z_k(s), \qquad C(s,s') = \sum_k w_k(s)w_k(s')\exp\!\big(-\|s-s'\|/\rho_k\big)$$

中心置於四個象限中心，range $\rho = (0.4, 1.0, 2.0, 4.0)$（**10 倍差距**），權重
$w_k(s)\propto\exp(-\|s-m_k\|^2/(2h^2))$、$h = 3$，並以 $\sum_k w_k(s)^2 = 1$ 正規化 → 單位對角。

與 PS 的差別：PS 是**平滑漸變**，Fuentes 是**分區突變**（西南角粗糙、東北角平滑）。
$C = \sum_k \mathrm{diag}(w_k)C_k\mathrm{diag}(w_k)$ 為 PSD 矩陣之和，加 `gamma` nugget 後嚴格正定。

### Setting 9：Tukey g-and-h 隨機場

$$\tau_{g,h}(z) = \frac{e^{gz}-1}{g}e^{hz^2/2}\ (g\neq 0),\qquad \tau_{0,h}(z) = z\,e^{hz^2/2}$$

$$g(s) = 0.8\,f_1(s)\in[-0.80,\ 0.80], \qquad h(s) = 0.12\big(1 + f_2(s)\big)\in[0,\ 0.24]$$

$g$ 控偏態（**在域內變號**），$h$ 控尾重（西側高斯、東側重尾）。注入：

$$u_t(s) = \sigma_u\cdot\frac{\tau_{g(s),h(s)}\big(Z_t(s)\big) - m(s)}{v(s)},\qquad Z_t\sim N(0, C_{\exp})$$

$m(s), v(s)$ 以一次大樣本 Monte Carlo（`nsim = 2e5`）估出，使 $u_t$ **逐點 mean-zero、單位變異數**
→ 均值與變異數維持平穩，**只有高階動差**（與少量的相關結構）非平穩。

> **兩個實作陷阱。**
> (1) $h < 1/2$ 才有有限變異數；我們壓到 0.24，讓**四階動差也有限**——否則少數極端值會主導整個場，
> SNR 校準失去意義（初版 $h \le 0.3$ 時經驗偏態衝到 36）。
> (2) `gh_transform()` 裡的 `g` / `h` 必須**先廣播到 `z` 的形狀**再進 `ifelse()`：`ifelse()` 的
> 回傳長度取決於**測試條件**的長度，`g` 是純量時會把整個結果截成長度 1。

**$g$ 由 $f_1$ 驅動、$h$ 由 $f_2$ 驅動，刻意與 setting 10 的 $\lambda(s) = 3 + 3f_2$ 相反。**
若 $g$ 也騎在 $f_2$ 上，兩個 setting 的高階動差就是同一個空間 pattern，
`cor(skewness, g)` 與 `cor(skewness, lambda)` 會因構造而恆等，任何診斷都分不開它們。

### Setting 10：空間變動偏態 ★

$$\lambda(s) = 3 + 3f_2(s) \in [0,\ 6]$$

西側近乎對稱（$\lambda\approx0$）、東側強烈右偏（$\lambda\approx6$）。因為 `nknots = 1`，
$z_t$ 是**每日一個純量**、$\mu_t = x'\beta + \lambda z_t$，所以可以完全 post-hoc 調整、不需動 generator：

```r
y_new[, t] <- y_old[, t] + (lambda.field - 3) * data$z[1, t]
```

**為什麼是「部分撿到」的中間態**：$E[Y(s)] = 10 + \lambda(s)E[z_t]$，所以 $\lambda(s)$ 會
**誘導出一個空間變動的均值**（實測範圍 [10.0, 24.4]，跨度 14.4，與 setting 1 的真均值面
[2.05, 17.9] 幾乎同量級）→ MRTS 撿得到這一塊，`mr` 會下降。但預測分布的**偏態異質性**
（高分位數的形狀）抓不到 → Brier/quantile 在尾端只會部分改善。這填補了 setting 1（全撿到）
與 setting 3（全撿不到）之間的空白。

---

## 診斷結果

由 [diagnose_nonsta.R](diagnose_nonsta.R) 產出。關鍵指標是 **`vgm_vs_ref`**：四象限經驗
variogram 在 lag ≈ 1.5 的分離度，除以**平穩參考組**的同一指標（噪音底線）。$\gg 1$ 才是真的非平穩。

| ID | `surf_type` | truemean_sd | pw_sd_ratio | **vgm_vs_ref** | min_corr | 判讀 |
|----|-------------|-------------|-------------|----------------|----------|------|
| 1 | invariant | 3.68 | 1.40 | **0.96** | −0.047 | 相依平穩 ✓（只有均值） |
| 2 | varying | 0.75 | 1.35 | **0.67** | −0.046 | 相依平穩 ✓ |
| 3 | ns_dependence | 0.31 | **3.12** | **3.73** | **−0.170** | ✓ 且**相關轉負** |
| 4 | ps_cov | **0.00** | 1.60 | **5.70** | −0.048 | ✓ |
| 5 | ps_add | 0.37 | 1.32 | **4.17** | −0.045 | ✓ |
| 6 | covdep_cov | **0.00** | 1.21 | **9.14** | −0.034 | ✓ |
| 7 | covdep_add | 0.42 | 1.12 | **6.42** | −0.039 | ✓ |
| 8 | fuentes_cov | **0.00** | 1.32 | **26.6** | −0.077 | ✓ 最強 |
| 9 | gh_marginal | 0.38 | 1.86 | **6.35** | −0.042 | ✓ |
| 10 | lambda_varying | **2.76** | 1.15 | **0.75** | −0.045 | 相依平穩 ✓（只有偏態） |
| — | **ref_stationary** | 0.00 | 1.28 | **1.00** | −0.048 | 噪音底線 |

診斷有鑑別力：該有相依非平穩的（3–9）全部 $\ge 3.7$，不該有的（1、2、10）全部落在噪音底線。

### 高階動差指紋

| ID | skew_range | kurt_range | `cor(skew, g)` | `cor(skew, λ)` | `cor(kurt, h)` |
|----|-----------|-----------|----------------|----------------|----------------|
| 9 `gh_marginal` | **5.36** | **25.7** | **0.952** | 0.236 | **0.584** |
| 10 `lambda_varying` | 0.77 | 2.04 | 0.309 | **0.970** | −0.230 |
| ref_stationary | 0.22 | 0.57 | 0.077 | 0.044 | 0.081 |

兩個 shape setting 乾淨分離：setting 9 的經驗偏態跟著 $g(s)$（0.952），setting 10 跟著 $\lambda(s)$（0.970），
對角優勢明顯。剩餘的交叉項（0.24 / 0.31）來自 $f_1$ 與 $f_2$ 本身的相關，可接受。

### 邊際不變性（`replace_C` 的驗證）

與**同 seed** 的平穩抽樣（相同 `tau`/`z`、相同 innovation，只有相關矩陣不同）比較：

| ID | mean（ns vs ref） | sd（ns vs ref） | max \|Δquantile\|（1%–99.9%） | KS p |
|----|-------------------|-----------------|------------------------------|------|
| 4 | 15.4600 / 15.4593 | 7.429 / 7.432 | 0.18 | 0.864 |
| 6 | 15.4149 / 15.4185 | 6.756 / 6.770 | 0.11 | 0.075 |
| 8 | 15.3040 / 15.2982 | 6.400 / 6.379 | 0.73 | 0.627 |

**確認 `replace_C` 只改了相依、沒動邊際。**

### 診斷腳本的殘差設計（容易踩雷的地方）

模型是 $y_t(s) = 10 + m_t(s) + \lambda(s)z_t + \sigma_t\varepsilon_t(s) + u_t(s)$，其中 $z_t$ 是
**每日一個純量**（`nknots = 1`）。診斷用兩種殘差，因為兩類指標要的東西相反：

- **`R.dep`**（相依殘差）= $y - 10 - m_t(s) - \lambda(s)z_t$。用存下的真值（`z.t`、`W.varying`、
  `lambda.field`）**精確**扣掉截距、均值面與偏態項。$\lambda(s)z_t$ **必須扣掉**——它在第 $t$ 天
  對所有站點是同一個常數，會在每一對站點之間製造一個**與距離無關的正相關**，把 variogram 整個淹掉。
  → 用於 variogram、correlogram、逐點標準差。
- **`R.shape`**（形狀殘差）= $y - \texttt{truemean.field}$。只扣**時間平均**的真均值，**保留**偏態項
  → 用於逐點偏態／峰度。

  扣**真均值**（而不是每日空間平均）是關鍵：初版只扣每日空間平均，殘差裡留下的 $g(s)$ 常數項
  與逐日 $\sigma_t$ 標準化交互作用，在 setting 1（$\lambda$ 是常數、根本沒有偏態結構）上
  **偽造出 `cor_skew_g = −0.70` 的假訊號**，且把 setting 10 的符號翻反（−0.98）。改用真均值後
  分別降到 −0.065 與翻正到 **+0.970**。

兩者都再逐日除以空間標準差以移除 $\sigma_t$（重尾，df 3，否則少數幾天會主導所有平均）。

---

## `simdata_nonsta.RData` 的物件結構

```
y                [144 × 50 × 50 × 10]  觀測 [ns, nt, nsets, nsettings]
tau.t / z.t / knots.t   list[[10]]     各 setting 的潛在變數
settings.nonsta  data.frame [10 × 8]   setting, dist, nknots, lambda, surf_type,
                                       moment, route, reference
truemean.field   list[[10]]            ★ 每個元素 [144 × 50]：該 setting／該 dataset
                                       在所有站點的「時間平均真均值」
ns, nt, s, nsets, ntest, x             空間／時間 metadata

# settings 1-2
f.basis, c1, c2, surf.scale, surf.coef.invariant, W.varying, M.var
# setting 3
F.re, W.nsdep, kappa.re, centers.re, tau.re
# settings 4-10
C.stat, sigma.add
Sigma.ps, C.ps, ps.par
rho.field, C.covdep, covdep.par
C.fuentes, fuentes.par
g.field, h.field, gh.par, gh.mom
lambda.field, lambda.par
```

### `truemean.field`：統一的真值契約

`scores.R` 的 `true_mean_rec()` 直接查這張表，取代舊的 `surf_type` 字串 if/else 鏈。

**舊寫法有一個會靜默生效的 bug**：任何它不認得的 `surf_type` 都會掉進 `else` 分支
→ `f.basis %*% colMeans(W.varying[[set]])`，**不報錯**，只是算出完全錯誤的 recovery。
實測：對 setting 6（真均值恆為 10），舊的 fall-through 會給出跨度 [8.49, 11.51] 的假真值，
RMSE 誤差 1.18。數字看起來完全合理，但是錯的。

`truemean.field[[setting]][, set]` 的內容：

| setting | 值 |
|---------|-----|
| 1 | $10 + f_{\text{basis}}\,a_{\text{fixed}}$ |
| 2 | $10 + f_{\text{basis}}\,\overline{W}_{set}$ |
| 3 | $10 + F_{re}\,\overline{W3}_{set}$ |
| 4, 6, 8（`replace_C`） | $10$（常數——**完全沒有均值結構**） |
| 5, 7, 9（`additive`） | $10 + \overline{u}_{set}$ |
| 10 | $10 + \lambda(s)\cdot\overline{z}_{set}$ |

`scores.R` 保留了舊格式的 fallback 分支，所以備份的 `simdata_nonsta.pre10.bak.RData` 仍可重跑。

---

## 向後相容

`results_nonsta/` 的 fit 檔名以 setting id 當 key，`output/results/scores{1,2,3}_nonsta.RData`
是既有的快取分數。因此 **settings 1–3 的 `y` 必須位元級不變**：它們的程式碼路徑與 seed
（`set.seed(setting * 1000 + set)`、`900000 + set`、`750000 + set`）完全沒動。

已驗證：`identical(y_new[, , , 1:3], y_old)` 為 `TRUE`，且 `truemean.field` 與舊 `scores.R`
公式的差為 `0.000e+00`。新 setting 的輔助 seed 用互不重疊的 offset
（`810000` / `820000` / `830000`），診斷腳本的平穩參考組用 `990000`。

---

## 實作說明

- **`rpotspatTS(cov.type = "precomputed", C.mat = ...)`**（[auxfunctions.R](../../../R/ar2/auxfunctions.R)）：
  新增的分支，接受呼叫端算好的 $n_s\times n_s$ **相關**矩陣（單位對角）。嚴格加法式改動，
  `match.arg` 預設仍是 `"matern"`、新參數預設 `NULL`，對 simstudy / simstudy_prop /
  time_block_forecast 的所有現有 caller 零影響。
  > 註：`rpotspatTS_arp()` 的 `cov.type` 仍只有 `c("matern", "deformed")`，落後一版。
  > nonsta 全部 `phi.* = 0`，用不到 arp 版。
- **`to_correlation(C, gamma)`**（[ns_cov_builders.R](ns_cov_builders.R)）：`cov2cor()` → 乘 `gamma`
  → `diag <- 1` → **正定性斷言**（`min(eigen(...)) > 1e-8`，否則 `stop()`）。三個取代式矩陣的
  最小特徵值分別是 0.18（PS）、0.19（covdep）、0.14（Fuentes）。
- 所有 $144\times144$ 稠密矩陣直接建，不需任何優化。

---

## 參考

> ⚠ 標題可靠；年份／期刊／卷期為記憶重建，**寫進 bibliography 前請逐筆查證**。

- Cressie, N. and Johannesson, G. (2008). Fixed rank kriging for very large spatial data sets. *JRSS-B*, **70**(1), 209–226. — setting 3
- Fuentes, M. (2001). A high frequency kriging approach for non-stationary environmental processes. *Environmetrics*, **12**, 469–483. — setting 8
- Fuentes, M. (2002). Spectral methods for nonstationary spatial processes. *Biometrika*, **89**, 197–210. — setting 8
- Higdon, D., Swall, J. and Kern, J. (1999). Non-stationary spatial modeling. *Bayesian Statistics 6*. — settings 4/5 的源頭
- Morris, S.A., Reich, B.J., Thibaud, E. and Cooley, D. (2017). A space-time skew-t model for threshold exceedances. *Biometrics*, **73**(3), 749–758. — 基礎模型；setting 10 破壞其「$\lambda$ 為常數」假設
- Paciorek, C.J. and Schervish, M.J. (2006). Spatial modelling using a new class of nonstationary covariance functions. *Environmetrics*, **17**(5), 483–506. — settings 4/5/6/7
- Risser, M.D. and Calder, C.A. (2015). Regression-based covariance functions for nonstationary spatial modeling. *Environmetrics*. — settings 6/7
- Schmidt, A.M., Guttorp, P. and O'Hagan, A. (2011). Considering covariates in the covariance structure of spatial processes. *Environmetrics*, **22**, 487–500. — settings 6/7
- Tzeng, S. and Huang, H.-C. (2018). Resolution adaptive fixed rank kriging. *Technometrics*. — settings 1/2/3
- Xu, G. and Genton, M.G. (2017). Tukey g-and-h random fields. *JASA*, **112**(519), 1236–1249. — setting 9
