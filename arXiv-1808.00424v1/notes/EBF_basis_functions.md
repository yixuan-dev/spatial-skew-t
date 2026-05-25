# Empirical Basis Functions (EBF) 筆記

> 對應論文:Morris, Reich & Thibaud (2018),*Exploration and inference in spatial extremes using empirical basis functions* (arXiv:1808.00424v1)
> 原始檔:[EBF_arxiv_v1.tex](../EBF_arxiv_v1.tex)

---

## 1. 模型假設(Low-rank max-stable model)

對空間位置 $\mathbf{s}\in\mathcal{S}\subset\mathbb{R}^2$、時間 $t$ 觀測 $Y_t(\mathbf{s})$,先把邊際轉成單位 Fréchet 得殘差過程 $Z_t(\mathbf{s})$,然後採用 Reich & Shaby (2012) 的 **低秩正穩定(positive stable)max-stable 模型**:

$$
Z_t(\mathbf{s})=\theta_t(\mathbf{s})\,\epsilon_t(\mathbf{s}),\qquad
\theta_t(\mathbf{s})=\Big\{\sum_{l=1}^L B_l(\mathbf{s})^{1/\alpha}A_{lt}\Big\}^{\alpha}
$$

關鍵假設:

- **有限秩** $L$:用 $L$ 個 basis functions 近似 spectral representation。
- **基底限制**:$B_l(\mathbf{s})\ge 0$ 且 $\sum_{l=1}^L B_l(\mathbf{s})=1$。
- **隨機效應**:$A_{lt}\overset{\text{iid}}{\sim}\text{PS}(\alpha)$(正穩定分布),$\alpha\in(0,1)$。
- **噪音項**:$\epsilon_t(\mathbf{s})\overset{\text{iid}}{\sim}\text{GEV}(1,\alpha,\alpha)$,使 $Z_t(\mathbf{s})$ 仍為單位 Fréchet max-stable。
- $\alpha$ 扮演 nugget:$\alpha\to 0$ 退化為 max-linear、$\alpha\to 1$ 退化為獨立。

成對 extremal coefficient:

$$
\vartheta(\mathbf{s}_1,\mathbf{s}_2)=\sum_{l=1}^{L}\big\{B_l(\mathbf{s}_1)^{1/\alpha}+B_l(\mathbf{s}_2)^{1/\alpha}\big\}^{\alpha},\qquad
\lim_{\|\mathbf{s}_1-\mathbf{s}_2\|\to 0}\vartheta=2^{\alpha}.
$$

---

## 2. Basis functions 是什麼?

- **角色**:$B_1(\mathbf{s}),\ldots,B_L(\mathbf{s})$ 是非負、和為 1 的空間函數,各自代表一種「主導空間模式」。
  - Smith (1990) storm 解釋:$B_l(\mathbf{s})$ 是第 $l$ 種風暴的空間範圍,$A_{lt}$ 是該年的強度。
- **與 PCA / EOF 對比**:同為 dimension reduction,但 EBF **不正交、非負**,精神接近 dictionary learning / NMF。
- **與 GKF 對比**:Reich & Shaby (2012) 使用固定形式的高斯核
  $$
  B_l(\mathbf{s})=\frac{\exp\{-(\|\mathbf{s}-\mathbf{k}_l\|/\rho)^2\}}{\sum_{j=1}^L \exp\{-(\|\mathbf{s}-\mathbf{k}_j\|/\rho)^2\}};
  $$
  本論文則 **完全由資料估出** $B_l$,所以叫 *empirical* basis functions,自然允許非平穩性。
- **重要性度量**:
  $$
  v_l=\frac{1}{n_s}\sum_{i=1}^{n_s}\hat B_{il},\qquad \sum_{l=1}^L v_l=1.
  $$
  按 $v_l$ 由大到小排序,作為類似主成分的詮釋。

---

## 3. 如何估計 $\hat{\alpha}$ 與 $\hat{B}_{il}$

三步演算法:

### Step 1 — 初估 extremal coefficient
對所有點對 $(i,j)$ 以 Cooley (2006) 的 F-madogram 估計

$$
\hat\nu^{F}(\mathbf{s}_i,\mathbf{s}_j)=\frac{1}{2n_t}\sum_{t=1}^{n_t}\big|\hat F_i(y_{it})-\hat F_j(y_{jt})\big|,\qquad
\hat\vartheta_{ij}=\frac{1+2\hat\nu^F}{1-2\hat\nu^F}.
$$

### Step 2 — 高斯核空間平滑

$$
\tilde\vartheta_{ij}=\frac{\sum_{u,v}w_{iu}w_{jv}\hat\vartheta_{uv}}{\sum_{u,v}w_{iu}w_{jv}},\qquad
w_{iu}=\exp\{-(\|\mathbf{s}_i-\mathbf{s}_u\|/\delta)^2\},\ w_{ii}=0.
$$

把 $w_{ii}$ 設為 0,以避免人造的 $\hat\vartheta_{ii}=1$ 拉低估計值,並消除原點不連續的影響。

### Step 3 — 兩階段估計

- **先估 $\alpha$**:取空間距離很近的點對集合 $\mathcal{N}$,利用 $\vartheta\approx 2^\alpha$,
  $$
  \hat\alpha=\log_2\!\Big(\tfrac{1}{|\mathcal{N}|}\textstyle\sum_{\mathcal{N}}\tilde\vartheta_{ij}\Big).
  $$
- **再估 $\mathbf{B}$**:固定 $\hat\alpha$,以 **梯度下降** 最小化
  $$
  \sum_{i<j}\Big\{\tilde\vartheta_{ij}-\sum_{l=1}^L\big(B_{il}^{1/\hat\alpha}+B_{jl}^{1/\hat\alpha}\big)^{\hat\alpha}\Big\}^2,
  $$
  受約束 $B_{il}\ge 0$、$\sum_l B_{il}=1$。

### 調參

- $\delta$:用 extremal coefficient 的交叉驗證選擇。
- $L$:可用誤差函數的 elbow,或以預測 CV 評估不同 $L$。

### 後續 Bayesian 推論

固定 $\hat{\mathbf{B}},\hat\alpha$,僅用 MCMC 估 $A_{lt}$ 與邊際參數;正穩定密度以一維數值積分(midpoint rule + Beta(0.5,0.5) 50 個分位點)近似。

---

## 4. 模擬實驗的設計與發現

### 設計

- $n_s=100$ 個位置均勻取自 $[1,10]\times[1,10]$。
- 真實模型為 Reich & Shaby (2012),用 GKF($\rho=2.5$)在 $\sqrt L\times\sqrt L$ 格點作為真基底。
- 因子設計:$L\in\{9,25\}$、$\alpha\in\{0.3,0.7\}$、$n_t\in\{50,200\}$,共 8 個情境,每情境 100 次重複。
- 每個情境用 $L_{\text{fit}}=9$ 與 $L_{\text{fit}}=25$ 兩種設定配適,比較 4 種 extremal coefficient 估計器的 MSE:
  1. F-madogram 原始估計 $\hat\vartheta_{ij}$
  2. 平滑後 $\tilde\vartheta_{ij}$
  3. EBF($L=9$)
  4. EBF($L=25$)

### 結果(節錄 Table 1)

| $L$ | $\alpha$ | $n_t$ | $\hat\alpha$ (sd) | Initial | Smoothed | EBF $L=9$ | EBF $L=25$ |
|----:|---------:|------:|------------------:|--------:|---------:|----------:|-----------:|
|  9 | 0.3 |  50 | 0.31 (0.02) | 1.08 | 0.87 | **0.84** | 0.88 |
|  9 | 0.3 | 200 | 0.31 (0.01) | 0.28 | 0.28 | **0.27** | 0.31 |
|  9 | 0.7 |  50 | 0.70 (0.04) | 1.19 | 0.49 | **0.50** | 0.51 |
|  9 | 0.7 | 200 | 0.70 (0.02) | 0.31 | 0.14 | **0.13** | 0.15 |
| 25 | 0.3 |  50 | 0.32 (0.02) | 1.12 | 0.89 | 2.08 | **0.95** |
| 25 | 0.3 | 200 | 0.32 (0.01) | 0.27 | 0.29 | 1.56 | **0.38** |
| 25 | 0.7 |  50 | 0.70 (0.03) | 1.12 | 0.55 | 0.84 | **0.60** |
| 25 | 0.7 | 200 | 0.70 (0.01) | 0.27 | 0.16 | 0.46 | **0.20** |

### 主要發現

1. **$\hat\alpha$ 在所有情境都很準**(偏差小、SD ≤ 0.04)。
2. **平滑步驟通常有幫助**;只有「強相依 $\alpha=0.3$ + 大樣本 $n_t=200$」時原始估計已夠好。
3. **指定正確 $L$**:EBF 的 MSE 與平滑估計相當 → 少量基底即可有效率地表徵相依結構。
4. **指定過多基底**($L_{\text{fit}}=25,\;L_{\text{true}}=9$):MSE 僅略上升,代價小。
5. **指定過少基底**($L_{\text{fit}}=9,\;L_{\text{true}}=25$):MSE 顯著上升(例如 $\alpha=0.3,n_t=200$ 由 0.38 → 1.56)。
6. **實務建議**:寧可多放一些基底,不要放太少。

### 與實證的呼應

對 NARCCAP 美東降水分析:$L=10$ 已足以掌握相依結構;在未來氣候的強相依情境,EBF 的預測 MAD 較 GKF 更低。
