# Huser & Wadsworth (2022)：空間極值相依模型與資料生成

> Huser, R. and Wadsworth, J.L. (2022). Advances in statistical modeling of spatial extremes.
> *Annual Review of Statistics and Its Application*, **9**: 401–431.
>
> 本文聚焦於**如何用該綜述所涵蓋的模型族生成非平穩極值相依資料**，而非估計。

---

## 核心問題：極值相依的分類

- 設 $Z(s_1), Z(s_2)$ 的邊際 CDF 均為 $F$，定義**尾部相依係數**（tail dependence coefficient）：

$$\chi(s_1, s_2) = \lim_{u \to 1^-} \Pr\!\left(F(Z(s_1)) > u \mid F(Z(s_2)) > u\right)$$

- $\chi > 0$：**漸進相依（Asymptotic Dependence, AD）**，極端事件傾向共現。
- $\chi = 0$：**漸進獨立（Asymptotic Independence, AI）**，極端事件趨向獨立，但仍可有中等相關。
- 補充量 $\bar{\chi} = \lim_{u\to 1^-} 2\log \Pr(F(Z(s))>u) / \log \Pr(F(Z(s_1))>u, F(Z(s_2))>u) - 1 \in [-1,1]$，在 AI 情況下刻畫相依強度。
- **Non-stationary**：若 $\chi(s_1, s_2)$ 不僅依賴 $s_1 - s_2$ 而是依賴 $s_1, s_2$ 各自的位置，即為非平穩極值相依。

---

## 一、Max-stable 過程（AD 框架）

### 譜表示（Spectral Representation）

- 任何 max-stable 過程有以下表示，其中 $\{(\xi_i, Y_i)\}$ 為 Poisson 過程（強度 $\xi^{-2} d\xi$）上的原子：

$$Z(s) = \max_{i \geq 1} \xi_i \cdot Y_i(s), \qquad \mathbb{E}\!\left[\max(Y(s), 0)\right] = 1$$

- 不同的 $Y_i(s)$（spectral function）給出不同的 max-stable 模型族，但都有 $\chi > 0$，屬 AD。

### Brown-Resnick 過程（最常用）

- 令 $\{W_i(s)\}$ 為 iid Gaussian 過程，具有半變異數 $\gamma(s_1, s_2) = \frac{1}{2}\text{Var}(W(s_1) - W(s_2))$，則：

$$Z(s) = \max_{i \geq 1} \xi_i \exp\!\left(W_i(s) - \gamma(s_0, s)\right)$$

- $\gamma$ 為**本質平穩（intrinsically stationary）**時，$Z$ 為平穩 max-stable；$\gamma$ 為非平穩時，$\chi(s_1,s_2)$ 隨位置變化：

$$\chi(s_1,s_2) = 2\left(1 - \Phi\!\left(\sqrt{\gamma(s_1,s_2)/2}\right)\right)$$

- 常見平穩半變異數：分數 Brown 運動（fBm），$\gamma(h) = \|h\|^\alpha / \ell^\alpha$，$\alpha \in (0,2]$，$\ell > 0$。

#### 生成演算法（截斷 Poisson 過程）

```r
# Brown-Resnick 生成（截斷近似）
library(MASS); library(fields)

rbrown_resnick <- function(s, gamma_fn, n_iter = 2000) {
  n  <- nrow(s)
  Z  <- rep(-Inf, n)
  xi_cum <- 0
  for (i in seq_len(n_iter)) {
    xi     <- rexp(1)            # Poisson 過程間距
    xi_cum <- xi_cum + xi
    zeta   <- 1 / xi_cum         # ξ_i ~ Pareto(1)
    # 模擬 W_i：以 Gaussian 過程（均值 0，協方差由半變異數推）
    Sigma_W <- outer(1:n, 1:n, Vectorize(
      function(i, j) gamma_fn(s[i,], s[n,]) + gamma_fn(s[j,], s[n,]) - gamma_fn(s[i,], s[j,])
    ))
    W <- mvrnorm(1, rep(0, n), Sigma_W)
    v <- gamma_fn(s[n,], s[n,])  # = 0，參考點
    Y <- zeta * exp(W - diag(Sigma_W)/2)
    Z <- pmax(Z, Y)
  }
  Z
}

# 半變異數（fBm，平穩）
gamma_fBm <- function(s1, s2, alpha = 1.5, ell = 3) {
  (sum((s1 - s2)^2))^(alpha/2) / (2 * ell^alpha)
}
```

### Extremal-t 過程（Opitz 2013）

- 令 $\{T_i(s)\}$ 為 iid Student-t 過程（自由度 $\nu$，相關函數 $\rho$），則：

$$Z(s) = c_\nu \max_{i \geq 1} \xi_i \max(0, T_i(s))^\nu, \quad c_\nu = \sqrt{\pi}\,\frac{2^{(\nu-2)/2}\,\Gamma\!\left(\frac{\nu+1}{2}\right)}{\Gamma\!\left(\frac{\nu}{2}\right)}$$

- 尾部相依係數：$\chi(s_1,s_2) = 2\,\bar{F}_{t_{\nu+1}}\!\left(\sqrt{(\nu+1)(1-\rho(s_1-s_2))/(1+\rho(s_1-s_2))}\right)$，其中 $\bar{F}_{t_k}$ 為 t 分布的生存函數。
- $\nu \to \infty$（$T_i \to \text{Gaussian}$）：趨向 Schlather 過程；$\nu \to 0$：趨向完全相依。

---

## 二、Scale Mixture 模型族（AD ↔ AI 橋接）

### 一般框架

- 設 $W(s)$ 為**標準 Gaussian 過程**，相關函數 $\rho(s_1, s_2)$；$R > 0$ 為獨立的隨機尺度：

$$Y(s) = \mu + \sigma \cdot R \cdot W(s)$$

- 邊際分布由 $R$ 決定；$R$ 的尾部重量決定極值相依類別：
  - $R \equiv 1$：Gaussian 過程 → AI（$\chi = 0$）
  - $R \sim \text{Inv-Gamma}(\nu/2, \nu/2)$：Student-t 過程 → AD（$\chi > 0$），$\chi = 2\bar{F}_{t_{\nu+1}}(\cdots)$
  - $R \sim \text{LogNormal}$：AI 但重尾邊際
  - $R \sim \text{Pareto}$：接近 max-stable 行為（強 AD）

### 控制 AD/AI 的連續橋接（Huser, Davison & Genton 2017）

- 引入參數 $\delta \in (0,1]$，構造：

$$Y(s) = R^{1/\delta} \cdot W(s), \qquad R \sim \text{Pareto}(1)$$

- $\delta = 1$：$Y$ 有 Pareto 邊際，AD；$\delta \to 0$：邊際變重但相依趨向 AI；$\delta$ 連續控制 AD/AI 過渡。
- 生成資料：先抽 $R$，再以 $R^{1/\delta}$ 縮放 Gaussian 過程的一個實現。

```r
# Scale mixture 生成（AD/AI 可控）
library(MASS)

r_scale_mixture <- function(s, rho_fn, delta = 0.5, n_rep = 1) {
  n     <- nrow(s)
  Sigma <- outer(1:n, 1:n, Vectorize(function(i,j) rho_fn(s[i,], s[j,])))
  W     <- mvrnorm(n_rep, rep(0, n), Sigma)   # Gaussian 過程
  R     <- 1 / runif(n_rep)^1                  # Pareto(1)
  R^(1/delta) * W                              # shape: n_rep x n
}
```

### Student-t 過程（特例，最常用）

- $R^2 \sim \text{Gamma}(\nu/2, \nu/2)$（即 $R^2$ 為 scaled chi-squared），生成：

```r
r_t_process <- function(s, rho_fn, nu = 4, n_rep = 1) {
  n     <- nrow(s)
  Sigma <- outer(1:n, 1:n, Vectorize(function(i,j) rho_fn(s[i,], s[j,])))
  R2    <- rgamma(n_rep, nu/2, nu/2)
  W     <- mvrnorm(n_rep, rep(0, n), Sigma)
  W / sqrt(R2)                # 等價於 t_ν 過程
}
```

---

## 三、廣義 Pareto 過程（Threshold Exceedance 框架）

- Max-stable 對應**區塊最大值**；廣義 Pareto 過程（GPP）對應**閾值超越**（peaks-over-threshold）：

$$\Pr\!\left(\frac{Z(s) - u}{\sigma(s)} \leq z \mid \max_s Z(s) > u\right) \to H(z; \xi)$$

- de Fondeville & Davison (2018) 建立了 GPP 的生成框架，以 Poisson 過程 + 正規化取得超越樣本。
- **生成策略**：① 生成 max-stable 過程一個實現，② 條件化在某站超越閾值，③ 以正規化取得 GPP 樣本；計算成本高，通常用於推論而非生成研究。

---

## 四、Non-stationary 相依的生成策略

### 策略 A：Non-stationary 相關函數（最直接）

- 將上述所有模型中的相關函數 $\rho(s_1, s_2)$ 替換為**非平穩**版本（如 Paciorek-Schervish）：

$$\rho(s_1, s_2) = |\Sigma(s_1)|^{1/4}|\Sigma(s_2)|^{1/4}\left|\frac{\Sigma(s_1)+\Sigma(s_2)}{2}\right|^{-1/2} \phi(Q_{12})$$

- $Q_{12} = (s_1-s_2)^\top\left[\frac{\Sigma(s_1)+\Sigma(s_2)}{2}\right]^{-1}(s_1-s_2)$，$\Sigma(s)$ 為位置依賴的局部 kernel 矩陣。
- **效果**：極值相依強度 $\chi(s_1,s_2)$ 因 $\rho(s_1,s_2)$ 的非平穩性而隨位置變化。

```r
# Paciorek-Schervish 非平穩相關
ps_corr <- function(s1, s2, Sigma_fn, phi_fn) {
  S1 <- Sigma_fn(s1); S2 <- Sigma_fn(s2)
  Sm <- (S1 + S2) / 2
  det_factor <- det(S1)^0.25 * det(S2)^0.25 / det(Sm)^0.5
  Q  <- as.numeric(t(s1-s2) %*% solve(Sm) %*% (s1-s2))
  det_factor * phi_fn(sqrt(Q))
}
```

### 策略 B：空間變形（Deformation）

- 設計非線性雙射 $D: \mathcal{G} \to \mathcal{F}$，在 F-space 中用**平穩**模型，拉回 G-space 即為非平穩。
- Brown-Resnick：$\gamma(s_1, s_2) = \|D(s_1) - D(s_2)\|^\alpha / (2\ell^\alpha)$（Blanchet & Davison 2011 的做法）。
- Scale mixture / t-process：$\rho(s_1, s_2) = \rho_0(\|D(s_1) - D(s_2)\|)$。
- $D$ 的設計詳見 [sampson_guttorp_1992.md](sampson_guttorp_1992.md)。

```r
# 以變形法建立非平穩 t-process
D <- make_D_axial(alpha = 0.15)   # 來自 sampson_guttorp_1992.md
f <- D(s)                          # 映射至 F-space
rho_fn_ns <- function(si, sj) exp(-sqrt(sum((D(t(si)) - D(t(sj)))^2)) / 3)
Y <- r_t_process(s, rho_fn = rho_fn_ns, nu = 4, n_rep = 50)
```

### 策略 C：空間變化的尾部參數（最彈性）

- 讓控制 AD/AI 的參數 $\delta(s)$ 或自由度 $\nu(s)$ 隨位置變化，生成**不同區域有不同相依類別**的過程。
- 例：$\nu(s) = \nu_0 \exp(\beta^\top s)$，使某些區域呈 AD，其他區域趨向 AI。
- **注意**：空間變化的 $\nu(s)$ 或 $\delta(s)$ 不再有乾淨的聯合分布閉合式，通常需以**條件模擬**（逐站條件化）或 copula 框架近似生成。

---

## 五、完整生成流程（以非平穩 t-process 為例）

```r
library(MASS); library(fields)

# 1. 站點
set.seed(42)
s <- cbind(runif(100, 0, 10), runif(100, 0, 10))

# 2. 設計非平穩相關（變形法）
D      <- make_D_axial(alpha = 0.12)
f      <- D(s)
distF  <- as.matrix(dist(f))
rho    <- exp(-distF / 2.5)           # 指數相關，在 F-space 平穩

# 3. 抽樣（t-process，nu = 5）
nu     <- 5
n_rep  <- 200
R2     <- rgamma(n_rep, nu/2, nu/2)   # 每個重複一個共同尺度
W      <- mvrnorm(n_rep, rep(0, nrow(s)), rho)
Y      <- W / sqrt(R2)                # shape: n_rep x n_sites

# 4. 驗證：chi 應隨 D(s_i) - D(s_j) 的距離單調遞減
# （但不應只隨 s_i - s_j 單調，因為有非平穩性）
```

---

## 六、模型族比較

| 模型 | 相依類別 | 非平穩擴展 | 生成難度 | 適合場景 |
|------|----------|-----------|---------|---------|
| Brown-Resnick | AD | 非平穩 $\gamma(s_1,s_2)$ | 中 | 降水最大值、風速 |
| Extremal-t | AD | 非平穩 $\rho(s_1,s_2)$ | 易 | 溫度極值 |
| Scale mixture（t-process） | AD | 同上 | 易 | 一般用途 |
| Scale mixture（$\delta$ 參數） | AD ↔ AI | $\delta(s)$ 空間變化 | 中 | 需橋接兩類時 |
| Gaussian process | AI | 非平穩 cov | 易 | 輕尾、作為對照組 |
| GPP | AD | 較複雜 | 高 | 閾值超越分析 |

- **本 project 的 spatial skew-t**：屬 scale mixture 族的變體，以 $\tau(s)$（局部尺度）取代全域的 $R$，本質上是**位置依賴的尺度混合**，自然地實現了非平穩極值相依。
- 生成非平穩極值資料的最低成本路線：**非平穩 $\rho$ 的 t-process**（策略 A 或 B）+ `MASS::mvrnorm` + 隨機尺度 $R$。

---

## 參考

- Huser, R. and Wadsworth, J.L. (2022). Advances in statistical modeling of spatial extremes. *Annual Review of Statistics and Its Application*, **9**: 401–431.
- Huser, R., Davison, A.C. and Genton, M.G. (2017). Modeling spatial processes with unknown extremal dependence class. *Journal of the American Statistical Association*, **112**(519): 1246–1261.
- Blanchet, J. and Davison, A.C. (2011). Spatial modeling of extreme snow depths. *Annals of Applied Statistics*, **5**(3): 1699–1725.
- de Fondeville, R. and Davison, A.C. (2018). High-dimensional peaks-over-threshold inference. *Biometrika*, **105**(3): 575–592.
- Opitz, T. (2013). Extremal-t processes: Elliptical domain of attraction and a spectral representation. *Journal of Multivariate Analysis*, **122**: 409–413.
- Paciorek, C.J. and Schervish, M.J. (2006). Spatial modelling using a new class of nonstationary covariance functions. *Environmetrics*, **17**(5): 483–506.
