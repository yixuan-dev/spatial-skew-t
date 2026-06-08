# 非線性變形法生成非平穩空間資料

> 基於 Sampson & Guttorp (1992) 框架。本文聚焦於**如何設計變形函數 $D$ 以生成資料**，而非估計問題。

---

## 框架與目標

- 設空間過程 $\{Z(s) : s \in \mathcal{G} \subset \mathbb{R}^2\}$，目標是讓它在 G-space 中呈現**非平穩**協方差結構。
- 核心思想：找雙射 $D: \mathcal{G} \to \mathcal{F}$，使得 Z 在 F-space 中是等向平穩的：

$$C(s_1, s_2) = \sigma^2\,\phi\!\left(\|D(s_1) - D(s_2)\|\right)$$

- 生成資料只需三步：① 設計 $D$，② 對觀測站計算 F-space 距離矩陣，③ 從 MVN 抽樣。
- **估計問題（MDS + TPS 回推 $D$）不在本文範圍內**；此處 $D$ 由研究者主動設計。

---

## $D$ 的必要條件

- **雙射性**：$D$ 須一對一，否則不同 G-space 位置在 F-space 重疊，協方差矩陣退化。
- **實用判準**：在 $\mathcal{G}$ 上處處有

$$\det(\nabla D(s)) > 0$$

- **平滑性**：$D \in C^1(\mathcal{G})$ 即足夠；TPS 自動滿足 $C^2$。
- **非線性**：若 $D(s) = As + b$（線性），退化為幾何各向異性（stationary）；非線性才有真正的 non-stationarity。

---

## 方法一：封閉形式參數化變形

### 1a. 軸向非線性縮放

$$D(s_1, s_2) = \bigl(s_1,\; g(s_1)\cdot s_2\bigr), \quad g: \mathbb{R} \to \mathbb{R}^+$$

- Jacobian 行列式 $= g(s_1) > 0$，**解析保證雙射**。
- $g(s_1) = e^{\alpha s_1}$：$s_1$ 越大，$s_2$ 方向在 F-space 中越壓縮 → 該區域的 $s_2$ 方向相關衰減更慢。
- $\alpha > 0$：右側壓縮；$\alpha < 0$：左側壓縮；$|\alpha|$ 控制非平穩程度。

```r
make_D_axial <- function(alpha) {
  function(s) cbind(s[, 1], exp(alpha * s[, 1]) * s[, 2])
}
D <- make_D_axial(alpha = 0.15)
```

### 1b. 正弦波橫向位移

$$D(s_1, s_2) = \bigl(s_1 + a\sin(b\, s_2),\; s_2\bigr)$$

- Jacobian 行列式 $= 1 > 0$，**無條件雙射**。
- $a$ 控制位移幅度，$b$ 控制週期；使等相關輪廓沿 $s_1$ 方向呈波浪狀扭曲。

```r
make_D_sine <- function(a, b) {
  function(s) cbind(s[, 1] + a * sin(b * s[, 2]), s[, 2])
}
D <- make_D_sine(a = 1.5, b = 0.4)
```

### 1c. 徑向推拉（Radial Push/Pull）

$$D(s) = s + h\!\left(\|s - c\|\right) \cdot \frac{s - c}{\|s - c\|}$$

- $c \in \mathcal{G}$ 為中心點，$h(r)$ 為平滑徑向函數（如 $h(r) = \beta r e^{-r/\ell}$）。
- 雙射條件：需 $1 + h'(r) > 0$ 且 $1 + h(r)/r > 0$（數值驗證，見下節）。
- 直觀效果：鄰近 $c$ 的點被推離（$h > 0$）或拉近（$h < 0$），使該區域的相關結構與外圍不同。

```r
make_D_radial <- function(c, beta, ell) {
  function(s) {
    d  <- sqrt((s[,1]-c[1])^2 + (s[,2]-c[2])^2)
    h  <- beta * d * exp(-d / ell)
    dx <- ifelse(d > 0, h * (s[,1]-c[1]) / d, 0)
    dy <- ifelse(d > 0, h * (s[,2]-c[2]) / d, 0)
    cbind(s[,1] + dx, s[,2] + dy)
  }
}
D <- make_D_radial(c = c(5, 5), beta = 0.2, ell = 3)
```

---

## 方法二：TPS + 控制點（最彈性）

- 在 G-space 中選 $K$ 個**控制點** $\{p_k\}_{k=1}^K$，手動指定它們在 F-space 的目標位置 $\{q_k\}_{k=1}^K$。
- 對 $D$ 的兩個座標分量分別擬合 TPS，使 $D(p_k) = q_k$：

$$\min_{D} \sum_{k=1}^K \|q_k - D(p_k)\|^2 + \lambda \int \|\nabla^2 D\|^2\, ds$$

- TPS 解析解：$D(s) = As + b + \sum_k c_k\,\psi(\|s - p_k\|)$，$\psi(r) = r^2 \ln r$。
- **設計原則**：將某區域的 $q_k$ 相對於 $p_k$ 集中（壓縮），該區域在 G-space 中的相關將更強；反之拉開則相關更弱。

```r
library(fields)

# 在 [0,10]^2 上設計 9 個控制點
P <- as.matrix(expand.grid(x = c(0, 5, 10), y = c(0, 5, 10)))
Q <- P                      # 從恆等映射出發
Q[5, ] <- c(7, 7)           # 把中心點推向右上 → 中心區域壓縮
Q[c(2,4,6,8), ] <- Q[c(2,4,6,8), ] * 0.9  # 邊中點向內收縮

tps1 <- Tps(P, Q[, 1])
tps2 <- Tps(P, Q[, 2])

D <- function(s) cbind(predict(tps1, s), predict(tps2, s))
```

---

## 方法三：串接簡單雙射

- 若 $D_1, D_2$ 各自雙射，則 $D = D_2 \circ D_1$ 亦雙射，且：

$$\det(\nabla D) = \det(\nabla D_2)\big|_{D_1(s)} \cdot \det(\nabla D_1)\big|_s > 0$$

- 策略：先用旋轉 $D_1$（調整主軸方向），再用軸向縮放 $D_2$（產生位置依賴的 range）。

```r
D1 <- function(s) {           # 旋轉 30 度
  theta <- pi / 6
  cbind(cos(theta)*s[,1] - sin(theta)*s[,2],
        sin(theta)*s[,1] + cos(theta)*s[,2])
}
D2 <- make_D_axial(alpha = 0.1)   # 軸向非線性縮放

D  <- function(s) D2(D1(s))       # 串接
```

---

## 完整生成流程

```r
library(MASS)

# 1. 觀測站位置
set.seed(42)
s <- cbind(runif(100, 0, 10), runif(100, 0, 10))

# 2. 選定 D（任一方法）
D <- make_D_axial(alpha = 0.15)

# 3. 映射至 F-space，計算距離矩陣
f     <- D(s)
distF <- as.matrix(dist(f))

# 4. 選定 φ（指數族）與參數
sigma2 <- 1
theta  <- 2
Sigma  <- sigma2 * exp(-distF / theta)

# 5. 模擬
Z <- mvrnorm(n = 50, mu = rep(0, nrow(s)), Sigma = Sigma)
#   Z[t, ] 為第 t 個重複的觀測向量
```

---

## 驗證雙射性：數值 Jacobian 檢查

- 對 $\mathcal{G}$ 上的網格點數值估計 $\det(\nabla D)$，確認全為正值。

```r
jac_det <- function(D, s, h = 1e-5) {
  e1 <- cbind(h, 0); e2 <- cbind(0, h)
  dD1 <- (D(s + e1) - D(s - e1)) / (2 * h)   # ∂D/∂s₁
  dD2 <- (D(s + e2) - D(s - e2)) / (2 * h)   # ∂D/∂s₂
  dD1[, 1] * dD2[, 2] - dD1[, 2] * dD2[, 1]  # det
}

grid <- as.matrix(expand.grid(s1 = seq(0, 10, 0.5), s2 = seq(0, 10, 0.5)))
jd   <- jac_det(D, grid)
all(jd > 0)   # 須為 TRUE；若有負值代表 D 在該區域折疊
range(jd)     # 觀察行列式的空間變化幅度
```

---

## 各方法比較

| 方法 | 雙射保證 | 非平穩彈性 | 直觀控制 | R 難度 |
|------|----------|-----------|---------|-------|
| 軸向非線性縮放 | 解析 | 低（單方向） | 高 | 易 |
| 正弦位移 | 解析 | 中 | 中 | 易 |
| 徑向推拉 | 需數值驗證 | 中（局部） | 高 | 中 |
| TPS + 控制點 | 需數值驗證 | 高 | 高（視覺化設計）| 中 |
| 串接組合 | 若各步保證則保證 | 可疊加 | 中 | 中 |

- **建議起點**：先用方法 1a（軸向縮放）確認 pipeline 正確，再換 TPS 控制點版本增加空間異質性的複雜度。
- 無論哪種方法，都應執行 Jacobian 檢查，尤其是 TPS 和徑向推拉在極端參數下容易折疊。

---

## 參考

- Sampson, P.D. and Guttorp, P. (1992). Nonparametric estimation of nonstationary spatial covariance structure. *Journal of the American Statistical Association*, **87**(417), 108–119.

