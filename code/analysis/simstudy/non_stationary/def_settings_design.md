# `simdata_def.RData` 實驗設計文件

> 對應實作：[setup_def.R](../setup_def.R)
> 理論背景：[sampson_guttorp_1992.md](sampson_guttorp_1992.md)

---

## 共同參數（全 6 個 settings 共用）

| 參數 | 值 | 說明 |
|------|-----|------|
| `beta.t` | `c(10, 0, 0)` | 迴歸係數（intercept + 兩個座標斜率） |
| `gamma.t` | `0.9` | 邊際相關縮放 |
| `rho.t` | `1` | F-space 中的指數相關 range |
| `tau.alpha.t` | `3` | τ 的 Gamma 形狀參數（除 2 後為內部參數） |
| `tau.beta.t` | `8` | τ 的 Gamma 率參數（除 2 後為內部參數） |
| `lambda` | `3` | 偏態參數（skew-t-1） |
| `dist` | `"t"` | 分布族（t 分布） |
| `nknots` | `1` | 空間分割 knot 數 |
| `phi.z = phi.w = phi.tau` | `0` | 時間 AR 係數（i.i.d. in time） |
| **站點** | 144 點，`runif([0,10]²)`，`set.seed(20)` | 與 `simdata.RData` 相同網格 |
| `nt` | `50` | 時間步數 |
| `nsets` | `50` | 每個 setting 的重複資料集數 |
| `ntest` | `44` | 測試站點數（用於下游評分） |

協方差函數統一為**指數族**（`nu = 0.5` 的 Matérn 退化形式）：

$$C(s_i, s_j) = \gamma \cdot \exp\!\left(-\frac{\|D(s_i) - D(s_j)\|}{\rho}\right)$$

在 F-space 中以歐氏距離計算，`diag(C) = 1`。

---

## Settings 1–3：線性變形（幾何各向異性）

由 `CorFxDef()` 實作。變形映射為線性：

$$D_{\text{lin}}(s) = A s, \quad A = R(\theta) \cdot \text{diag}(1,\, \text{ratio})$$

$$R(\theta) = \begin{pmatrix} \cos\theta & \sin\theta \\ -\sin\theta & \cos\theta \end{pmatrix}$$

線性映射使 $\|D(s_i) - D(s_j)\| = \|A(s_i - s_j)\|$ 仍只依賴 $s_i - s_j$，故協方差**在 G-space 仍為平穩**（各向異性但平穩）。

| Setting | `deform_type` | θ | ratio | 說明 |
|---------|------------|--------|-------|------|
| 1 | `linear` | 0 | 1.00 | 等向基準（isotropic） |
| 2 | `linear` | π/4 | 0.50 | 中度各向異性，主軸 45° |
| 3 | `linear` | π/6 | 0.25 | 強度各向異性，主軸 30° |

Setting 1（θ=0, ratio=1）退化為純指數等向協方差，是其他所有 settings 的比較基準。

---

## Settings 4–6：非線性變形（真正的非平穩協方差）

由 `CorFxNL()` 實作。非線性 $D$ 使 $\|D(s_i) - D(s_j)\|$ 依賴 $s_i, s_j$ 各自的位置而非只依賴差值，協方差在 G-space 為**非平穩**。

RNG seed：`set.seed(setting * 1000 + set)`，保持與 settings 1–3 一致的 seed 結構。

---

### Setting 4：軸向指數縮放（`nl_axial`）

$$D_4(s_1, s_2) = \bigl(s_1,\; e^{\alpha s_1} \cdot s_2\bigr), \quad \alpha = 0.15$$

**Jacobian**：

$$\nabla D_4 = \begin{pmatrix} 1 & 0 \\ \alpha e^{\alpha s_1} s_2 & e^{\alpha s_1} \end{pmatrix}, \quad \det(\nabla D_4) = e^{\alpha s_1} > 0 \;\text{（解析保證雙射）}$$

**F-space 幾何效果**：$s_2$ 方向的座標被 $e^{0.15 s_1}$ 放大，$s_1$ 越大放大越多。
- 在 F-space 中，右側（$s_1$ 大）的站點在 $s_2$ 方向被拉遠 → 右側站點的 $s_2$ 方向 F-space 距離大 → G-space 中右側的 $s_2$ 方向相關衰減**更慢**（相關更強）。
- $s_1$ 方向的相關不受影響（F-space $s_1$ 座標與 G-space 相同）。
- 非平穩強度：`exp(0.15×10) ≈ 4.5`——右側 $s_2$ 方向的有效 range 約為左側的 4.5 倍。

```r
D_fns[[4]] <- function(s) cbind(s[, 1], exp(0.15 * s[, 1]) * s[, 2])
```

---

### Setting 5：正弦波橫向位移（`nl_sine`）

$$D_5(s_1, s_2) = \bigl(s_1 + a\sin(b s_2),\; s_2\bigr), \quad a = 1.5,\; b = 0.4$$

**Jacobian**：

$$\nabla D_5 = \begin{pmatrix} 1 & ab\cos(b s_2) \\ 0 & 1 \end{pmatrix}, \quad \det(\nabla D_5) = 1 > 0 \;\text{（無條件雙射）}$$

**F-space 幾何效果**：$s_1$ 座標依 $s_2$ 做正弦位移，使等相關輪廓在 G-space 中呈**波浪狀**扭曲。
- 在 $s_2 \approx \pi/(2b) \approx 3.9$ 處位移最大（$+1.5$），$s_2 \approx 3\pi/(2b) \approx 11.8$ 處反向最大（但超出域，故域內週期 ≈ 15.7，比域寬 10 大）。
- 相關在 F-space 中仍等向，但在 G-space 中，沿 $s_1$ 方向距離相同的兩點因 $s_2$ 不同而有不同的 F-space 距離 → 非平穩。
- 位移幅度 $a=1.5$ 佔域寬 15%，提供明顯但不極端的扭曲。

```r
D_fns[[5]] <- function(s) cbind(s[, 1] + 1.5 * sin(0.4 * s[, 2]), s[, 2])
```

---

### Setting 6：旋轉 ∘ 軸向縮放組合（`nl_composed`）

$$D_6 = D_{\text{ax}} \circ D_{\text{rot}}, \quad D_{\text{rot}} = R(\pi/6),\; D_{\text{ax}}(u_1, u_2) = (u_1,\; e^{0.10 u_1} u_2)$$

展開：

$$D_6(s) = D_{\text{ax}}\!\left(\cos\tfrac{\pi}{6}\, s_1 - \sin\tfrac{\pi}{6}\, s_2,\; \sin\tfrac{\pi}{6}\, s_1 + \cos\tfrac{\pi}{6}\, s_2\right)$$

**Jacobian**（鏈式法則）：

$$\det(\nabla D_6) = \det(\nabla D_{\text{ax}})\big|_{D_{\text{rot}}(s)} \cdot \det(\nabla D_{\text{rot}}) = e^{0.10 u_1} \cdot 1 > 0 \;\text{（解析保證）}$$

其中 $u_1 = \cos(\pi/6)\, s_1 - \sin(\pi/6)\, s_2$。

**F-space 幾何效果**：
- $D_{\text{rot}}$（π/6 旋轉）先改變非平穩梯度的主軸方向，使梯度不再沿座標軸而是沿 30° 方向。
- $D_{\text{ax}}$（α=0.10 軸向縮放）再在旋轉後的座標系中施加指數拉伸。
- 兩種效果疊加：旋轉後 $s_2$ 方向的相關梯度轉 30°，疊加軸向的位置依賴 range。

```r
D_rot  <- function(s) {
  th <- pi / 6
  cbind(cos(th)*s[,1] - sin(th)*s[,2], sin(th)*s[,1] + cos(th)*s[,2])
}
D_fns[[6]] <- function(s) cbind(s[,1], exp(0.10 * s[,1]) * s[,2])[
  # 實際 closure：先旋轉再縮放
]
D_fns[[6]] <- function(s) {
  u <- D_rot(s)
  cbind(u[, 1], exp(0.10 * u[, 1]) * u[, 2])
}
```

---

## 三個 NL Settings 的非平穩強度比較

| Setting | 最大 Jacobian 比值 | 非平穩方向 | 效果類型 |
|---------|------------------|-----------|---------|
| 4 | $e^{0.15 \times 10} \approx 4.5$ | $s_1$ 梯度（單調） | 全域漸進梯度 |
| 5 | 1（det 恆等於 1） | $s_2$ 週期（正弦） | 週期性扭曲 |
| 6 | $e^{0.10 \times u_{1,\max}} \approx 2$–$3$ | 旋轉 30° 後的梯度 | 複合（梯度 + 旋轉） |

Setting 4 的非平穩最強（Jacobian 變化幅度最大），setting 5 的 Jacobian 恆為 1（空間「面積」不變，但相關輪廓形狀扭曲），setting 6 居中但方向最難被平穩模型捕捉。

---

## `simdata_def.RData` 的物件結構

```
y             [144 × 50 × 50 × 6]   觀測陣列 [ns, nt, nsets, nsettings]
tau.t         list[[6]]             各 setting 的 τ：[nknots=1, nt=50, nsets=50]
z.t           list[[6]]             各 setting 的 z：[nknots=1, nt=50, nsets=50]
knots.t       list[[6]]             各 setting 的 knot：[1, 50, 2, 50]
settings.def  data.frame [6 × 7]   columns: setting, dist, nknots, lambda,
                                             theta, ratio, deform_type
f_coords_list list[[6]]            F-space 座標 [144 × 2]；settings 1–3 為 NULL
ns, nt, s, nsets, ntest, x         空間 / 時間 metadata
```

`settings.def` 內容：

| setting | dist | nknots | lambda | theta | ratio | deform_type |
|---------|------|--------|--------|-------|-------|-------------|
| 1 | t | 1 | 3 | 0 | 1.00 | linear |
| 2 | t | 1 | 3 | π/4 | 0.50 | linear |
| 3 | t | 1 | 3 | π/6 | 0.25 | linear |
| 4 | t | 1 | 3 | NA | NA | nl_axial |
| 5 | t | 1 | 3 | NA | NA | nl_sine |
| 6 | t | 1 | 3 | NA | NA | nl_composed |

---

## 實作說明

- `CorFxNL(f_coords, gamma, rho)`（`code/R/ar2/auxfunctions.R`）：接受預先計算的 F-space 座標，直接計算歐氏距離矩陣後建指數協方差。
- `rpotspatTS(..., cov.type = "nl_deformed", f_coords = f_coords)`：新增的 `nl_deformed` 分支，將 `f_coords` 傳入 `CorFxNL()`；`f_coords` 在 setting 迴圈外預先計算（`D_fns[[setting]](s)`），所有 `nsets` 重複共用同一個 F-space 座標。
- Jacobian 驗證在主迴圈前執行，`stopifnot(all(jac_det(...) > 0))` 保證雙射性。

---

## 參考

- Sampson, P.D. and Guttorp, P. (1992). Nonparametric estimation of nonstationary spatial covariance structure. *JASA*, **87**(417), 108–119.
- [setup_def.R](../setup_def.R) — 生成腳本
- [sampson_guttorp_1992.md](sampson_guttorp_1992.md) — D 函數設計理論
