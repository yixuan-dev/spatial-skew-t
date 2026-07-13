# time_block_forecast / simstudy

Simulation study for the **time-block holdout** redesign of
`tex/time_block_strategy/ar2_rethink.tex`. It tests whether AR(2)
temporal pooling improves prediction of extremes once the evaluation is
no longer structurally blind to temporal dependence.

## Why this study exists

The Morris baseline study (`code/analysis/simstudy/`) scored AR(2) and
i.i.d. skew-t models with the Brier score under a **spatial-only**
holdout and found a null result. `ar2_rethink.tex` proves this is an
artifact of the experiment, not of temporal pooling:

- **Kernel locality** (`prop:kernel-decomposition`) — the predictive at
  a spatially held-out cell `(s*, t)` conditions on the time-`t` latents
  alone, so the Brier score depends on the predictive only through the
  marginal exceedance probability at that cell.
- **Marginal invariance** (`prop:invariance`, `prop:phi-free`) — the
  standardised construction leaves every one-point marginal at `N(0,1)`
  for all `phi`, including the i.i.d. case `phi = 0`; `phi` moves only
  the latent autocovariance. The AR(2) and i.i.d. specifications thus
  induce *identical prior marginals* at `(s*, t)`, and `phi` can reach
  the Brier score only through the posterior — a residual benefit
  bounded at `O(q_t^-2)` (`prop:frozen-bound`), small enough to be
  swamped, so the lower-dimensional model wins the bias–variance
  contest.
- **Section 4** — the cure is to hold out **time blocks** and forecast
  forward from the posterior of the seam state, scoring by lead time.

This study implements that redesign.

## Pipeline (mirrors `code/analysis/simstudy/`)

1. **Data generation** — `setup.R` → `simdata.RData`
2. **Fit + forecast** — `run-settings.R` → `results/<setting>-<method>-<dataset>.RData`
3. **Post-fit pipeline**
   - `scores.R` (Stage 1) → `output/results/scores<setting>.RData`
   - `tables.R` (Stage 2) → CSV tables + `simresults<setting>.RData`
   - `plots.R`  (Stage 3) → PDF lead-time curves

`time_block_helpers.R` carries the shared CLI / filename / seed / catalog
helpers and the forecasting function. The AR(2) model code is imported
through the Morris study's loader, `../../simstudy/ar2_load.R` — this
study keeps no private copy of the backend.

## Scripts

| script                 | role                                                                   |
| ---------------------- | ---------------------------------------------------------------------- |
| `setup.R`              | generate `simdata.RData` (**implemented + run**)                       |
| `time_block_helpers.R` | CLI / filename / seed / catalogs / block geometry / `forecast_block()` |
| `run-settings.R`       | per-block prefix fit + Algorithm 1 forecast driver                     |
| `scores.R`             | Stage 1 — lead-time CRPS/Brier + energy + variogram scores             |
| `tables.R`             | Stage 2 — lead-time curve tables, relative curve, crossing lead        |
| `plots.R`              | Stage 3 — lead-time curve PDFs (SE band, h*, crossing-lead marker)     |

## Data settings (`--setting`, the data-generating axis)

Every setting is the same family — skew-t, `dist = "t"`, `K = 1` knot,
`lambda = 3` — and differs only in the AR(2) coefficient pair `phi`
shared by the latent processes `tau*, z*, w*` (eq. 2 of the note):

| setting | label    | `phi = (phi1, phi2)` | `rho(F)` | `h*(0.05)` | role                                               |
| ------- | -------- | -------------------- | -------- | ---------- | -------------------------------------------------- |
| 1       | iid      | (0.00, 0.00)         | 0        | 0          | null control — AR(2) must **not** beat i.i.d. here |
| 2       | weak     | (0.12, -0.05)        | 0.224    | 3          | small spectral radius, short decorrelation time    |
| 3       | moderate | (0.60, -0.30)        | 0.548    | 5          | the note's reference case (Remark 5)               |
| 4       | strong   | (0.80, -0.35)        | 0.592    | 6          | largest spectral radius of the four                |

`rho(F)` is the companion-matrix spectral radius (Definition 5) and
`h*(0.05)` the effective memory horizon (eq. 14); both come from
`ar2_spectral_radius()` / `ar2_memory_horizon()` in
[time_block_helpers.R](time_block_helpers.R).

All four settings are **short-memory** processes: a stationary AR(2) has
an autocorrelation function that decays geometrically, `rho(h) ~ rho(F)^h`,
so `sum_h |rho(h)| < Inf`. "Strong" here means a longer decorrelation
time, *not* long memory in the Hosking/Granger sense (hyperbolic decay
`rho(h) ~ C h^(2d-1)`, which no finite-order AR can reproduce). A genuine
long-memory data-generating process — e.g. ARFIMA(0, d, 0) — would be a
separate setting, and is the natural DGP for an AR(2)-vs-AR(1)
misspecification study, since it lies outside both model classes.

All three non-i.i.d. pairs have **complex** characteristic roots
(`phi1^2 + 4 phi2 < 0`), so their ACFs are damped sinusoids rather than
monotone decays.

With `K = 1` the Voronoi membership is constant, so `w` carries no
signal and the temporal effect lives entirely in `tau` and `z` — this
isolates the AR(2) signal exactly as Remark 8 recommends.

## Method catalog (`--methods`, the analysis axis)

| method | label                                                                                  |
| ------ | -------------------------------------------------------------------------------------- |
| 1      | Skew-t, K=1, i.i.d. in time (baseline)                                                 |
| 2      | Skew-t, K=1, AR(2) temporal (tau, z, w)                                                |
| 3      | Skew-t, K=1, AR(2) temporal, fixed membership (Remark 8 ablation; ≡ method 2 when K=1) |

## Time-block holdout geometry

Stored in `simdata.RData`:

- `nt = 200` — long record so several near-independent seams fit.
- `block_H = 15` — forecast horizon (lead times `1..H`), within the
  practical range `H ∈ [5,15]` of Definition 7.
- `block_seams = c(50, 80, 110, 140, 170)` — five expanding-window
  blocks. Block `b` fits on the contiguous prefix `y[, 1:T_o]` and
  forecasts the window `(T_o, T_o + H]`. Seams are spread across the
  record so the five seam states are near-independent.

## Forecasting (`forecast_block()`, Algorithm 1 of the note)

For each MCMC draw `m`:

1. extract the seam state `(X_{T_o-1}, X_{T_o})` from draw `m`'s imputed
   latent trajectory — **never** from the stationary distribution
   (Corollary 1: that erases the AR(2) signal);
2. recurse the AR(2) on the latent Gaussian scale with Yule-Walker
   innovation variance (eq. 5), so the marginal variance stays 1;
3. apply the copula transforms (eqs. 8–9) **after** the recursion;
4. draw the spatial field.

Averaging the predictive over draws propagates parameter *and*
seam-state uncertainty (Proposition 3). The i.i.d. baseline instead
draws each latent slice from the `N(0,1)` marginal, reproducing the
stationary predictive of Corollary 1.

## 預測流程（中文版，對應 `forecast_block()`）

對每一個 held-out block `b`（seam time $T_o$，預測範圍 $H$），每個 MCMC
後驗 draw $m = 1, \dots, M$ 都會走過下面**五個階段**，最後產出
`yhat[m, s, h]` 陣列，shape 為 $M \times n_s \times H$。完整英文版見
`tex/time_block_strategy/ar2_rethink.tex` §4.2「Forecast pipeline at a glance」。

### 階段 1：在 prefix 上擬合

- 用觀察資料 $y(\cdot, t \le T_o)$ 跑一次 MCMC（`run-settings.R` 裡呼叫 `mcmc()`）
- 得到 $M$ 個後驗 draws：$\theta^{(m)} = (\beta, \lambda, \alpha, \beta_\tau, \rho, \nu, \gamma, \phi_\tau, \phi_z)^{(m)}$
- 連帶 imputed observable latent 軌跡 $\tau^{(m)}_{k, t}, z^{(m)}_{k, t}$（$t \le T_o$）
- 144 個 sites 全部都觀察到 — holdout **純粹是時間軸**的

### 階段 2：用 forward copula 把 seam state 還原回 Gaussian scale

- 從每個 draw 自己的軌跡裡，取 $t \in \{T_o - 1, T_o\}$ 那兩欄
- 轉回 Gaussian 潛在尺度：
  - $\tau^{\star(m)}_{k, t} = \Phi^{-1}\!\big(G_{\alpha^{(m)}/2,\,\beta^{(m)}/2}(\tau^{(m)}_{k, t})\big)$
  - $z^{\star(m)}_{k, t} = \Phi^{-1}\!\big(H_{\sigma_z}(z^{(m)}_{k, t})\big),\quad \sigma_z = (\tau^{(m)}_{k, t})^{-1/2}$
- ⚠️ **不可以**從 stationary distribution $\mathcal{N}(\mathbf{0}, \Sigma)$ 抽 — 那會把 AR(2) 訊號抹掉（Corollary 1）

### 階段 3：AR(2) latent 遞迴（$h = 1, \dots, H$）

- 在 Gaussian scale 上遞迴：
  - $X^{\star(m)}_{k, T_o + h} = \phi^{(m)}_1 X^{\star(m)}_{k, T_o + h - 1} + \phi^{(m)}_2 X^{\star(m)}_{k, T_o + h - 2} + \sigma^{(m)} \xi,\quad \xi \sim \mathcal{N}(0, 1)$
- 創新標準差 $\sigma^{(m)}$ 由 **Yule–Walker** 鎖定（讓 marginal variance 維持 1）：
  - $\gamma_1 = \phi_1 / (1 - \phi_2),\quad \gamma_2 = \phi_1 \gamma_1 + \phi_2,\quad \sigma^2 = 1 - \phi_1 \gamma_1 - \phi_2 \gamma_2$
- $\tau^\star$ 和 $z^\star$ **各自獨立遞迴**，各自用各自的 $\phi^{(m)}$
- **i.i.d. baseline（method 1）跳過這步**，直接 $\tau^\star, z^\star \overset{\text{iid}}{\sim} \mathcal{N}(0, 1)$（Corollary 1 的 $h \to \infty$ limit）

### 階段 4：用 inverse copula 轉回 observable scale

- $\tau^{(m)}_{k, T_o + h} = G^{-1}_{\alpha^{(m)}/2,\,\beta^{(m)}/2}\!\big(\Phi(\tau^{\star(m)}_{k, T_o + h})\big)$
- $z^{(m)}_{k, T_o + h} = H^{-1}_{\sigma_z}\!\big(\Phi(z^{\star(m)}_{k, T_o + h})\big),\quad \sigma_z = (\tau^{(m)}_{k, T_o + h})^{-1/2}$
- $\tau$–$z$ coupling **在 draw $m$ 內部要保留**：$\sigma_z$ 用同一個 draw 預測出來的 $\tau$，不能用獨立 marginal

### 階段 5：抽空間場

- 對每個 lead $h$（時間 $t = T_o + h$）：
  - $\hat y^{(m)}(s, t) = x(s, t)^\top \beta^{(m)} + \lambda^{(m)} z^{(m)}_{g(s, t), t} + \varepsilon^{(m)}(s, t)$
  - $\varepsilon^{(m)}(\cdot, t) \sim \mathcal{N}\!\big(\mathbf{0},\,(\tau^{(m)}_{g(s, t), t})^{-1} C^{(m)}\big)$
- 實作：先取 Cholesky $L L^\top = C^{(m)}$，再抽 $\varepsilon = L \zeta$，其中 $\zeta \sim \mathcal{N}(0, \tau^{-1} I)$
- $K = 1$ 時 membership $g(s, t) \equiv 1$，一個 knot 服務所有 144 個 sites

### 輸出與不確定性傳遞

- 每個 block 產出 `yhat` shape `[M, n_s, H]`，配對 truth `y_val` shape `[n_s, H]`
- **對 $m$ 取平均**自動把 **parameter uncertainty** 與 **seam-state uncertainty** 同時積分掉（Proposition 3）
- 完全不用 plug-in posterior mean — 所以 predictive 不會 under-disperse

### 5 個 block 之間的關係

- 5 個 seams 設在 $\{50, 80, 110, 140, 170\}$，相鄰間距 30 步 ≫ $H = 15$
- 不同 block 的 forecast window 不會重疊，seam state **接近獨立**
- 這個近似獨立就是 lead-time curve 的 cross-block standard error 能成立的關鍵
  （Definition 9）：$\mathrm{SE}(h) = \frac{1}{\sqrt{B}} \mathrm{sd}_b\{\bar S_b(h)\}$
- 一個 block 給點估計，五個 block 給**誤差帶** — AR(2) 與 i.i.d. 的差異才有辦法跟 sampling noise 分開

## Scoring (Section 4.3)

- **Lead-time curve** `S_bar(h)` — univariate proper scores (CRPS, and
  Brier at a grid of high thresholds) evaluated separately at each lead
  `h`, with a standard error from the dispersion of the `B` per-block
  means (Definition 9). Expected signature of a genuine effect: the
  AR(2) curve lies below the i.i.d. baseline at short leads and
  converges to it near the memory horizon `h*` — the **crossing lead**
  is the headline quantity.
- **Joint summary** — energy score and variogram score over the
  per-site length-`H` time vector, one number per method.

## Usage

```powershell
$R = "C:\Program Files\R\R-4.5.1\bin\Rscript.exe"

# Stage 0: generate data (already done)
& $R .\setup.R

# Stage A: fit + forecast (smoke test: setting 3, dataset 1, both methods)
& $R .\run-settings.R --setting=3 "1" 1 "1:2"

# Stage 1 + 2 + 3: scores + tables + plots
& $R .\scores.R --setting=3 --datasets="1" --methods="1:2"
& $R .\tables.R --setting=3
& $R .\plots.R  --setting=3
```

## Status

- `setup.R` — implemented and executed; `simdata.RData` generated
  (144 sites × 200 times × 50 datasets × 4 settings).
- `run-settings.R`, `scores.R`, `tables.R`, `plots.R`,
  `time_block_helpers.R` — implemented to the design above; not yet
  validated by a run. Run the
  smoke test before scaling out, and check the three diagnostics of
  Remark 9 (`H=0` parity, `phi=0` collapse, long-horizon `N(0,1)`).
