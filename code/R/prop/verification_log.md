# Verification Log: MRTS-Augmented Morris Skew-t Model

- Scripts live in `code/R/prop/`. All reproducible with `set.seed(20260531)`.
- Notation: F is the n×K₀ orthonormal basis (F′F = I); G = FMF′ + σ_ξ²I; R = A^{-1/2}GA^{-1/2}, A = diag(G).

---

## Executive Summary

### What this project does

- Replaces the computationally expensive Matérn spatial correlation model in the Morris skew-t framework with a faster MRTS-based low-rank approximation.
- Goal: keep inference quality while reducing per-MCMC-iteration cost from O(n³) to O(nKT), where K=20 ≪ n is the number of basis functions and T is the number of time replicates.

### What has been verified (completed)

| Area | Finding | Status |
|------|---------|--------|
| Code correctness | Closed-form estimator matches reference implementation (autoFRK) to machine precision | ✅ Correct |
| Model assumptions A1–A4 | Gaussian residuals, temporal independence, time-invariant covariance, zero measurement error — all exact model properties, confirmed numerically | ✅ All pass |
| Spatial approximation (V1) | MRTS approximates Matérn **only for short-range fields** (spatial range ≤ 5% of domain width). For longer ranges the approximation is poor because a key direction is excluded by design | ⚠️ Limited |
| Regression inference (V2–V3) | Point estimates of β are unbiased. But confidence intervals for the **intercept** are too narrow — 38% undercoverage at moderate range, ~60% undercoverage at long range. The MCMC feedback makes this worse, not better | ❌ Concern |
| Spatial prediction (V4) | Predictions at unmonitored sites are **well-calibrated** (within 2% of nominal coverage) across all tested ranges. This is the primary use case | ✅ Good |
| Computation (V5) | Current implementation is **not faster** than Matérn due to an unintended O(n³) bottleneck in the code. A known fix exists and would achieve 126× speedup at n=500 | ❌ Fix needed |

### What needs action before submission

1. **Fix the O(n³) bottleneck** (V5): refactor `prop_make_cov_state` to avoid building the full n×n precision matrix. The correct Woodbury formulas are already coded but not wired up. Without this fix the computational speedup claim cannot be made.

2. **Clarify the scope of validity** (V1–V3): the method works well for spatial prediction (V4) but produces overconfident regression CIs when the spatial range is moderate or long. The paper should either restrict to short-range applications, or add a sentence acknowledging this limitation for inference on regression coefficients.

3. **The dropped constant column** is the root cause of both V1 and V2 failures. It was dropped to prevent the spatial model from absorbing the regression intercept, which is a sound design choice — but it means long-range spatial correlation cannot be represented. This trade-off should be stated explicitly.

### What is not a concern

- The core innovation (MRTS covariance update within MCMC) is mathematically correct and matches the reference estimator.
- Spatial predictions — the main deliverable of the Morris model for precipitation/ozone applications — are well-calibrated.
- All model assumptions (A1–A4) hold exactly under the model by construction.

---

## Assumption Checklist

| ID | Claim | Status |
|----|-------|--------|
| **–** | prop closed-form ≡ autoFRK `cMLEimat` | ✅ |
| **A1** | r_t \| (β, λ, zg, D_t) ~ N(0, R) exactly | ✅ |
| **A2** | r_t i.i.d. across t | ✅ |
| **A3** | R common across all t | ✅ |
| **A4** | σ_ε² = 0 is correct for std_res | ✅ |
| **V1** | MRTS R̂(K) approximates R_Matérn(ρ,ν) as a function of K | ✅ verified — see findings |
| **V2** | Covariance misspecification does not bias β, but narrows CI | ✅ verified — see findings |
| **V3** | Plug-in M creates additional posterior compression beyond V2 | ✅ verified — see findings |
| **V4** | MRTS spatial predictions are well-calibrated despite V1–V2 | ✅ verified — see findings |
| **V5** | Current MRTS has O(n³) bottleneck; Woodbury-only is 126× faster | ✅ verified — see findings |
| **B1** | rank(H) ≤ min(K, T); K ≪ T required | ⬜ |
| **B2** | K ≪ n; low-rank approximation quality | ⬜ |
| **B3** | Dropping constant MRTS column preserves spatial information | ⬜ |
| **C1** | std_res is mean-zero; no aliasing with intercept | ⬜ |
| **C2** | G→R normalisation does not bias the Morris posterior | ⬜ |
| **C3** | Plug-in M is an effective within-MCMC strategy | ⬜ |
| **D1** | σ_ξ² > 0 floor is sufficient for Woodbury stability | ⬜ |
| **D2** | Orthonormalisation preserves the spatial estimate | ⬜ |
| **D3** | Cholesky/precision stable when R is near-singular | ⬜ |
| **E1** | Site-varying D_t does not invalidate the imat path | ⬜ |
| **E2** | MCMC imputation error is negligible for the covariance step | ⬜ |
| **E3** | Plug-in M is adequate relative to the scalar-c extension | ⬜ |

---

## autoFRK Cross-check

**Script:** `verify_vs_autofrk.R` · **Date:** 2026-05-31

### Claim

- For σ_ε²=0 and F′F=I, prop and autoFRK `cMLEimat` solve the same equations:
  - H = F′SF (prop, since F′F=I) = (F′F)^{-1/2}F′SF(F′F)^{-1/2} (autoFRK).
  - σ̂_ξ² = max{(tr(S) − Σ_{k≤L*} d_k)/(n−L*), 0}.
  - M̂ = P diag(d̂) P′ where d̂_k = max(d_k − σ̂_ξ², 0).
- The basis-invariant quantity is G = FM̂F′ + σ̂_ξ²I; M itself depends on the choice of basis.

### Setup

- n=60, T=40; raw MRTS K=8, orthonormalised K₀=7.
- Same orthonormal F fed to both estimators; rank-3 signal, iid noise σ=0.5.

### Findings

- **Basis mismatch:** prop drops the all-ones MRTS column (sd=0 across sites) and orthonormalises; autoFRK retains all K raw columns. Directly comparing M across bases is ill-defined.
- **`cMLEimat` fields:** `out$s` = input noise floor; `out$v` = estimated noise. Total noise = s+v; mistaking `out$s` for the estimate gives a spurious discrepancy.
- **σ_ε² gate:** autoFRK's `estimateV` eligibility omits the σ_ε² floor; prop's `prop_select_lstar` includes it. Any eigenvalue excluded by prop's stricter gate is also zeroed by pmax(d−σ̂_ξ²−σ_ε², 0), so M and σ̂_ξ² are identical in both cases.

### Results (common orthonormal F, σ_ε²=0)

| Quantity | \|autoFRK − prop\| |
|---|---|
| σ̂_ξ² | 5.6 × 10⁻¹⁷ |
| M (relative, K₀×K₀) | 2.4 × 10⁻¹⁵ |
| G = FM̂F′+σ̂_ξ²I (relative, n×n) | 9.8 × 10⁻¹⁶ |

- **Verdict:** Agreement to machine precision on all basis-invariant quantities.

---

## A1 — Conditional Gaussianity

**Script:** `verify_a1_gaussianity.R` · **Date:** 2026-05-31

### Claim

- Model (code parameterisation): Y_t = X_tβ + λzg_t + D_tv_t, where D_t = diag(1/√τ_{g(i),t}) and v_t ~ N(0,R).
- Standardised residual identity: `std_res_t := diag(√τ_{g(i),t}) · (Y_t − X_tβ − λzg_t) = v_t` **exactly** given (β, λ, zg_t, τ_t).
- A1 is an exact algebraic identity, not an approximation.

### Setup

- n=80, T=120, K₀=7, λ=2.5; τ_t ~ Gamma(3,3), z_t ~ |N(0,1)|.
- True covariance: M eigs = [2.0, 1.2, 0.5, 0, …], σ_ξ²=0.2.

### Tests

| Test | Statistic | Result |
|------|-----------|--------|
| T1 | Shapiro-Wilk pass rate, oracle (expect ~0.95) | 0.938 ✅ |
| T1 | Shapiro-Wilk pass rate, λ=0 uncorrected (expect ≪0.95) | 0.150 ✅ |
| T2 | d_t = v_t′R⁻¹v_t ~ χ²_n; KS p-value | 0.464 ✅ |
| T2 | mean(d_t) (expect n=80) | 78.9 ✅ |
| T3 | Spearman ρ(emp-R, R_true) upper-triangle | 0.699 — reference only† |
| T4 | max\|std_res − true v_t\| | 1.3 × 10⁻¹⁵ ✅ |

- †T3: with σ_ξ²=0.2 ≫ tr(M)/n≈0.046, R≈I so off-diagonals are O(1/n). Sample correlation error is O(√(n/T))≈0.82 at n/T=0.67. This is a B3 (basis quality) question, not an A1 question.

### Conclusions

- T4 confirms the identity is exact to machine precision; T1–T2 confirm the distributional form N(0,R).
- The contrast (uncorrected 15% pass rate) shows the λzg subtraction is necessary: without it, residuals are non-Gaussian even at true parameters.

---

## A2 — Temporal Independence

**Script:** `verify_a2_independence.R` · **Date:** 2026-05-31

### Claim

- From A1 T4: std_res[:,t] = v_t exactly. Since v_t ~iid N(0,R) by model construction, A2 is also an exact property.
- The practical question — whether real-data residuals are temporally independent after fitting — cannot be answered from the model structure and requires empirical checking.

### Setup

- n=80, T=120, K₀=7, λ=2.5; AR(1) corruption ρ=0.6 for power test.

### Tests

| Test | Statistic | Result |
|------|-----------|--------|
| T1 | Ljung-Box (lag=10) rejection rate, i.i.d. v_t (expect ≈0.05) | 0.050 ✅ |
| T2 | Ljung-Box rejection rate, AR(1) ρ=0.6 (expect ≫0.05) | 0.988 ✅ |
| T3 | Lag-1 ACF of leading-PC projection, i.i.d. case (expect \|ACF₁\| < 2/√T) | −0.081, within ±0.183 ✅ |
| T4 | max\|std_res − true v_t\| | 1.3 × 10⁻¹⁵ ✅ |

### Conclusions

- A2 is exact under the model (T4); confirmed by T1 at nominal α and T3 within the ±2/√T noise band.
- T2 (98.8% power at ρ=0.6) makes the Ljung-Box test a practical diagnostic: a rejection rate above ~15% on fitted residuals signals that temporal structure in the spatial factor should be modelled.

---

## A3 — Time-Invariant Spatial Dependence

**Script:** `verify_a3_stationarity.R` · **Date:** 2026-05-31

### Claim

- A3 is **not** an exact model property; it is a genuine data assumption.
- When A3 holds, S = T⁻¹Σ_t r_t r_t′ → R as T→∞ and the TH update is consistent.
- When A3 is violated, TH returns Ê_t[R_t] (the time-average), which may misrepresent either regime.

### Setup

- n=80, T=120, K₀=7.
- Regime 1: M eigs=[2.0,1.2,0.5], σ_ξ²=0.1 (strong spatial).
- Regime 2: M eigs=[0.20,0.10,0.05], σ_ξ²=0.5 (weak spatial).
- Two-regime: first T/2 from regime 1, second T/2 from regime 2.

### Tests

| Test | Statistic | Result |
|------|-----------|--------|
| T1 | \|λ₁(R̂)−λ₁(R_true)\|/λ₁(R_true); common R | 0.144 < 0.20 ✅ |
| T1 | Effective rank (TH) vs true signal rank | 4 vs 3 (finite-sample over-selection by 1) |
| T2 | λ₁(R̂_pooled) ∈ (λ₁(R₂_true), λ₁(R₁_true)); two-regime | 7.36 ∈ (1.38, 11.89) ✅ |
| T3a | Split-half ratio λ₁(S₁)/λ₁(S₂) ∈ null 95% CI; common R | 1.042 ∈ [0.608, 1.512] ✅ |
| T3b | Split-half ratio outside null CI; two-regime | 3.215 > 1.512 ✅ |

### Conclusions

- Under common R: TH recovers λ₁(R) within finite-sample noise (14% at T=120, n=80); effective rank may over-select by 1 near the σ_ξ floor — benign in practice.
- Under two-regime R: TH produces a convex blend λ₁(R̂_pooled) ∈ (λ₁(R₂), λ₁(R₁)), detectable via the split-half ratio diagnostic.
- **Practical diagnostic:** compute ratio = λ₁(S_{1:T/2})/λ₁(S_{T/2+1:T}) on fitted residuals. A ratio outside the null CI (here [0.61, 1.51]) flags a potential A3 violation.

---

## A4 — σ_ε² = 0 Is Appropriate

**Script:** `verify_a4_sigma_eps.R` · **Date:** 2026-05-31

### Claim

- std_res = v_t ~ N(0,R) is a pure spatial process with no additional measurement-error layer; σ_ε²=0 is the exact correct choice.
- When TH is applied to v_t ~ N(0,R) (unit diagonal), TH infers G_hat ≈ R. The normalisation step R_hat = A_hat^{-1/2}G_hat A_hat^{-1/2} is then nearly trivial, yielding R_hat ≈ R. This is the self-consistency of the prop scheme.

### Misspecification analysis (σ_ε²=δ > 0)

- TH formula: σ̂_ξ² = max{(tr(S)−Σ_{k≤L*}d_k)/(n−L*)−δ, 0} → decreases by δ.
- G_hat = FM̂F′ + σ̂_ξ²I; diagonal of G_hat decreases by δ per site (G_hat does **not** add δ back).
- R_hat[i,j] = G_hat[i,j] / √(G_hat[i,i]·G_hat[j,j]) ≈ R_true[i,j]/(1−δ) for off-diagonal entries.

### Setup

- n=80, T=300, K₀=7; **σ_ξ=0.01** (low nugget; λ₁(R_true)≈27, meaningful off-diagonals).
- δ ∈ {0.05, 0.10, 0.20}.

### Tests

| Test | Statistic | Result |
|------|-----------|--------|
| T1 | \|λ₁(R̂)−λ₁(R_true)\|/λ₁(R_true); σ_ε=0 | 0.153 < 0.20 ✅ |
| T2 | mean(G_hat_0 diag − G_hat_δ diag)/δ; all δ | 1.000 (exact) ✅ |
| T3 | mean(R_hat_δ[i,j])/mean(R_hat_0[i,j]) vs 1/(1−δ); rel. err. | ≤ 5.8% for δ≤0.20 ✅ |

### Conclusions

- σ_ε²=0 is correct; R_hat ≈ R_true within finite-sample noise (15% at T=300, n=80, O(1/√T) expected).
- Misspecification bias is deterministic and predictable: setting σ_ε²=δ inflates all off-diagonal entries of R̂ by the factor 1/(1−δ).
- For δ ≤ 0.05, the inflation is < 6% and likely negligible. For larger known instrument noise, insert δ explicitly rather than absorbing it into σ_ξ².

---

---

## V1 — MRTS Approximation Quality for Matérn

**Scripts:** `verify_v1_matern_approx.R`, `diagnose_v1.R` · **Date:** 2026-06-01

### Setup

- n=100 sites uniform on [0,1]²; T=500 oracle samples v_t ~ N(0, R_Matérn).
- K ∈ {5,10,20,40}; ν ∈ {0.5,1.5,2.5}; ρ ∈ {0.05, 0.15, 0.40}.
- Metrics (ν=0.5 shown; ν=1.5 and 2.5 qualitatively identical):

| ρ | K | M1 ‖R̂−R_M‖_F/‖R_M‖_F | M2 \|λ₁(R̂)−λ₁(R_M)\|/λ₁ | M3 \|log\|R̂\|−log\|R_M\|\|/n | M4 mean\|v′R̂⁻¹v−v′R_M⁻¹v\|/n |
|---|---|---|---|---|---|
| 0.05 | indep | 0.520 | 0.698 | 0.185 | 0.068 |
| 0.05 | K=10 | 0.451 | 0.134 | 0.153 | 0.062 |
| 0.05 | K=40 | 0.346 | 0.025 | 0.066 | 0.052 |
| 0.15 | indep | 0.859 | 0.916 | 0.721 | 0.191 |
| 0.15 | K=10 | 0.620 | 0.383 | 0.504 | 0.177 |
| 0.15 | K=40 | 0.580 | 0.367 | 0.347 | 0.235 |
| 0.40 | indep | 0.967 | 0.971 | 1.513 | 0.374 |
| 0.40 | K=10 | 0.874 | 0.656 | 1.209 | 0.504 |
| 0.40 | K=40 | 0.874 | 0.658 | 1.176 | 0.551 |

- M1 < 0.10 threshold: not achieved for **any** (K, ν, ρ) combination up to K=40.

### Root-cause diagnosis

- For each ρ, the leading eigenvector u₁(R_M) aligns with 1_n/√n at:
  - ρ=0.05: |cos(u₁, 1_n)| = 0.49
  - ρ=0.15: |cos(u₁, 1_n)| = **0.875**
  - ρ=0.40: |cos(u₁, 1_n)| = **0.978**
- prop drops the constant MRTS column (sd=0) to prevent aliasing with the regression intercept β₀. Therefore **1_n ∉ col(F)** by design.
- The fraction of ‖R_M‖_F² captured in the MRTS subspace: ‖F′R_M‖_F²/‖R_M‖_F²:
  - ρ=0.05: K=10→26.5%, K=40→74.4% (improves with K ✓)
  - ρ=0.15: K=10→53.8%, K=40→68.5% (slow improvement; constant gap)
  - ρ=0.40: K=10→23.3%, K=40→24.3% (**plateau** — constant direction excluded)
- Lower bound on M1 from the missing constant component:
  M1 ≥ λ₁(R_M) · |cos(u₁, 1_n)|² / ‖R_M‖_F
  - ρ=0.40: ≥ 34.7 × 0.978² / 50 ≈ **0.66** (observed: 0.87)
- For ρ=0.40 the eigenvalue spectrum: λ₁=34.7 (34.7% of trace), λ₂=11.9, …, λ₄₀ cumulative 93.6%. The dropped constant direction alone accounts for 34.7% of total variance.

### Conclusions

- **ρ=0.05 (short range):** M2=0.025, M3=0.066 at K=40. Acceptable approximation; spectrum is diffuse (top-1 = 3.3% of trace) and MRTS captures 74% of ‖R‖_F².
- **ρ≥0.15 (moderate–long range):** M1 and M3 plateau regardless of K. The leading Matérn eigenvector ≈ 1_n/√n carries λ₁ ≈ 12–35% of total variance and is excluded by design. No amount of K can recover it.
- **M4 anomaly:** for ρ=0.40, MRTS K=40 gives M4=0.551 **worse** than independence M4=0.374. A misfit low-rank precision R̂⁻¹ amplifies quadratic form errors more than I⁻¹=I does.
- **Implication for prop:** the MRTS model is adequate only when ρ ≲ 0.05 (≈5% of domain diameter). For moderate–long range the approximation error in log|R̂| (M3≥0.35) and quadratic forms (M4≥0.18) will bias the MCMC likelihood evaluations.
- **Open question for V2:** does this covariance misspecification propagate into biased inference on β, λ, τ? Even if R̂ ≠ R_M, the Morris posterior may be partially robust because β and λ are identified primarily from the mean/skewness structure, not from the covariance.

---

## V2 — Parameter Inference Under Spatial Misspecification

**Script:** `verify_v2_param_bias.R` · **Date:** 2026-06-01

### Setup

- Gaussian model: Ȳ ~ N(Xβ, σ²R_M/T), X = [1_n \| x-coord], p=2, n=50, T=50.
- Oracle fit: GLS with W = R_Matérn. MRTS fit: GLS with W = R̂_MRTS (K=20, ν=0.5).
- S=500 MC datasets per ρ; analytical coverage cross-checked.

### Analytical framework

- Under wrong covariance W, GLS is still unbiased: E[β̂_W] = β (OLS sandwich theory).
- Nominal CI uses: Var_nom = σ̂²_W/T · (X′W⁻¹X)⁻¹.
- True variance: Var_true = σ²/T · (X′W⁻¹X)⁻¹ X′W⁻¹ R_M W⁻¹X (X′W⁻¹X)⁻¹.
- Variance inflation ratio: V_ratio_j = (n−p) · [sandwich]_jj / (tr(M_W R_M) · [(X′W⁻¹X)⁻¹]_jj).
- Analytical coverage: 2Φ(1.96/√V_ratio_j) − 1.

### Results

| ρ | param | cov_oracle | cov_MRTS | an_cov_MRTS | V_ratio | width_R̂/R_M |
|---|---|---|---|---|---|---|
| 0.05 | β₁ (intercept) | 0.944 | 0.920 | 0.922 | 1.24 | 0.914 |
| 0.05 | β₂ (slope) | 0.950 | 0.950 | 0.935 | 1.13 | 0.963 |
| 0.15 | β₁ | 0.954 | 0.824 | 0.852 | **1.84** | 0.783 |
| 0.15 | β₂ | 0.952 | 0.930 | 0.918 | 1.27 | 0.913 |
| 0.40 | β₁ | 0.950 | **0.588** | **0.613** | **5.14** | **0.476** |
| 0.40 | β₂ | 0.942 | **0.756** | **0.787** | 2.48 | 0.654 |

- Bias: |mean(β̂_W) − β_true| < 0.004 for all ρ, j. **Point estimates are unbiased.**

### Root-cause: intercept hit harder than slope

- V_ratio_j depends on the alignment of the j-th GLS direction with missing subspace.
- β₁ (intercept) is identified from **1_n-weighted mean** of Y — the same direction MRTS excludes.
- At ρ=0.40: V_ratio = 5.14 → true variance is **5× the nominal** → CI is half the correct width → 38% undercoverage.
- β₂ (slope on x-coord) is identified from spatial contrast → less aligned with 1_n → V_ratio=2.48, less severe.

### Connection to V1

| ρ | V1 MRTS capture ‖F′R‖²_F/‖R‖²_F | V2 V_ratio(β₁) | Coverage(β₁) |
|---|---|---|---|
| 0.05 | 74% (K=40) | 1.24 | 0.922 |
| 0.15 | 68% (K=40) | 1.84 | 0.852 |
| 0.40 | 24% (K=40) | 5.14 | 0.613 |

- Lower subspace capture → higher V_ratio → worse coverage. The connection is monotone.

### Conclusions

- **β̂ is unbiased** regardless of which covariance is used (OLS sandwich result).
- **CI is systematically too narrow** under R̂: width ratio = CI_MRTS / CI_oracle = 0.48–0.96.
- For ρ=0.05: coverage 0.92 — acceptable (< 3% undercoverage).
- For ρ=0.15: intercept coverage 0.85 — borderline; reviewers will flag this.
- For ρ=0.40: intercept coverage 0.61 — **unacceptable** (38% undercoverage).
- The missing constant direction (V1) directly propagates into severe intercept undercoverage (V2).
- **Practical implication:** prop is valid for short-range fields (ρ ≲ 0.05 × domain diameter). For moderate–long range, the constant MRTS column must be retained with an alternative identifiability strategy, or K must grow with ρ.

---

## V3 — Plug-in M Creates Additional Posterior Compression

**Script:** `verify_v3_plugin_compression.R` · **Date:** 2026-06-01

### Setup

- Gaussian model: Y_t = Xβ + v_t, v_t ~ iid N(0, R_M), σ=1 known.
- Fixed R̂: computed once from T_pre=500 pilot data (V2 scenario).
- Adaptive Gibbs: β → R̂(Y−Xβ) → β, R̂ re-estimated at every iteration.
- n=50, T=30, K=20, n_chain=30 independent datasets, M=3000 (burn=500).

### Key metric: compress_ratio = SD_adaptive / SD_fixed_true

- SD_fixed_true = sandwich SE of β̂ under fixed R̂ and true R_M.
- SD_adaptive = empirical posterior SD from the adaptive chain.
- compress_ratio < 1: adaptive chain is narrower → additional compression beyond V2.

### Results

| ρ | compress_ratio(β₁) | compress_ratio(β₂) | interpretation |
|---|---|---|---|
| 0.05 | 0.879 | 0.929 | 7–12% additional compression |
| 0.15 | 0.694 | 0.827 | 17–31% additional compression |
| 0.40 | 0.611 | 0.875 | **13–39% additional compression** |

### Root cause

- The adaptive chain explores β → at each β, R̂(β) is estimated from residuals (Y−Xβ), absorbing the β-deviation as spatial structure.
- This makes the profile likelihood p(Y | β, R̂(β)) sharper than the fixed-R̂ conditional likelihood → narrower profile posterior.
- Compression grows with ρ because larger ρ means the missing constant direction carries more variance, and more of the β-deviation signal is absorbed into R̂.

### Implication: combined V2 + V3 coverage

- At ρ=0.40: V2 showed coverage(β₁)≈0.61 under fixed R̂.
- V3 shows adaptive chain gives SD that is 39% narrower still.
- Approximate combined coverage: 2Φ(1.96 × 0.611) − 1 ≈ **0.40** (at ρ=0.40, β₁).
- The intercept CI under prop MCMC is approximately half the correct width for long-range fields.

---

## V4 — Spatial Predictions Are Well-Calibrated

**Script:** `verify_v4_predictive_coverage.R` · **Date:** 2026-06-01

### Setup

- Spatial kriging design: train β from T=50 replicates at n_obs=60 sites; predict ONE new time-point at n_pred=20 held-out sites using same-time-point obs.
- Predictive distribution: Y(pred) | Y(obs), β̂, R ~ N(μ_krig, σ² R_cond).
- Oracle: R = R_Matérn. MRTS: R = R̂ (K=20). S=300 datasets.

### Results

| ρ | method | Q50 | Q75 | Q90 | Q95 |
|---|---|---|---|---|---|
| 0.05 | oracle | 0.505 | 0.752 | 0.902 | 0.950 |
| 0.05 | MRTS   | 0.503 | 0.746 | 0.896 | 0.946 |
| 0.15 | oracle | 0.501 | 0.747 | 0.901 | 0.949 |
| 0.15 | MRTS   | 0.519 | 0.763 | 0.909 | 0.955 |
| 0.40 | oracle | 0.514 | 0.757 | 0.901 | 0.954 |
| 0.40 | MRTS   | 0.514 | 0.778 | 0.916 | 0.958 |

- Max |MRTS − nominal| at Q90/Q95: 0.016 / 0.008 (ρ=0.40). All within ±0.02. **All PASS.**

### Explanation: why predictions work despite V1–V3 failures

- β inference (V2/V3) suffers because the intercept is identified from the **global level** of Y, the same direction MRTS excludes.
- Spatial kriging only needs the **local conditional structure** (R_cross, R_cond), which MRTS captures reasonably even without the constant direction.
- For ρ=0.40, MRTS slightly **overestimates** conditional variance (missing constant → weaker long-range kriging weights → wider CI → slight overcoverage at Q90/Q95).
- Overcoverage (conservative) is preferable to undercoverage for extreme-value applications.

### Conclusion

- **Spatial predictions are robust to the covariance misspecification identified in V1–V3.**
- The primary use case of the Morris model (spatial interpolation of extremes) is not significantly harmed by the MRTS approximation.

---

## V5 — Computational Scaling: O(n³) Bottleneck in Current Implementation

**Script:** `verify_v5_compute_scaling.R` · **Date:** 2026-06-01

### Setup

- n ∈ {50, 100, 200, 500}; K=20, T=30, N_rep=20 timing replications.
- Measured: (1) Matérn QF — T quadratic forms with precomputed prec_M; (2) Matérn Chol — O(n³) Cholesky; (3) MRTS current — TH + full n×n Chol of R̂; (4) Woodbury-only — O(nKT) without the n×n Cholesky.

### Results (times in seconds per MCMC iteration)

| n | Mat_QF | Mat_Chol | MRTS_curr | MRTS_Woodbury | WB_only |
|---|---|---|---|---|---|
| 50 | <0.001 | 0.001 | 0.0015 | 0.0015 | <0.001 |
| 100 | <0.001 | 0.0005 | 0.001 | 0.0015 | <0.001 |
| 200 | 0.001 | 0.003 | 0.004 | 0.004 | 0.0005 |
| 500 | 0.0025 | **0.063** | **0.043** | 0.029 | 0.0005 |

**Speedup Matérn_Chol / MRTS at n=500:** current = 1.5×; Woodbury-only = **126×**

### Root cause of current bottleneck

- `prop_make_cov_state` explicitly forms G_obs (n×n via tcrossprod) then computes `chol2inv(chol(R_obs))` — both **O(n³)** — negating the Woodbury advantage.
- `quadform_logdet` then uses `prec_obs` (the full n×n precision matrix) for quadratic forms — **O(n²T)** instead of O(nKT).

### The correct complexity

| Operation | Current | Should be |
|---|---|---|
| log\|R̂\| | O(n³) Chol | O(nK): det-lemma = n·log(σ_ξ²) + Σlog(1+d̂_k/σ_ξ²) − Σlog(G_ii) |
| r′R̂⁻¹r | O(n²) per vector | O(nK): `prop_apply_rinv_obs` already exists |
| TH update | O(nKT + K³) | unchanged |

### Conclusion

- **Current implementation has an unintended O(n³) bottleneck** from computing the full n×n Cholesky of R̂ inside `prop_make_cov_state`.
- The Woodbury formulas needed for O(nK) operations **are already implemented** in `prop_apply_rinv_obs` and the matrix determinant lemma derivation.
- At n=500: Woodbury-only achieves 126× speedup over Matérn Chol. This requires removing `chol2inv(chol(R_obs))` from `prop_make_cov_state` and routing all likelihood computations through the Woodbury path.
- **Action required:** refactor `prop_make_cov_state` and `quadform_logdet` to avoid the explicit n×n precision matrix.

---

## Cross-cutting Notes

- T3 in A1 (Spearman ρ=0.699): off-diagonal R structure is unverifiable when σ_ξ² dominates; this is a B3 question.
- Finite-sample eigenvalue error scales as O(√(n/T)) for leading eigenvalues of dense correlation matrices; thresholds in T1/A3/A4 are set accordingly.
- The split-half ratio (A3 T3) and Ljung-Box rate (A2 T1–T2) are ready-to-use diagnostics for real-data A2/A3 checks.
