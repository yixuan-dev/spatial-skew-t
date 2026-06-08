# A1 Verification: r_t | (beta, lambda, zg, D_t) ~ N(0, R)?
#
# Model (code parameterisation):
#   Y_t = X_t beta + lambda * zg_t + D_t v_t,
#   zg_t[i] = z[g(i), t]  (partition half-normal, tau-independent),
#   D_t = diag(1/sqrt(taug_t)),  v_t ~ N(0, R)
#
# => std_res_t = sqrt(taug_t) * (Y_t - X_t beta - lambda * zg_t) = v_t ~ N(0,R)
#
# Tests
#   T1  Oracle residuals: site-level marginals ~ N(0,1)  [Shapiro-Wilk on T values]
#   T2  Oracle residuals: Mahalanobis distances ~ chi^2_n  [KS test]
#   T3  Oracle residuals: empirical correlation ~ R  [Frobenius error]
#   T4  Uncorrected (skew not removed): site marginals NOT N(0,1)  [should fail T1]

suppressMessages({
  library(autoFRK)
})
prop_dir <- "d:/Github/spatial-skew-t/code/R/prop"
source(file.path(prop_dir, "prop_utils.R"))
source(file.path(prop_dir, "prop_basis.R"))
source(file.path(prop_dir, "prop_covariance.R"))
source(file.path(prop_dir, "prop_modules.R"))

set.seed(20260531)

## ---- setup ----------------------------------------------------------------
n  <- 80L   # spatial locations
TT <- 120L  # replicates (large enough for per-site SW tests)
K  <- 8L    # basis dimension

loc <- cbind(runif(n), runif(n))
basis <- build_basis_matrix(s_obs = loc, k = K)
Fo <- basis$F_obs  # orthonormal: F'F = I

## ---- true covariance parameters -------------------------------------------
# rank-3 signal + noise
true_eigs <- c(2.0, 1.2, 0.5, rep(0, basis$rank_kept - 3))
true_sigma_xi_sq <- 0.2

m_upd <- list(
  rotation       = diag(basis$rank_kept),
  sigma_xi_sq    = true_sigma_xi_sq,
  lowrank_eigs   = true_eigs,
  effective_rank = sum(true_eigs > 0),
  projected_eigs = true_eigs,
  trace_s        = NA_real_,
  lstar          = 3L,
  M_hat          = diag(true_eigs),
  F_rot_obs      = Fo,
  F_rot_pred     = NULL
)
cov_state <- build_R_from_G(F_obs = Fo, M_update = m_upd)
R_true <- cov_state$r_obs   # n x n correlation matrix

## ---- true regression / skew parameters ------------------------------------
true_beta   <- c(2.0, -0.8)           # intercept + slope
true_lambda <- 2.5                     # strong positive skew

X <- cbind(1, runif(n))               # n x 2 design (fixed across t)
Xbeta <- as.vector(X %*% true_beta)   # n-vector

## ---- generate Morris model data -------------------------------------------
# nknots = 1 -> one z_t and one tau_t per replicate (same for all sites)
true_z_t   <- abs(rnorm(TT))          # z_t ~ HalfNormal  (partition-level)
true_tau_t <- rgamma(TT, 3, 3)        # tau_t ~ Gamma(3,3), E[tau]=1

# zg matrix: all sites same z_t within replicate
zg_mat   <- matrix(true_z_t,   n, TT, byrow = TRUE)
taug_mat <- matrix(true_tau_t, n, TT, byrow = TRUE)

# Draw v_t ~ N(0, R) using factor structure
true_v <- matrix(NA_real_, n, TT)
for (tt in seq_len(TT)) {
  eta <- rnorm(basis$rank_kept) * sqrt(true_eigs)  # N(0, M)
  xi  <- rnorm(n) * sqrt(true_sigma_xi_sq)
  raw <- as.vector(Fo %*% eta) + xi               # N(0, G)
  g_diag <- cov_state$a_obs_diag                  # diag(G)
  true_v[, tt] <- raw / sqrt(g_diag)              # N(0, R)
}

# Y = X beta + lambda * zg + (1/sqrt(tau)) * v
D_inv_v <- sweep(true_v, 2, 1 / sqrt(true_tau_t), "*")
Y <- matrix(Xbeta, n, TT) + true_lambda * zg_mat + D_inv_v

## ---- oracle residuals (TRUE parameters) ------------------------------------
res_oracle <- update_residuals(
  y      = Y,
  x_beta = matrix(Xbeta, n, TT),
  lambda = true_lambda,
  zg     = zg_mat,
  taug   = taug_mat
)
std_oracle <- res_oracle$std    # should equal true_v

## ---- uncorrected residuals (lambda = 0, no skew removal) ------------------
res_uncorr <- update_residuals(
  y      = Y,
  x_beta = matrix(Xbeta, n, TT),
  lambda = 0,
  zg     = zg_mat,
  taug   = taug_mat
)
std_uncorr <- res_uncorr$std

## ---- T1: site-level Shapiro-Wilk -------------------------------------------
sw_oracle <- apply(std_oracle, 1, function(r) shapiro.test(r)$p.value)
sw_uncorr <- apply(std_uncorr, 1, function(r) shapiro.test(r)$p.value)
alpha_sw  <- 0.05
pass_oracle <- mean(sw_oracle > alpha_sw)
pass_uncorr <- mean(sw_uncorr > alpha_sw)

cat("=== A1 Gaussianity verification ===\n")
cat(sprintf("n=%d  T=%d  K0=%d  lambda=%.1f\n\n", n, TT, basis$rank_kept, true_lambda))

cat("--- T1: site-level Shapiro-Wilk (each row ~ N(0,1) across T) ---\n")
cat(sprintf("Oracle   fraction p>0.05 : %.3f  (expect ~0.95)\n", pass_oracle))
cat(sprintf("Uncorr   fraction p>0.05 : %.3f  (expect <<0.95 when lambda!=0)\n\n", pass_uncorr))

## ---- T2: Mahalanobis distances ~ chi^2_n -----------------------------------
# d_t = v_t' R^{-1} v_t ~ chi^2_n under H0
mahal <- sapply(seq_len(TT), function(tt) {
  vt <- std_oracle[, tt]
  as.numeric(t(vt) %*% prop_apply_rinv_obs(cov_state, vt))
})
ks_p <- ks.test(mahal, "pchisq", df = n)$p.value
cat("--- T2: Mahalanobis d_t = v_t'R^{-1}v_t ~ chi^2_n  (KS test) ---\n")
cat(sprintf("KS p-value : %.4f  (expect >0.05)\n", ks_p))
cat(sprintf("mean(d_t)  : %.2f  (expect ~%d)\n\n", mean(mahal), n))

## ---- T3: empirical correlation vs R_true -----------------------------------
# Frobenius error is O(sqrt(n/T)) and is not a tight check at n=80, T=120.
# Use Spearman rank correlation of the upper-triangle elements instead:
# a high Spearman rho confirms the correlation *structure* (ranking) is right.
sigma_hat <- apply(std_oracle, 1, sd)
emp_cov   <- tcrossprod(std_oracle) / TT
emp_R     <- sweep(sweep(emp_cov, 1, sigma_hat, "/"), 2, sigma_hat, "/")
uptri     <- upper.tri(R_true)
sp_rho    <- cor(emp_R[uptri], R_true[uptri], method = "spearman")
frob_err  <- norm(emp_R - R_true, "F") / norm(R_true, "F")

cat("--- T3: empirical correlation vs true R ---\n")
cat(sprintf("mean(site SD)              : %.4f  (expect ~1)\n", mean(sigma_hat)))
cat(sprintf("Spearman rho (upper-tri)   : %.4f  (expect >0.90)\n", sp_rho))
cat(sprintf("relative Frobenius err     : %.4f  (O(sqrt(n/T))=%.2f, reference)\n\n",
            frob_err, sqrt(n / TT)))

## ---- T4: verify true_v == std_oracle (exact equality check) ---------------
max_diff <- max(abs(std_oracle - true_v))
cat("--- T4: std_res == true v_t (exact by construction) ---\n")
cat(sprintf("max |std_res - true_v| : %.3e  (expect <1e-12)\n\n", max_diff))

## ---- summary ---------------------------------------------------------------
cat("=== VERDICT ===\n")
t1_ok <- pass_oracle > 0.85 && pass_uncorr < 0.5
t2_ok <- ks_p > 0.05
t4_ok <- max_diff < 1e-10
cat(sprintf("T1 oracle Gaussian / uncorr non-Gaussian : %s\n",
            ifelse(t1_ok, "PASS", "FAIL")))
cat(sprintf("T2 Mahalanobis ~ chi^2_n                 : %s\n",
            ifelse(t2_ok, "PASS", "FAIL")))
cat(sprintf("T3 [reference only, see B3]              : sp_rho=%.3f\n", sp_rho))
cat(sprintf("T4 exact equality with true v_t          : %s\n",
            ifelse(t4_ok, "PASS", "FAIL")))
cat(sprintf("OVERALL                                  : %s\n",
            ifelse(all(t1_ok, t2_ok, t4_ok), "PASS", "FAIL")))
