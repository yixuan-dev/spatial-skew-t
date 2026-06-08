# A2 Verification: r_t conditionally i.i.d. across t?
#
# From A1 T4: std_res[:,t] = v_t exactly at true parameters.
# Since v_t are i.i.d. N(0,R) by model construction, A2 is an
# exact model property, not an approximation.
#
# This script confirms:
#   T1  i.i.d. v_t:  site-level Ljung-Box rejection rate ≈ α = 0.05
#   T2  AR(1) v_t:   site-level Ljung-Box rejection rate >> 0.05  (power)
#   T3  leading-PC autocorrelation: ≈ 0 for i.i.d., detectable for AR(1)
#
# NOTE: A2 is an ASSUMPTION about the data, not derivable from the model.
# For real data (precipitation, ozone), temporal independence of the
# spatial-factor residuals must be assessed separately.

suppressMessages(library(autoFRK))
prop_dir <- "d:/Github/spatial-skew-t/code/R/prop"
source(file.path(prop_dir, "prop_utils.R"))
source(file.path(prop_dir, "prop_basis.R"))
source(file.path(prop_dir, "prop_covariance.R"))
source(file.path(prop_dir, "prop_modules.R"))

set.seed(20260531)

## ---- shared setup (same as A1) --------------------------------------------
n  <- 80L
TT <- 120L
K  <- 8L
lb_lag   <- 10L    # Ljung-Box lag
alpha_lb <- 0.05

loc   <- cbind(runif(n), runif(n))
basis <- build_basis_matrix(s_obs = loc, k = K)
Fo    <- basis$F_obs

true_eigs       <- c(2.0, 1.2, 0.5, rep(0, basis$rank_kept - 3))
true_sigma_xi   <- 0.2
m_upd <- list(
  rotation       = diag(basis$rank_kept),
  sigma_xi_sq    = true_sigma_xi,
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

true_beta   <- c(2.0, -0.8)
true_lambda <- 2.5
X     <- cbind(1, runif(n))
Xbeta <- as.vector(X %*% true_beta)

true_z_t   <- abs(rnorm(TT))
true_tau_t <- rgamma(TT, 3, 3)
zg_mat     <- matrix(true_z_t,   n, TT, byrow = TRUE)
taug_mat   <- matrix(true_tau_t, n, TT, byrow = TRUE)

## helper: draw T iid replicates from N(0,R) via factor model ----------------
draw_iid_v <- function(TT) {
  v <- matrix(NA_real_, n, TT)
  for (tt in seq_len(TT)) {
    eta   <- rnorm(basis$rank_kept) * sqrt(true_eigs)
    xi    <- rnorm(n) * sqrt(true_sigma_xi)
    raw   <- as.vector(Fo %*% eta) + xi
    v[, tt] <- raw / sqrt(cov_state$a_obs_diag)
  }
  v
}

## ---- T1: i.i.d. case -------------------------------------------------------
true_v_iid <- draw_iid_v(TT)
D_inv_v    <- sweep(true_v_iid, 2, 1 / sqrt(true_tau_t), "*")
Y_iid      <- matrix(Xbeta, n, TT) + true_lambda * zg_mat + D_inv_v

res_iid <- update_residuals(
  y = Y_iid, x_beta = matrix(Xbeta, n, TT),
  lambda = true_lambda, zg = zg_mat, taug = taug_mat
)

lb_pvals_iid <- apply(res_iid$std, 1, function(r)
  Box.test(r, lag = lb_lag, type = "Ljung-Box")$p.value
)
reject_iid <- mean(lb_pvals_iid < alpha_lb)

## ---- T2: AR(1) case (rho = 0.6) -------------------------------------------
rho <- 0.6
true_v_ar1 <- matrix(NA_real_, n, TT)
true_v_ar1[, 1] <- draw_iid_v(1L)
for (tt in 2:TT) {
  true_v_ar1[, tt] <- rho * true_v_ar1[, tt - 1] +
                      sqrt(1 - rho^2) * draw_iid_v(1L)
}

D_inv_v_ar1 <- sweep(true_v_ar1, 2, 1 / sqrt(true_tau_t), "*")
Y_ar1       <- matrix(Xbeta, n, TT) + true_lambda * zg_mat + D_inv_v_ar1

res_ar1 <- update_residuals(
  y = Y_ar1, x_beta = matrix(Xbeta, n, TT),
  lambda = true_lambda, zg = zg_mat, taug = taug_mat
)

lb_pvals_ar1 <- apply(res_ar1$std, 1, function(r)
  Box.test(r, lag = lb_lag, type = "Ljung-Box")$p.value
)
reject_ar1 <- mean(lb_pvals_ar1 < alpha_lb)

## ---- T3: autocorrelation in leading PC ------------------------------------
# Project oracle residuals onto 1st basis vector to get a scalar time series
pc1_iid <- as.vector(t(Fo[, 1, drop = FALSE]) %*% res_iid$std)
pc1_ar1 <- as.vector(t(Fo[, 1, drop = FALSE]) %*% res_ar1$std)

ac1_iid <- acf(pc1_iid, lag.max = 1, plot = FALSE)$acf[2]
ac1_ar1 <- acf(pc1_ar1, lag.max = 1, plot = FALSE)$acf[2]
ac_se   <- 1 / sqrt(TT)   # approximate SE for lag-1 ACF

## ---- T4: exact equality v_t == std_res (sanity, from A1) ------------------
max_diff_iid <- max(abs(res_iid$std - true_v_iid))

## ---- output ----------------------------------------------------------------
cat("=== A2 Temporal independence verification ===\n")
cat(sprintf("n=%d  T=%d  K0=%d  lb_lag=%d  AR(1) rho=%.1f\n\n",
            n, TT, basis$rank_kept, lb_lag, rho))

cat("--- T1: i.i.d. v_t  (Ljung-Box rejection rate, expect ~0.05) ---\n")
cat(sprintf("rejection rate : %.3f  (alpha=%.2f)\n\n", reject_iid, alpha_lb))

cat("--- T2: AR(1) v_t  (Ljung-Box rejection rate, expect >>0.05) ---\n")
cat(sprintf("rejection rate : %.3f  (rho=%.1f)\n\n", reject_ar1, rho))

cat("--- T3: lag-1 ACF of leading-PC projection ---\n")
cat(sprintf("i.i.d. case : ACF1 = %+.4f  (|ACF1| < 2*SE=%.3f ?  %s)\n",
            ac1_iid, 2 * ac_se, ifelse(abs(ac1_iid) < 2 * ac_se, "yes", "no")))
cat(sprintf("AR(1) case  : ACF1 = %+.4f  (rho_true=%.2f)\n\n", ac1_ar1, rho))

cat("--- T4: exact equality std_res == true v_t (i.i.d. case) ---\n")
cat(sprintf("max |std_res - true_v| : %.3e  (expect <1e-12)\n\n", max_diff_iid))

## ---- verdict ---------------------------------------------------------------
cat("=== VERDICT ===\n")
t1_ok <- reject_iid < 0.10         # LB rate near nominal alpha
t2_ok <- reject_ar1 > 0.50         # good power against AR(1)
t3_ok <- abs(ac1_iid) < 2 * ac_se  # i.i.d. lag-1 ACF within noise band
t4_ok <- max_diff_iid < 1e-10

cat(sprintf("T1 i.i.d. rejection rate ≈ alpha   : %s  (%.3f)\n",
            ifelse(t1_ok, "PASS", "FAIL"), reject_iid))
cat(sprintf("T2 AR(1) power > 0.50              : %s  (%.3f)\n",
            ifelse(t2_ok, "PASS", "FAIL"), reject_ar1))
cat(sprintf("T3 leading-PC ACF1 within 2*SE     : %s  (%+.4f, SE=%.3f)\n",
            ifelse(t3_ok, "PASS", "FAIL"), ac1_iid, ac_se))
cat(sprintf("T4 exact std_res == true v_t       : %s\n",
            ifelse(t4_ok, "PASS", "FAIL")))
cat(sprintf("OVERALL                            : %s\n",
            ifelse(all(t1_ok, t2_ok, t3_ok, t4_ok), "PASS", "FAIL")))
