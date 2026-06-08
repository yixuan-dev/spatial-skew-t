# V2: Parameter inference under spatial misspecification (Gaussian model)
#
# Model: Y_t = X beta + sigma * v_t,  v_t ~ iid N(0, R_Matern),  t=1..T
# Fitted: GLS using W = R_Matern (oracle) or W = R_hat (MRTS)
#
# For T iid replicates, sufficient statistic Ybar ~ N(X beta, sigma^2 R / T).
# GLS under W:
#   beta_hat   = (X'W^{-1}X)^{-1} X'W^{-1} Ybar
#   sigma2_hat = T/(n-p) * Ybar' M_W Ybar   [M_W = W^{-1} - hat matrix]
#   Var_beta   = sigma2_hat/T * (X'W^{-1}X)^{-1}
#   95% CI     = beta_hat +/- 1.96 * sqrt(diag(Var_beta))
#
# Analytical coverage (no MC noise):
#   V_ratio_j = Var_true(beta_hat_j) / E[Var_nominal(beta_hat_j)]
#   coverage  = 2*Phi(1.96 / sqrt(V_ratio_j)) - 1
# where Var_true uses the sandwich matrix under R_Matern.
#
# MC coverage (S=500 datasets): empirical check of analytical result.
#
# Metrics:
#   M1  coverage rate for beta_1 (intercept) — affected by missing 1_n direction
#   M2  coverage rate for beta_2 (slope on x-coord)
#   M3  CI width ratio: CI_R_hat / CI_R_Matern (efficiency loss)
#   M4  mean absolute bias of beta_hat under R_hat vs oracle
#
# Expected from V1:
#   rho=0.05: V1 showed adequate approximation   -> coverage ~0.95 expected
#   rho=0.15: V1 showed moderate error           -> some undercoverage expected
#   rho=0.40: V1 showed catastrophic M3 error   -> severe undercoverage expected

suppressMessages(library(autoFRK))
prop_dir <- "d:/Github/spatial-skew-t/code/R/prop"
for (f in c("prop_utils.R","prop_basis.R","prop_covariance.R","prop_modules.R"))
  source(file.path(prop_dir, f))

set.seed(20260601)

## ---- fixed design ---------------------------------------------------------
n      <- 50L
T_rep  <- 50L       # replicates per dataset
S      <- 500L      # Monte Carlo datasets
K      <- 20L       # MRTS basis size (realistic)
nu     <- 0.5       # Matern smoothness
rhos   <- c(0.05, 0.15, 0.40)
T_pre  <- 500L      # pilot samples for estimating R_hat

beta_true  <- c(2.0, 1.5)
sigma_true <- 1.0
p          <- length(beta_true)

loc <- cbind(runif(n), runif(n))
X   <- cbind(1, loc[, 1])          # intercept + x-coordinate

## ---- Matern correlation ---------------------------------------------------
matern_exp <- function(d, rho) exp(-d / rho)
d_mat <- as.matrix(dist(loc))

## ---- helpers: GLS under covariance W -------------------------------------
# Returns list(beta_hat, sigma2_hat, var_beta, ci_lo, ci_hi)
gls_fit <- function(Ybar, prec_W, T_rep) {
  XtWX    <- crossprod(X, prec_W %*% X)          # p x p
  XtWy    <- crossprod(X, prec_W %*% Ybar)       # p-vector
  beta_h  <- solve(XtWX, XtWy)                   # p-vector
  resid   <- Ybar - X %*% beta_h
  # sigma2 per-obs (accounting for T_rep replicates in Ybar)
  sigma2_h <- T_rep / (n - p) * as.numeric(
    crossprod(resid, prec_W %*% resid))
  var_b <- sigma2_h / T_rep * solve(XtWX)
  se    <- sqrt(pmax(diag(var_b), 0))
  list(beta = beta_h, sigma2 = sigma2_h,
       ci_lo = beta_h - 1.96 * se,
       ci_hi = beta_h + 1.96 * se,
       se = se)
}

## ---- analytical coverage --------------------------------------------------
# V_ratio_j = Var_true(beta_hat_j) / E[SE^2(beta_hat_j)]
# = n * [sandwich_jj] / (tr(M_W R_M) * [(X'W^{-1}X)^{-1}]_jj * (n-p)/n)
# where M_W = W^{-1} - W^{-1}X(X'W^{-1}X)^{-1}X'W^{-1}
analytical_coverage <- function(prec_W, R_M) {
  XtWX    <- crossprod(X, prec_W %*% X)
  XtWXinv <- solve(XtWX)
  # sandwich: actual Var_beta * T = (X'W^{-1}X)^{-1} X'W^{-1} R_M W^{-1}X (X'W^{-1}X)^{-1}
  WX       <- prec_W %*% X
  sandwich <- XtWXinv %*% crossprod(WX, R_M %*% WX) %*% XtWXinv
  # annihilator M_W = W^{-1} - W^{-1}X(X'W^{-1}X)^{-1}X'W^{-1}
  M_W     <- prec_W - WX %*% XtWXinv %*% t(WX)
  tr_M_R  <- sum(M_W * R_M)   # tr(M_W R_M) = sum(M_W * R_M) (Frobenius inner)
  # E[sigma2_hat] = sigma_true^2 * tr(M_W R_M) / (n-p)
  # Var_nominal = E[sigma2_hat] / T * (X'W^{-1}X)^{-1}
  #             = sigma_true^2 * tr(M_W R_M) / ((n-p)*T) * (X'W^{-1}X)^{-1}
  V_ratio  <- (n - p) * diag(sandwich) / (tr_M_R * diag(XtWXinv))
  coverage <- 2 * pnorm(1.96 / sqrt(V_ratio)) - 1
  list(V_ratio = V_ratio, coverage = coverage,
       tr_M_R = tr_M_R)
}

## ---- main loop ------------------------------------------------------------
results <- vector("list", length(rhos))

for (ri in seq_along(rhos)) {
  rho <- rhos[ri]
  cat(sprintf("\nrho=%.2f ...\n", rho))

  R_M      <- matern_exp(d_mat, rho); diag(R_M) <- 1
  L_M      <- t(prop_stable_chol(R_M))     # lower Cholesky of R_M
  prec_M   <- chol2inv(t(L_M))             # R_M^{-1}

  ## --- estimate R_hat from pilot data ------------------------------------
  v_pilot  <- L_M %*% matrix(rnorm(n * T_pre), n, T_pre)
  m_upd    <- update_M_closed_form(
    F_obs = build_basis_matrix(loc, k = K)$F_obs,
    residuals_std = v_pilot, sigma_eps_sq = 0, sigma_floor = 1e-6)
  bs       <- build_basis_matrix(loc, k = K)
  cs       <- build_R_from_G(bs$F_obs, m_upd)
  R_hat    <- cs$r_obs
  prec_hat <- chol2inv(prop_stable_chol(R_hat))

  ## --- analytical coverage ------------------------------------------------
  ac_M   <- analytical_coverage(prec_M,   R_M)  # oracle (should be ~1.00)
  ac_hat <- analytical_coverage(prec_hat, R_M)  # MRTS (the key result)

  ## --- Monte Carlo coverage -----------------------------------------------
  cov_M   <- matrix(0L, S, p)
  cov_hat <- matrix(0L, S, p)
  bias_hat <- matrix(0, S, p)
  width_ratio <- matrix(0, S, p)

  for (ss in seq_len(S)) {
    # Draw Ybar ~ N(X beta_true, sigma_true^2 R_M / T_rep)
    Ybar <- X %*% beta_true +
            sigma_true / sqrt(T_rep) * as.vector(L_M %*% rnorm(n))

    fit_M   <- gls_fit(Ybar, prec_M,   T_rep)
    fit_hat <- gls_fit(Ybar, prec_hat, T_rep)

    cov_M[ss, ]    <- (fit_M$ci_lo   <= beta_true) & (beta_true <= fit_M$ci_hi)
    cov_hat[ss, ]  <- (fit_hat$ci_lo <= beta_true) & (beta_true <= fit_hat$ci_hi)
    bias_hat[ss, ] <- fit_hat$beta - beta_true
    width_ratio[ss, ] <- fit_hat$se / fit_M$se
  }

  results[[ri]] <- list(
    rho          = rho,
    mc_cov_M     = colMeans(cov_M),
    mc_cov_hat   = colMeans(cov_hat),
    an_cov_M     = ac_M$coverage,
    an_cov_hat   = ac_hat$coverage,
    V_ratio      = ac_hat$V_ratio,
    bias_hat     = colMeans(bias_hat),
    width_ratio  = colMeans(width_ratio)
  )
}

## ---- output ---------------------------------------------------------------
cat("\n=== V2: GLS coverage under spatial misspecification ===\n")
cat(sprintf("n=%d  T=%d  S=%d  K=%d  nu=%.1f\n\n", n, T_rep, S, K, nu))

cat(sprintf("%-6s %-6s | %-14s %-14s | %-14s %-14s | %-10s %-10s\n",
            "rho", "param",
            "cov_oracle", "cov_MRTS",
            "an_cov_orc", "an_cov_MRTS",
            "V_ratio", "width_R"))
cat(strrep("-", 90), "\n")

for (ri in seq_along(rhos)) {
  r <- results[[ri]]
  for (j in 1:p) {
    pname <- if (j==1) "beta_1" else "beta_2"
    cat(sprintf("%-6.2f %-6s | %-14.3f %-14.3f | %-14.3f %-14.3f | %-10.3f %-10.3f\n",
                r$rho, pname,
                r$mc_cov_M[j], r$mc_cov_hat[j],
                r$an_cov_M[j], r$an_cov_hat[j],
                r$V_ratio[j], r$width_ratio[j]))
  }
  cat(sprintf("  bias(beta_1)=%.4f  bias(beta_2)=%.4f\n",
              r$bias_hat[1], r$bias_hat[2]))
  cat("\n")
}

## ---- verdict --------------------------------------------------------------
cat("=== VERDICT ===\n")
for (ri in seq_along(rhos)) {
  r <- results[[ri]]
  drop1 <- r$an_cov_hat[1]   # intercept coverage under MRTS
  drop2 <- r$an_cov_hat[2]   # slope coverage under MRTS
  flag <- if (min(drop1, drop2) < 0.85) "WARN" else "OK"
  cat(sprintf("rho=%.2f  intercept_cov=%.3f  slope_cov=%.3f  [%s]\n",
              r$rho, drop1, drop2, flag))
}
