# V4: Predictive coverage at held-out sites under spatial misspecification
#
# Setup (proper spatial kriging):
#   Train:  Y_1,...,Y_T at n_obs sites  -> estimate beta
#   Test:   one new time point Y_{T+1} observed at ALL sites;
#           hold out pred sites to check coverage
#
# Predictive distribution for Y_{T+1}(pred) | Y_{T+1}(obs), beta_hat, R:
#   mean  =  X_pred beta_hat  +  R_cross prec_oo  (Y_{T+1}(obs) - X_obs beta_hat)
#   var   =  sigma^2 * (R_pp - R_cross prec_oo R_cross')   [= sigma^2 R_cond]
#
# Prediction depends on the within-time-point spatial covariance, not on T.
# T only affects beta estimation quality.
#
# Coverage metric: fraction of Y_{T+1}(pred) <= Q_q of predictive dist
# Oracle: R = R_Matern (correctly specified)
# MRTS:   R = R_hat   (TH estimate, K=20)
#
# Expected: oracle coverage ~ nominal q; MRTS deviates based on V1 error.

suppressMessages(library(autoFRK))
prop_dir <- "d:/Github/spatial-skew-t/code/R/prop"
for (f in c("prop_utils.R","prop_basis.R","prop_covariance.R","prop_modules.R"))
  source(file.path(prop_dir, f))

set.seed(20260601)

n_obs  <- 60L; n_pred <- 20L; T_rep <- 50L; K <- 20L
S      <- 300L; sigma  <- 1.0; nu <- 0.5
rhos   <- c(0.05, 0.15, 0.40)
quants <- c(0.50, 0.75, 0.90, 0.95)
p      <- 2L

loc_all  <- cbind(runif(n_obs + n_pred), runif(n_obs + n_pred))
loc_obs  <- loc_all[seq_len(n_obs), , drop = FALSE]
loc_pred <- loc_all[(n_obs + 1L):(n_obs + n_pred), , drop = FALSE]
X_obs    <- cbind(1, loc_obs[, 1])
X_pred   <- cbind(1, loc_pred[, 1])
beta_true <- c(2.0, 1.5)

matern_exp <- function(d, rho) exp(-d / rho)
d_all  <- as.matrix(dist(loc_all))

basis  <- build_basis_matrix(s_obs = loc_obs, s_pred = loc_pred, k = K)
Fo     <- basis$F_obs
Fp     <- basis$F_pred

results <- vector("list", length(rhos))

for (ri in seq_along(rhos)) {
  rho  <- rhos[ri]
  R_all <- matern_exp(d_all, rho); diag(R_all) <- 1
  R_oo  <- R_all[seq_len(n_obs),   seq_len(n_obs),   drop = FALSE]
  R_op  <- R_all[seq_len(n_obs),   (n_obs+1L):(n_obs+n_pred), drop = FALSE]
  R_pp  <- R_all[(n_obs+1L):(n_obs+n_pred), (n_obs+1L):(n_obs+n_pred), drop = FALSE]

  L_all    <- t(prop_stable_chol(R_all))
  L_oo     <- t(prop_stable_chol(R_oo))
  prec_oo  <- chol2inv(t(L_oo))
  krig_wt  <- t(R_op) %*% prec_oo          # n_pred x n_obs
  R_cond   <- R_pp - krig_wt %*% R_op      # n_pred x n_pred (conditional cov)
  cond_sd_oracle <- sqrt(pmax(diag(R_cond), 0))  # per-site predictive SD

  ## estimate R_hat from T_pre=500 pilot oracle residuals
  v_pilot <- L_oo %*% matrix(rnorm(n_obs * 500L), n_obs, 500L)
  m_upd   <- update_M_closed_form(Fo, v_pilot, sigma_eps_sq = 0, sigma_floor = 1e-6)
  m_upd$F_rot_pred <- Fp %*% m_upd$rotation
  cs      <- build_R_from_G(Fo, m_upd, F_pred = Fp)
  prec_hat     <- cs$prec_obs
  R_cond_mrts  <- cs$r_pred -
    cs$r_cross_pred_obs %*% (prec_hat %*% t(cs$r_cross_pred_obs))
  cond_sd_mrts <- sqrt(pmax(diag(R_cond_mrts), 0))

  below_oracle <- matrix(0, length(quants), S)
  below_mrts   <- matrix(0, length(quants), S)

  for (ss in seq_len(S)) {
    ## train: estimate beta from T_rep obs replicates at obs sites
    V_train  <- L_oo %*% matrix(rnorm(n_obs * T_rep), n_obs, T_rep)
    Y_train  <- as.vector(X_obs %*% beta_true) + sigma * V_train
    Ybar     <- rowMeans(Y_train)
    XtWX     <- crossprod(X_obs, prec_oo %*% X_obs)
    beta_hat <- solve(XtWX, crossprod(X_obs, prec_oo %*% Ybar))

    ## test: new time point observed at ALL sites
    v_test_all <- L_all %*% rnorm(n_obs + n_pred)
    Y_test_obs  <- as.vector(X_obs  %*% beta_hat) +
                   sigma * v_test_all[seq_len(n_obs)]
    Y_test_pred <- as.vector(X_pred %*% beta_true) +
                   sigma * v_test_all[(n_obs+1L):(n_obs+n_pred)]

    ## oracle kriging at test time point
    r_obs_oracle <- Y_test_obs - as.vector(X_obs %*% beta_hat)
    pred_mean_o  <- as.vector(X_pred %*% beta_hat) +
                    as.vector(krig_wt %*% r_obs_oracle)

    ## MRTS kriging at test time point
    r_obs_mrts   <- Y_test_obs - as.vector(X_obs %*% beta_hat)
    pred_mean_h  <- as.vector(X_pred %*% beta_hat) +
                    as.vector(cs$r_cross_pred_obs %*% (prec_hat %*% r_obs_mrts))

    for (qi in seq_along(quants)) {
      q <- quants[qi]
      cutoff_o <- pred_mean_o + qnorm(q) * sigma * cond_sd_oracle
      cutoff_h <- pred_mean_h + qnorm(q) * sigma * cond_sd_mrts
      below_oracle[qi, ss] <- mean(Y_test_pred <= cutoff_o)
      below_mrts[qi, ss]   <- mean(Y_test_pred <= cutoff_h)
    }
  }

  results[[ri]] <- list(
    rho        = rho,
    cov_oracle = rowMeans(below_oracle),
    cov_mrts   = rowMeans(below_mrts)
  )
  cat(sprintf("rho=%.2f done\n", rho))
}

## output -------------------------------------------------------------------
cat("\n=== V4: Predictive coverage at held-out sites ===\n")
cat(sprintf("n_obs=%d  n_pred=%d  T=%d  S=%d  K=%d  nu=%.1f\n\n",
            n_obs, n_pred, T_rep, S, K, nu))
cat(sprintf("%-6s %-6s | %s\n", "rho", "method",
            paste(sprintf("Q%-4.0f", quants * 100), collapse = " ")))
cat(strrep("-", 54), "\n")

for (ri in seq_along(rhos)) {
  r <- results[[ri]]
  cat(sprintf("%-6.2f %-6s | %s\n", r$rho, "oracle",
              paste(sprintf("%.3f", r$cov_oracle), collapse = " ")))
  cat(sprintf("%-6.2f %-6s | %s\n", r$rho, "MRTS",
              paste(sprintf("%.3f", r$cov_mrts), collapse = " ")))
  cat("\n")
}

cat("=== VERDICT ===\n")
for (ri in seq_along(rhos)) {
  r   <- results[[ri]]
  dev90_o <- abs(r$cov_oracle[quants == 0.90] - 0.90)
  dev90_h <- abs(r$cov_mrts[quants   == 0.90] - 0.90)
  dev95_o <- abs(r$cov_oracle[quants == 0.95] - 0.95)
  dev95_h <- abs(r$cov_mrts[quants   == 0.95] - 0.95)
  flag_o <- if (max(dev90_o, dev95_o) > 0.05) "WARN" else "OK"
  flag_h <- if (max(dev90_h, dev95_h) > 0.05) "WARN" else "OK"
  cat(sprintf("rho=%.2f  oracle [%s]: Q90_err=%.3f  Q95_err=%.3f\n",
              r$rho, flag_o, dev90_o, dev95_o))
  cat(sprintf("         MRTS   [%s]: Q90_err=%.3f  Q95_err=%.3f\n",
              flag_h, dev90_h, dev95_h))
}
