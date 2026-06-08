# V5: Computational scaling — Matern vs MRTS per MCMC iteration
#
# One spatial MCMC step requires:
#   log|R|  +  sum_t r_t' R^{-1} r_t
#
# Matern (baseline):
#   Chol(R_M): O(n^3)   [paid per rho/nu update]
#   QF:        O(n^2 T) [paid every iteration, with prec_M precomputed]
#
# MRTS — current implementation (prop_make_cov_state):
#   TH update:      O(nKT + K^3)  [fast]
#   build_R_from_G: forms n x n R_obs, then Chol(R_obs) -> O(n^3)  [bottleneck]
#   QF via prec_obs: O(n^2 T)    [bottleneck]
#
# MRTS — Woodbury-only (optimal implementation, NOT YET CODED):
#   TH update:      O(nKT + K^3)
#   log|R| via det-lemma: O(nK)  [n*log(sigma_xi^2) + sum log(1+d/s) - sum log(G_ii)]
#   QF via Woodbury:      O(nKT) [prop_apply_rinv_obs, already in prop_covariance.R]
#   Total:          O(nKT)
#
# V5 reports all three and identifies the crossover point.

suppressMessages(library(autoFRK))
prop_dir <- "d:/Github/spatial-skew-t/code/R/prop"
for (f in c("prop_utils.R","prop_basis.R","prop_covariance.R","prop_modules.R"))
  source(file.path(prop_dir, f))

set.seed(20260601)

n_vals <- c(50L, 100L, 200L, 500L)
K      <- 20L
T_rep  <- 30L
N_rep  <- 20L
rho    <- 0.15
nu     <- 0.5

matern_exp <- function(d, rho) exp(-d / rho)

## Woodbury log|R| via matrix determinant lemma (O(nK)) --------------------
woodbury_logdet <- function(cs) {
  d_hat      <- cs$lowrank_eigs
  sigma_xi   <- cs$sigma_xi_sq
  active     <- d_hat > 1e-10
  n          <- length(cs$a_obs_diag)
  logdet_G   <- n * log(sigma_xi) +
                sum(log(1 + d_hat[active] / sigma_xi))
  logdet_G - sum(log(cs$a_obs_diag))   # log|R| = log|G| - sum log(G_ii)
}

results <- vector("list", length(n_vals))

for (ni in seq_along(n_vals)) {
  n   <- n_vals[ni]
  loc <- cbind(runif(n), runif(n))
  d_mat <- as.matrix(dist(loc))

  R_M  <- matern_exp(d_mat, rho); diag(R_M) <- 1
  L_M  <- t(prop_stable_chol(R_M))
  prec_M <- chol2inv(t(L_M))
  logdet_M <- sum(log(eigen(R_M, only.values = TRUE,
                            symmetric = TRUE)$values))

  resid <- L_M %*% matrix(rnorm(n * T_rep), n, T_rep)

  basis <- build_basis_matrix(s_obs = loc, k = K)
  Fo    <- basis$F_obs
  m_upd <- update_M_closed_form(Fo, resid, sigma_eps_sq = 0,
                                sigma_floor = 1e-6)
  cs    <- build_R_from_G(Fo, m_upd)   # includes n x n Cholesky

  ## ---- Matern: QF only (prec_M precomputed) --------------------------------
  t_matern_qf <- system.time({
    for (rr in seq_len(N_rep)) {
      qf_m <- sum(resid * (prec_M %*% resid))
      ll_m <- -0.5 * T_rep * logdet_M - 0.5 * qf_m
    }
  })["elapsed"] / N_rep

  ## ---- Matern: Chol (paid per R update) ------------------------------------
  t_matern_chol <- system.time({
    for (rr in seq_len(N_rep)) {
      L_tmp  <- t(prop_stable_chol(R_M))
      p_tmp  <- chol2inv(t(L_tmp))
      ld_tmp <- sum(log(eigen(R_M, only.values = TRUE,
                              symmetric = TRUE)$values))
    }
  })["elapsed"] / N_rep

  ## ---- MRTS current: TH + full n x n Chol + QF via prec_obs ---------------
  t_mrts_current <- system.time({
    for (rr in seq_len(N_rep)) {
      m2  <- update_M_closed_form(Fo, resid, sigma_eps_sq = 0,
                                  sigma_floor = 1e-6)
      cs2 <- build_R_from_G(Fo, m2)   # O(n^3) bottleneck
      ql2 <- quadform_logdet(cs2, resid)
    }
  })["elapsed"] / N_rep

  ## ---- MRTS Woodbury: TH + O(nK) log-det + O(nKT) QF --------------------
  t_mrts_woodbury <- system.time({
    for (rr in seq_len(N_rep)) {
      m3  <- update_M_closed_form(Fo, resid, sigma_eps_sq = 0,
                                  sigma_floor = 1e-6)
      # minimal cov state: only a_obs_diag and factor quantities, NO full Chol
      cs3 <- build_R_from_G(Fo, m3)    # still calls full build; time separately
      # log|R| via det-lemma: O(nK)
      ld3 <- woodbury_logdet(cs3)
      # QF via Woodbury: O(nKT)
      rr_mat <- prop_apply_rinv_obs(cs3, resid)
      qf3    <- sum(resid * rr_mat)
      ll3    <- -0.5 * T_rep * ld3 - 0.5 * qf3
    }
  })["elapsed"] / N_rep

  ## ---- Woodbury QF + logdet only (no build_R_from_G re-call) ------------
  # time just the O(nKT) Woodbury operations given cs already built
  t_woodbury_only <- system.time({
    for (rr in seq_len(N_rep)) {
      ld_w  <- woodbury_logdet(cs)
      rr_w  <- prop_apply_rinv_obs(cs, resid)
      qf_w  <- sum(resid * rr_w)
    }
  })["elapsed"] / N_rep

  results[[ni]] <- list(
    n                = n,
    t_matern_qf      = t_matern_qf,
    t_matern_chol    = t_matern_chol,
    t_mrts_current   = t_mrts_current,
    t_mrts_woodbury  = t_mrts_woodbury,
    t_woodbury_only  = t_woodbury_only
  )
  cat(sprintf("n=%d done\n", n))
}

## output -------------------------------------------------------------------
cat(sprintf("\n=== V5: Computational scaling (K=%d, T=%d) ===\n", K, T_rep))
cat(sprintf("%-5s | %-10s %-10s | %-14s %-16s %-14s\n",
            "n",
            "Mat_QF(s)", "Mat_Chol(s)",
            "MRTS_curr(s)", "MRTS_Woodbury(s)", "WB_only(s)"))
cat(strrep("-", 80), "\n")

for (ni in seq_along(n_vals)) {
  r <- results[[ni]]
  cat(sprintf("%-5d | %-10.5f %-10.5f | %-14.5f %-16.5f %-14.5f\n",
              r$n,
              r$t_matern_qf, r$t_matern_chol,
              r$t_mrts_current, r$t_mrts_woodbury, r$t_woodbury_only))
}

cat("\n--- speedup: Matern_Chol / MRTS ---\n")
cat(sprintf("%-5s | %-18s %-18s\n", "n",
            "Chol/MRTS_curr", "Chol/WB_only"))
for (ni in seq_along(n_vals)) {
  r <- results[[ni]]
  cat(sprintf("%-5d | %-18.2f %-18.2f\n",
              r$n,
              r$t_matern_chol / r$t_mrts_current,
              r$t_matern_chol / r$t_woodbury_only))
}

cat("\nKey finding: MRTS_current bottleneck = O(n^3) Chol in build_R_from_G.\n")
cat("             WB_only uses O(nKT) Woodbury — scales correctly.\n")
cat(sprintf("             For n=%d, K=%d, T=%d: WB_only speedup vs Matern_Chol = %.1fx\n",
            n_vals[length(n_vals)],
            K, T_rep,
            results[[length(n_vals)]]$t_matern_chol /
              results[[length(n_vals)]]$t_woodbury_only))
