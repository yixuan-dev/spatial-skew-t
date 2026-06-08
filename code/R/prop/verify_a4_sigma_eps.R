# A4 Verification: sigma_eps_sq = 0 is appropriate?
#
# In prop, std_res[:,t] = v_t ~ N(0, R) exactly (A1 T4).  v_t is a pure
# spatial process with NO separate measurement-error layer, so sigma_eps = 0
# is the exact correct choice.
#
# When TH is applied to v_t ~ N(0, R) (correlation matrix, unit diagonal):
#   G_hat  ≈  R   (TH treats R as the fitted covariance)
#   R_hat  =  normalize(G_hat)  ≈  R          [sigma_eps = 0]
#
# If sigma_eps = delta > 0 is misspecified:
#   sigma_xi_hat  ≈  sigma_xi_hat_0  -  delta  (reduced by delta)
#   G_hat diagonal  ≈  1  -  delta             (underestimated by delta/site)
#   R_hat off-diagonal  ≈  R_true[i,j] / (1 - delta)  (inflated by 1/(1-delta))
#
# Tests
#   T1  sigma_eps=0 : R_hat leading eigenvalue ≈ R_true leading eigenvalue
#   T2  sigma_eps=delta : G_hat diagonal shift ≈ -delta  (algebraic check)
#   T3  sigma_eps=delta : R_hat off-diagonal scaling ≈ 1/(1-delta)

suppressMessages(library(autoFRK))
prop_dir <- "d:/Github/spatial-skew-t/code/R/prop"
source(file.path(prop_dir, "prop_utils.R"))
source(file.path(prop_dir, "prop_basis.R"))
source(file.path(prop_dir, "prop_covariance.R"))
source(file.path(prop_dir, "prop_modules.R"))

set.seed(20260531)

## ---- shared setup ----------------------------------------------------------
n  <- 80L
TT <- 300L   # large T for clean algebraic relationships
K  <- 8L

loc   <- cbind(runif(n), runif(n))
basis <- build_basis_matrix(s_obs = loc, k = K)
Fo    <- basis$F_obs
k0    <- basis$rank_kept

# Low nugget so R has meaningful off-diagonal correlations
true_eigs     <- c(2.0, 1.2, 0.5, rep(0, k0 - 3))
true_sigma_xi <- 0.01    # small nugget -> spatial signal dominates

m_upd_true <- list(
  rotation       = diag(k0),
  sigma_xi_sq    = true_sigma_xi,
  lowrank_eigs   = true_eigs,
  effective_rank = 3L,
  projected_eigs = true_eigs,
  trace_s        = NA_real_,
  lstar          = 3L,
  M_hat          = diag(true_eigs),
  F_rot_obs      = Fo,
  F_rot_pred     = NULL
)
cs_true    <- build_R_from_G(F_obs = Fo, M_update = m_upd_true)
R_true     <- cs_true$r_obs          # n x n, unit diagonal
G_true_diag <- cs_true$a_obs_diag   # diag(G_true)

## ---- draw oracle v_t ~ N(0, R) from factor model --------------------------
draw_v <- function(TT) {
  v <- matrix(NA_real_, n, TT)
  for (tt in seq_len(TT)) {
    eta     <- rnorm(k0) * sqrt(true_eigs)
    xi      <- rnorm(n) * sqrt(true_sigma_xi)
    raw     <- as.vector(Fo %*% eta) + xi
    v[, tt] <- raw / sqrt(G_true_diag)
  }
  v
}
v_oracle <- draw_v(TT)

## ---- helper: TH + build R_hat ---------------------------------------------
run_th <- function(v, sigma_eps_sq) {
  m_upd <- update_M_closed_form(
    F_obs         = Fo,
    residuals_std = v,
    sigma_eps_sq  = sigma_eps_sq,
    sigma_floor   = 0            # disable floor to expose pure sigma_eps effect
  )
  cs <- build_R_from_G(Fo, m_upd)
  list(
    g_diag      = cs$a_obs_diag,
    r_obs       = cs$r_obs,
    sigma_xi_sq = m_upd$sigma_xi_sq
  )
}

## ---- T1: sigma_eps = 0 (correct) ------------------------------------------
th0     <- run_th(v_oracle, 0)
eig1_true <- eigen(R_true,      only.values = TRUE, symmetric = TRUE)$values[1]
eig1_hat  <- eigen(th0$r_obs,   only.values = TRUE, symmetric = TRUE)$values[1]
t1_relerr <- abs(eig1_hat - eig1_true) / eig1_true

## ---- T2: G_hat diagonal shift ≈ -delta ------------------------------------
deltas    <- c(0.05, 0.10, 0.20)
th_d      <- lapply(deltas, function(d) run_th(v_oracle, d))
# shift = G_hat_0 diagonal - G_hat_delta diagonal
shifts_obs  <- sapply(th_d, function(x) mean(th0$g_diag - x$g_diag))
shifts_pred <- deltas

## ---- T3: R_hat off-diagonal scaling ≈ 1/(1-delta) ------------------------
uptri <- upper.tri(R_true)
od_true <- R_true[uptri]
od_hat0 <- th0$r_obs[uptri]

# For each delta: ratio of off-diagonal R_hat to R_hat_0 should ≈ 1/(1-delta)
od_ratios  <- sapply(th_d, function(x) mean(x$r_obs[uptri]) / mean(od_hat0))
od_pred    <- 1 / (1 - deltas)   # predicted scaling

## ---- output ----------------------------------------------------------------
cat("=== A4 Verification: sigma_eps_sq = 0 is appropriate? ===\n")
cat(sprintf("n=%d  T=%d  K0=%d  true sigma_xi2=%.3f  (low nugget)\n\n",
            n, TT, k0, true_sigma_xi))

cat("--- T1: sigma_eps=0 — R_hat leading eigenvalue vs R_true ---\n")
cat(sprintf("True  leading eig of R : %.4f\n", eig1_true))
cat(sprintf("Est   leading eig of R : %.4f\n", eig1_hat))
cat(sprintf("Relative error         : %.4f  (expect <0.20)\n\n", t1_relerr))

cat("--- T2: sigma_eps=delta — G_hat diagonal shift vs predicted ---\n")
cat(sprintf("%-8s  %-14s  %-14s  %-10s\n",
            "delta", "shift_obs", "shift_pred(-delta)", "ratio"))
for (ii in seq_along(deltas)) {
  cat(sprintf("%-8.2f  %-14.4f  %-14.4f  %-10.4f\n",
              deltas[ii], shifts_obs[ii], shifts_pred[ii],
              shifts_obs[ii] / shifts_pred[ii]))
}

cat("\n--- T3: R_hat off-diagonal scaling vs 1/(1-delta) ---\n")
cat(sprintf("mean off-diagonal R_true    : %+.6f\n", mean(od_true)))
cat(sprintf("mean off-diagonal R_hat(0)  : %+.6f\n\n", mean(od_hat0)))
cat(sprintf("%-8s  %-14s  %-14s  %-10s\n",
            "delta", "ratio_obs", "1/(1-delta)", "rel_err"))
for (ii in seq_along(deltas)) {
  cat(sprintf("%-8.2f  %-14.4f  %-14.4f  %-10.4f\n",
              deltas[ii], od_ratios[ii], od_pred[ii],
              abs(od_ratios[ii] - od_pred[ii]) / od_pred[ii]))
}

## ---- verdict ---------------------------------------------------------------
cat("\n=== VERDICT ===\n")
t1_ok <- t1_relerr < 0.20
t2_ok <- all(abs(shifts_obs / shifts_pred - 1) < 0.05)
t3_ok <- all(abs(od_ratios - od_pred) / od_pred < 0.10)

cat(sprintf("T1 R_hat leading eig within 20%%      : %s  (%.4f)\n",
            ifelse(t1_ok, "PASS", "FAIL"), t1_relerr))
cat(sprintf("T2 G diag shift ≈ -delta (within 5%%): %s\n",
            ifelse(t2_ok, "PASS", "FAIL")))
cat(sprintf("T3 off-diag scaling ≈ 1/(1-delta)   : %s\n",
            ifelse(t3_ok, "PASS", "FAIL")))
cat(sprintf("OVERALL                              : %s\n",
            ifelse(all(t1_ok, t2_ok, t3_ok), "PASS", "FAIL")))
