#########################################################################
# fit_diag_utils.R -- shared self-consistency assertions for time-block fits.
#
# Extracted from fit_diagnostics.R so that fitting drivers (e.g. the
# block1_positive_control study) can reuse the same checks as guards.
#
#   A  z self-consistency : mean(z) ~ mean(1/sqrt(tau)) * sqrt(2/pi)
#                           (the model IS z ~ HalfNormal(0, 1/sqrt(tau)))
#   A' z temporal-sd      : sd_t(posterior-mean z path) ~ sigma*sqrt(1-2/pi)
#                           (transverse to the reflected ridge: a reflected
#                           fit keeps the LEVEL of z but shrinks its temporal
#                           variation by 1/|lambda|; see
#                           tex/lambda_phiz_ridge/lambda_phiz_ridge.tex S6)
#   B  truth recovery     : lambda/beta0 near truth, lambda sign correct
#                           (simulation only)
#   C  predictive spread  : sd(yhat[,,h]) <= k * marginal_sd for all h
#                           (a stationary AR forecast converges to climatology)
#   D  ridge check        : report beta0 + lambda*mean(z) vs the data mean
#########################################################################

# Returns a one-row data.frame: numeric summaries + logical pass flags.
# A' is REPORTED (Aprime_sdz + sdz_ratio) but not intended as a hard gate yet.
check_fit_consistency <- function(fit, yhat, truth, marginal_sd,
                                  data_mean, tol_z = 0.40, tol_lam = 0.60,
                                  k_spread = 3) {
  a <- mean(fit$tau.alpha) / 2          # internal gamma shape
  b <- mean(fit$tau.beta) / 2
  lam <- mean(fit$lambda)
  beta0 <- mean(fit$beta[, 1])
  zbar <- mean(fit$z)
  # non-temporal fits (e.g. the iid method) carry no phi chains
  get_phi2 <- function(p) {
    if (is.null(p)) return(c(NA_real_, NA_real_))
    cm <- colMeans(as.matrix(p))
    c(cm, NA_real_)[1:2]
  }
  phi.z <- get_phi2(fit$phi.z)
  phi.tau <- get_phi2(fit$phi.tau)

  # A: z self-consistency. E[z] = E[1/sqrt(tau)] * sqrt(2/pi) on the SAME cells.
  z_pred <- mean(1 / sqrt(fit$tau)) * sqrt(2 / pi)
  z_ratio <- zbar / z_pred
  A <- abs(z_ratio - 1) < tol_z

  # A': temporal variation of the fitted z path vs the model's own marginal.
  # fit$z is iters x nt (K = 1); the posterior-mean path is colMeans.
  z_path <- colMeans(fit$z)
  sdz_pred <- mean(1 / sqrt(fit$tau)) * sqrt(1 - 2 / pi)
  sdz_ratio <- sd(z_path) / sdz_pred
  Aprime <- sdz_ratio > 0.5             # informative flag, not a gate

  # B: truth recovery (skip if truth is NULL, i.e. real data).
  if (!is.null(truth)) {
    B <- (sign(lam) == sign(truth$lambda)) &&
      (abs(lam / truth$lambda - 1) < tol_lam) &&
      (abs(beta0 - truth$beta0) < 0.5 * abs(truth$beta0) + 2)
  } else {
    B <- NA
  }

  # C: predictive spread never exceeds k * marginal.
  H <- dim(yhat)[3]
  sd_h <- vapply(seq_len(H), function(h) sd(as.numeric(yhat[, , h])), numeric(1))
  C <- max(sd_h) <= k_spread * marginal_sd

  # D: ridge -- does the fit reproduce the data mean despite bad parts?
  mu_recon <- beta0 + lam * zbar

  data.frame(
    a = a, b = b, lambda = lam, beta0 = beta0, zbar = zbar,
    z_pred = z_pred, z_ratio = z_ratio,
    sdz_ratio = sdz_ratio,
    phi1.z = phi.z[1], phi2.z = phi.z[2],
    phi1.tau = phi.tau[1], phi2.tau = phi.tau[2],
    sd_lead1 = sd_h[1], sd_lead_max = max(sd_h),
    mu_recon = mu_recon, data_mean = data_mean,
    A_zconsist = A, Aprime_sdz = Aprime, B_truth = B, C_spread = C,
    stringsAsFactors = FALSE
  )
}
