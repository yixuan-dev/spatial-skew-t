#########################################################################
# fit_diag_utils.R -- numeric per-fit diagnostic summaries for time-block
# fits. DESCRIPTIVE ONLY: this file contains no assertions and no pass /
# fail flags.
#
# fit_diag_summary() reduces a live `fit` (plus its forecast draws) to the
# one-row numeric record that has to be taken while the fit is still in
# memory -- run-settings.R discards `fit` immediately and the chunked
# driver deletes the whole file after scoring, so anything not summarised
# here needs a full refit to recover.
#
# WHERE THE ASSERTIONS LIVE. The A / A' / B / C self-consistency checks
# that used to be computed here now live with the study that gates on
# them:
#
#   block1_positive_control/fit_assertions.R   check_fit_consistency()
#   block1_positive_control/fit_diagnostics.R  the refit harness
#
# They were moved out on 2026-08-05. This study stopped recording the
# flags when the lambda ~ HN(0, 20) prior removed the reflected ridge they
# guarded against (commit bcc2d39); it never gated on them, and keeping
# their definition here made a diagnostic-only study look like it ran
# assertions it did not run. block1_positive_control DOES gate on B and C,
# so that is where they belong. fit_assertions.R sources this file and
# builds the flags on top of the summary below -- the arithmetic is
# unchanged and the emitted column order is preserved.
#
# The columns z_ratio and sdz_ratio remain here. They are ratios, not
# verdicts: A and A' are the thresholded versions of them, and the
# threshold is what moved.
#########################################################################

# Returns a one-row data.frame of numeric summaries. No logical columns.
fit_diag_summary <- function(fit, yhat, data_mean) {
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

  # z self-consistency ratio. The model IS z ~ HalfNormal(0, 1/sqrt(tau)),
  # so E[z] = E[1/sqrt(tau)] * sqrt(2/pi) on the SAME cells; z_ratio is the
  # realised over the predicted.
  z_pred <- mean(1 / sqrt(fit$tau)) * sqrt(2 / pi)
  z_ratio <- zbar / z_pred

  # Temporal variation of the fitted z path against the model's own
  # marginal. fit$z is iters x nt (K = 1), so the posterior-mean path is
  # colMeans. This ratio is transverse to the reflected lambda ridge: a
  # reflected fit keeps the LEVEL of z but shrinks its temporal variation
  # by 1/|lambda| (tex/lambda_phiz_ridge S6).
  z_path <- colMeans(fit$z)
  sdz_pred <- mean(1 / sqrt(fit$tau)) * sqrt(1 - 2 / pi)
  sdz_ratio <- sd(z_path) / sdz_pred

  # Predictive spread by lead. A stationary AR forecast converges to
  # climatology, so this should stay within a small multiple of the
  # marginal SD -- but the comparison is the caller's business.
  H <- dim(yhat)[3]
  sd_h <- vapply(seq_len(H), function(h) sd(as.numeric(yhat[, , h])), numeric(1))

  # Ridge reconstruction: does the fit reproduce the data mean?
  mu_recon <- beta0 + lam * zbar

  data.frame(
    a = a, b = b, lambda = lam, beta0 = beta0, zbar = zbar,
    z_pred = z_pred, z_ratio = z_ratio,
    sdz_ratio = sdz_ratio,
    phi1.z = phi.z[1], phi2.z = phi.z[2],
    phi1.tau = phi.tau[1], phi2.tau = phi.tau[2],
    sd_lead1 = sd_h[1], sd_lead_max = max(sd_h),
    mu_recon = mu_recon, data_mean = data_mean,
    stringsAsFactors = FALSE
  )
}
