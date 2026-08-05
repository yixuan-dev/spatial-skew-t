#########################################################################
# fit_assertions.R -- the A / A' / B / C self-consistency assertions for
# time-block fits, and the guard this study gates on.
#
#   A  z self-consistency : mean(z) ~ mean(1/sqrt(tau)) * sqrt(2/pi)
#                           (the model IS z ~ HalfNormal(0, 1/sqrt(tau)))
#   A' z temporal-sd      : sd_t(posterior-mean z path) ~ sigma*sqrt(1-2/pi)
#                           (transverse to the reflected ridge: a reflected
#                           fit keeps the LEVEL of z but shrinks its
#                           temporal variation by 1/|lambda|; see
#                           tex/lambda_phiz_ridge/lambda_phiz_ridge.tex S6)
#   B  truth recovery     : lambda/beta0 near truth, lambda sign correct
#                           (simulation only)
#   C  predictive spread  : sd(yhat[,,h]) <= k * marginal_sd for all h
#                           (a stationary AR forecast converges to climatology)
#
# B AND C are the pre-registered guard of run_block1_guarded.R. A and A'
# are recorded but are not gates.
#
# HISTORY. These assertions used to live in simstudy/fit_diag_utils.R and
# were sourced across the directory boundary. They were moved here on
# 2026-08-05 so that the assertions sit with the only study that acts on
# them: the simstudy arm stopped recording the flags when the
# lambda ~ HN(0, 20) prior removed the reflected ridge (commit bcc2d39)
# and never gated on them. The arithmetic is unchanged, and
# check_fit_consistency() keeps its signature and its emitted column
# order, so the blk1_diag*.csv schema is unaffected.
#
# The numeric half of the record still comes from the shared summariser,
# which is descriptive and stays with the simstudy pipeline that writes
# it into every score cache.
#########################################################################

# Callers run from two different working directories -- this directory
# (run_block1_guarded.R, fit_diagnostics.R) and hn_prior_experiment/ --
# so the summariser is located rather than assumed. ar2_load.R does
# rm(list = ls()), so on a re-source the exists() guard is correctly
# false and this reloads.
if (!exists("fit_diag_summary", mode = "function")) {
  .fda_cand <- c("../simstudy/fit_diag_utils.R",
                 "../../simstudy/fit_diag_utils.R")
  .fda_hit <- .fda_cand[file.exists(.fda_cand)]
  if (!length(.fda_hit)) {
    stop("fit_assertions.R: cannot find simstudy/fit_diag_utils.R from ",
         getwd(), call. = FALSE)
  }
  source(.fda_hit[1])
  rm(.fda_cand, .fda_hit)
}

# Returns a one-row data.frame: the numeric summaries of
# fit_diag_summary(), followed by the four logical pass flags.
# A' is REPORTED (Aprime_sdz + sdz_ratio) but not intended as a hard gate.
check_fit_consistency <- function(fit, yhat, truth, marginal_sd,
                                  data_mean, tol_z = 0.40, tol_lam = 0.60,
                                  k_spread = 3) {
  d <- fit_diag_summary(fit, yhat, data_mean = data_mean)

  # A: the z self-consistency ratio, thresholded.
  A <- abs(d$z_ratio - 1) < tol_z

  # A': informative flag, not a gate.
  Aprime <- d$sdz_ratio > 0.5

  # B: truth recovery (skip if truth is NULL, i.e. real data).
  if (!is.null(truth)) {
    B <- (sign(d$lambda) == sign(truth$lambda)) &&
      (abs(d$lambda / truth$lambda - 1) < tol_lam) &&
      (abs(d$beta0 - truth$beta0) < 0.5 * abs(truth$beta0) + 2)
  } else {
    B <- NA
  }

  # C: predictive spread never exceeds k * marginal.
  C <- d$sd_lead_max <= k_spread * marginal_sd

  cbind(d, data.frame(
    A_zconsist = A, Aprime_sdz = Aprime, B_truth = B, C_spread = C,
    stringsAsFactors = FALSE
  ))
}
