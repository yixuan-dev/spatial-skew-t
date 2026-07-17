#########################################################################
# refit_diag.R -- refit ONE time-block case, SAVE the fit, and diagnose the
# AR(2) forecast divergence ON THE LATENT SCALE.
#
# Why: the original run saved only `forecast`, so fit$tau.alpha / lambda /
# the latent trajectory were never inspectable. Three hypotheses have already
# been REFUTED on the Morris fits (gamma-shape collapse to the grid floor;
# the 1e-6 tau clamp firing; non-stationary phi draws) -- but those are a
# different study. This refits the exact time-block case that blew up
# (setting 4, dataset 1, block 3, seam 110: predictive SD reached 104 vs a
# marginal SD of 4.38) so the mechanism can be observed rather than inferred.
#
# KEY METHODOLOGICAL POINT: the invariant "predictive variance <= marginal
# variance" holds on the LATENT GAUSSIAN scale X only. A convex monotone map
# (here tau = qgamma(pnorm(X))) can legitimately inflate a bounded latent
# spread into a huge one on the y scale. So a big y-scale SD does NOT by
# itself prove the recursion is broken. All checks below are therefore on the
# latent scale.
#
#   & $R refit_diag.R
#########################################################################

# ar2_load.R does rm(list = ls()) and loads the MORRIS simdata; source it first
# (chdir so its relative paths resolve), then overwrite with THIS study's data.
source("../../simstudy/ar2_load.R", chdir = TRUE)
source("./time_block_helpers.R")
load("./simdata.RData", envir = .GlobalEnv)   # y, s, x, nt, block_H, block_seams
options(warn = 1)

setting <- 4L
dataset_id <- 1L
block_id <- 3L      # the exploded block (seam 110)
iters <- 20000
burn <- 10000
update <- 2000

blocks <- tbf_blocks(block_seams, block_H, nt)
blk <- blocks[[block_id]]
cat(sprintf("setting %d, dataset %d, block %d: seam = %d, train = 1:%d, H = %d\n",
            setting, dataset_id, block_id, blk$seam, max(blk$train_times), block_H))

y.d <- y[, , dataset_id, setting]
y.train <- y.d[, blk$train_times, drop = FALSE]
x.train <- x[, blk$train_times, , drop = FALSE]
x.block <- x[, blk$test_times, , drop = FALSE]

set.seed(get_tbf_seed(setting, 2L, dataset_id))   # method 2 = AR(2)
fit <- mcmc(
  y = y.train, x = x.train, s = s,
  method = "t", skew = TRUE,
  thresh.all = 0, thresh.quant = TRUE,
  nknots = 1,
  iterplot = FALSE, iters = iters, burn = burn, update = update,
  min.s = c(0, 0), max.s = c(10, 10),
  temporalw = TRUE, temporaltau = TRUE, temporalz = TRUE,
  ar2_w = TRUE, ar2_tau = TRUE, ar2_z = TRUE,
  z.init = 0.01,
  rho.upper = 15, nu.upper = 10
)
save(fit, blk, file = "refit_diag_fit.RData")
cat("\nfit saved -> refit_diag_fit.RData\n")

## ================= DIAGNOSTICS (latent scale) =========================
a <- fit$tau.alpha / 2      # internal shape actually used by gamma.cop
b <- fit$tau.beta / 2
lam <- fit$lambda
tau <- fit$tau              # [iters, T_o]  observable scale
z <- fit$z

cat("\n===== 1. Does the gamma shape collapse? (true internal a = 3) =====\n")
cat(sprintf("  a : mean %.3f  min %.3f  max %.3f  | %% at grid floor 0.05: %.2f%%  | %% at ceiling 5.0: %.2f%%\n",
            mean(a), min(a), max(a), 100 * mean(a <= 0.051), 100 * mean(a >= 4.999)))
cat(sprintf("  b : mean %.3f  min %.3f  max %.3f\n", mean(b), min(b), max(b)))
cat(sprintf("  lambda : mean %.3f  min %.3f  max %.3f  | %% negative: %.2f%%\n",
            mean(lam), min(lam), max(lam), 100 * mean(lam < 0)))

cat("\n===== 2. Does tau hit the 1e-6 clamp? =====\n")
cat(sprintf("  tau: min %.3e  q0.1%% %.3e  median %.3f  | AT CLAMP: %d / %d (%.4f%%)\n",
            min(tau), quantile(tau, .001), median(tau),
            sum(tau <= 1e-6 + 1e-12), length(tau), 100 * mean(tau <= 1e-6 + 1e-12)))
cat(sprintf("  z  : min %.3e  | AT CLAMP: %d / %d (%.4f%%)\n",
            min(z), sum(z <= 1e-6 + 1e-12), length(z), 100 * mean(z <= 1e-6 + 1e-12)))

cat("\n===== 3. The SEAM STATE on the Gaussian copula scale =====\n")
cat("  (this is what ar2_recurse starts from; it must look like N(0,1))\n")
seam_cols <- c(blk$seam - 1L, blk$seam)
ts_all <- zs_all <- numeric(0)
for (m in seq_len(nrow(tau))) {
  a.m <- a[m]; b.m <- b[m]
  tau.obs <- tau[m, seam_cols]
  z.obs <- z[m, seam_cols]
  ts_all <- c(ts_all, gamma.cop(tau.obs, a.m, b.m))
  zs_all <- c(zs_all, hn.cop(z.obs, sig = 1 / sqrt(tau.obs)))
}
rep_state <- function(nm, v) {
  cat(sprintf("  %-10s mean %+7.3f  sd %6.3f  min %+8.2f  max %+7.2f  | %%|.|>4 : %.3f%%  | non-finite: %d\n",
              nm, mean(v[is.finite(v)]), sd(v[is.finite(v)]),
              min(v[is.finite(v)]), max(v[is.finite(v)]),
              100 * mean(abs(v[is.finite(v)]) > 4), sum(!is.finite(v))))
}
rep_state("tau.seam", ts_all)
rep_state("z.seam", zs_all)

cat("\n===== 4. Is ar2_recurse's OUTPUT marginally N(0,1)? =====\n")
cat("  (the true invariant: stationary AR(2) => marginal variance 1 at every lead)\n")
K <- 1L
TS <- ZS <- matrix(NA_real_, nrow(tau), block_H)
for (m in seq_len(nrow(tau))) {
  a.m <- a[m]; b.m <- b[m]
  tau.obs <- matrix(tau[m, seam_cols], K, 2L)
  z.obs <- matrix(z[m, seam_cols], K, 2L)
  tsm <- gamma.cop(tau.obs, a.m, b.m)
  zsm <- hn.cop(z.obs, sig = 1 / sqrt(tau.obs))
  TS[m, ] <- ar2_recurse(tsm, fit$phi.tau[m, ], block_H, K)
  ZS[m, ] <- ar2_recurse(zsm, fit$phi.z[m, ], block_H, K)
}
cat(sprintf("  %4s | %-28s | %-28s\n", "lead", "tau.star", "z.star"))
cat(sprintf("  %4s | %8s %8s %9s | %8s %8s %9s\n", "", "mean", "sd", "%|.|>4", "mean", "sd", "%|.|>4"))
for (h in c(1, 2, 3, 5, 8, 12, 15)) {
  t1 <- TS[, h]; z1 <- ZS[, h]
  cat(sprintf("  %4d | %8.3f %8.3f %8.2f%% | %8.3f %8.3f %8.2f%%\n", h,
              mean(t1[is.finite(t1)]), sd(t1[is.finite(t1)]), 100 * mean(abs(t1) > 4, na.rm = TRUE),
              mean(z1[is.finite(z1)]), sd(z1[is.finite(z1)]), 100 * mean(abs(z1) > 4, na.rm = TRUE)))
}
cat("\n  VERDICT: if sd -> 1 and mean -> 0, the recursion is FINE and the blow-up\n")
cat("           is created by the convex copula map, not by the AR recursion.\n")

cat("\n===== 5. What sigma = 1/sqrt(tau.fc) does that imply? =====\n")
tau.fc <- qgamma(pnorm(TS), rep(a, block_H), rep(b, block_H))
sig <- 1 / sqrt(tau.fc)
cat(sprintf("  data marginal sd of y = %.2f\n", sd(y[, , , setting])))
for (h in c(1, 5, 15)) {
  s1 <- sig[, h]
  cat(sprintf("  lead %2d : sigma  median %8.3f  q99 %10.2f  max %12.2e  | %% > 20 : %.2f%%\n",
              h, median(s1, na.rm = TRUE), quantile(s1, .99, na.rm = TRUE),
              max(s1, na.rm = TRUE), 100 * mean(s1 > 20, na.rm = TRUE)))
}
