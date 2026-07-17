#########################################################################
# guard_effect.R -- quantify what the guard (reseed-on-assertion-failure)
# did to the scores. For every case where the kept fit is NOT attempt 0
# (or was flagged), refit the attempt-0 chain (same seed -> bit-for-bit
# the fit an unguarded pipeline would have used), score it, and compare
# against the kept fit:
#
#   unguarded study = attempt-0 fit for every case
#   guarded study   = the fits run_block1_guarded.R kept
#
# Output: per-case and per-(setting, method) Brier q95 / CRPS comparison.
#########################################################################

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1]),
    winslash = "/", mustWork = FALSE))
  if (dir.exists(script_dir)) setwd(script_dir)
}
source("../../simstudy/ar2_load.R", chdir = TRUE)
source("../simstudy/time_block_helpers.R")
load("../simstudy/simdata.RData", envir = .GlobalEnv)
options(warn = 1)
if (!requireNamespace("scoringRules", quietly = TRUE)) {
  stop("needs scoringRules", call. = FALSE)
}

iters <- 20000
burn <- 10000
update <- 20001
blk <- tbf_blocks(block_seams, block_H, nt)[[1]]
catalog <- get_tbf_method_catalog()
probs <- c(0.90, 0.92, 0.94, 0.95, 0.96, 0.98, 0.99)
qi95 <- which(probs == 0.95)
lead_short <- 1:5

score_case <- function(yhat, y_val) {   # mirror of analyze_blk1.R
  Hh <- dim(yhat)[3]
  crps_h <- vapply(seq_len(Hh), function(h) {
    mean(scoringRules::crps_sample(y = y_val[, h], dat = t(yhat[, , h])),
         na.rm = TRUE)
  }, numeric(1))
  thr <- quantile(y_val, probs = probs, na.rm = TRUE)
  b95_h <- vapply(seq_len(Hh), function(h) {
    phat <- colMeans(yhat[, , h] > thr[qi95])
    mean((as.numeric(y_val[, h] > thr[qi95]) - phat)^2, na.rm = TRUE)
  }, numeric(1))
  c(crps_s = mean(crps_h[lead_short]), crps_a = mean(crps_h),
    b95_s = mean(b95_h[lead_short]), b95_a = mean(b95_h))
}

# cases whose kept fit != attempt 0 (from the run log), incl. flagged ones
cases <- rbind(
  data.frame(setting = 4, method = c(2, 1, 2, 2, 1), dataset = c(4, 5, 5, 6, 10)),
  data.frame(setting = 5, method = c(2, 1, 2, 1, 1, 2, 1, 2),
             dataset = c(3, 5, 6, 7, 8, 8, 2, 2))
)

rows <- list()
for (i in seq_len(nrow(cases))) {
  st <- cases$setting[i]; m <- cases$method[i]; d <- cases$dataset[i]
  spec <- catalog[catalog$method_id == m, , drop = FALSE]
  y.train <- y[, blk$train_times, d, st]
  x.train <- x[, blk$train_times, , drop = FALSE]
  x.block <- x[, blk$test_times, , drop = FALSE]
  y.val <- y[, blk$test_times, d, st]

  set.seed(get_tbf_seed(st, m, d))          # attempt 0, bit-for-bit
  fit <- mcmc(
    y = y.train, x = x.train, s = s,
    method = "t", skew = isTRUE(spec$skew[1]),
    thresh.all = 0, thresh.quant = TRUE, nknots = spec$nknots[1],
    iterplot = FALSE, iters = iters, burn = burn, update = update,
    min.s = c(0, 0), max.s = c(10, 10),
    temporalw = isTRUE(spec$temporal[1]), temporaltau = isTRUE(spec$temporal[1]),
    temporalz = isTRUE(spec$temporal[1]),
    ar2_w = isTRUE(spec$ar2[1]), ar2_tau = isTRUE(spec$ar2[1]),
    ar2_z = isTRUE(spec$ar2[1]),
    rho.upper = 15, nu.upper = 10
  )
  yhat0 <- forecast_block(fit = fit, seam = blk$seam, H = block_H,
                          x_block = x.block, s = s, ar2 = isTRUE(spec$ar2[1]))
  lam0 <- mean(fit$lambda)
  rm(fit); gc()
  s0 <- score_case(yhat0, y.val)
  rm(yhat0); gc()

  env <- new.env(parent = emptyenv())
  load(sprintf("results/blk1-%d-%d-%d.RData", st, m, d), envir = env)
  s1 <- score_case(env$res$yhat, env$res$y_val)
  kept <- env$res$chk
  rm(env); gc()

  rows[[i]] <- data.frame(setting = st, method = m, dataset = d,
                          lam_att0 = lam0, lam_kept = kept$lambda,
                          kept_attempt = kept$attempt, kept_pass = kept$pass,
                          crps_s0 = s0["crps_s"], crps_s1 = s1["crps_s"],
                          b95_s0 = s0["b95_s"], b95_s1 = s1["b95_s"],
                          crps_a0 = s0["crps_a"], crps_a1 = s1["crps_a"],
                          b95_a0 = s0["b95_a"], b95_a1 = s1["b95_a"])
  cat(sprintf(
    "s%d d%d m%d: lam %+6.2f -> %+6.2f | CRPS(1-5) %6.3f -> %6.3f | Brier95(1-5) %7.4f -> %7.4f\n",
    st, d, m, lam0, kept$lambda, s0["crps_s"], s1["crps_s"],
    s0["b95_s"], s1["b95_s"]))
}
eff <- do.call(rbind, rows)
save(eff, file = "output/guard_effect.RData")

cat("\n=== guard effect on the affected cases (0 = unguarded att0, 1 = kept) ===\n")
for (st in c(4, 5)) {
  for (m in 1:2) {
    sub <- eff[eff$setting == st & eff$method == m, ]
    if (!nrow(sub)) next
    cat(sprintf(
      "setting %d method %d (%d cases): CRPS(1-5) %6.3f -> %6.3f | Brier95(1-5) %7.4f -> %7.4f\n",
      st, m, nrow(sub), mean(sub$crps_s0), mean(sub$crps_s1),
      mean(sub$b95_s0), mean(sub$b95_s1)))
  }
}
cat("\nwrote output/guard_effect.RData\n")
