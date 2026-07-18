#########################################################################
# reflect_score_channels.R -- does the reflected ridge damage Brier/CRPS,
# or only the multivariate energy/variogram scores?
#
# For the strongest reflected cells (|lambda_att0| >= 5), refit the
# attempt-0 chain bit-for-bit (the reflected fit an unguarded pipeline
# keeps) and score BOTH it and the kept healthy fit on ALL FOUR scores:
#   Brier q95, CRPS (univariate, per-lead marginal)
#   energy, variogram (multivariate, over the length-15 time vector)
# The kept fit's yhat is already on disk; only the reflected chain is refit.
#
# If reflection hurt only the joint structure, Brier/CRPS would be roughly
# flat and only energy/variogram would move. The tables here decide it.
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
stopifnot(requireNamespace("scoringRules", quietly = TRUE))

iters <- 20000; burn <- 10000; update <- 20001
blk <- tbf_blocks(block_seams, block_H, nt)[[1]]
catalog <- get_tbf_method_catalog()
probs95 <- 0.95

four_scores <- function(yhat, y_val) {
  Hh <- dim(yhat)[3]; ns <- dim(yhat)[2]
  crps_h <- vapply(seq_len(Hh), function(h)
    mean(scoringRules::crps_sample(y = y_val[, h], dat = t(yhat[, , h])),
         na.rm = TRUE), numeric(1))
  thr <- quantile(y_val, probs = probs95, na.rm = TRUE)
  b95_h <- vapply(seq_len(Hh), function(h) {
    phat <- colMeans(yhat[, , h] > thr)
    mean((as.numeric(y_val[, h] > thr) - phat)^2, na.rm = TRUE)
  }, numeric(1))
  es <- vs <- numeric(ns)
  for (i in seq_len(ns)) {
    dat_i <- t(yhat[, i, ])                       # H x iters
    es[i] <- scoringRules::es_sample(y = y_val[i, ], dat = dat_i)
    vs[i] <- scoringRules::vs_sample(y = y_val[i, ], dat = dat_i, p = 0.5)
  }
  c(crps = mean(crps_h[1:5]), brier = mean(b95_h[1:5]),
    energy = mean(es), vario = mean(vs))
}

# strongest reflected cells (lambda_att0 <= -5 from the guard log)
cells <- data.frame(
  setting = c(5, 4, 4, 4, 4),
  method  = c(2, 2, 2, 1, 1),
  dataset = c(3, 6, 5, 5, 10),
  lam0    = c(-15.18, -8.90, -7.28, -7.23, -4.97))

rows <- list()
for (i in seq_len(nrow(cells))) {
  st <- cells$setting[i]; m <- cells$method[i]; d <- cells$dataset[i]
  spec <- catalog[catalog$method_id == m, , drop = FALSE]
  y.train <- y[, blk$train_times, d, st]
  x.train <- x[, blk$train_times, , drop = FALSE]
  x.block <- x[, blk$test_times, , drop = FALSE]
  y.val <- y[, blk$test_times, d, st]

  set.seed(get_tbf_seed(st, m, d))                # attempt 0 = reflected
  fit <- mcmc(y = y.train, x = x.train, s = s, method = "t",
              skew = isTRUE(spec$skew[1]), thresh.all = 0, thresh.quant = TRUE,
              nknots = spec$nknots[1], iterplot = FALSE,
              iters = iters, burn = burn, update = update,
              min.s = c(0, 0), max.s = c(10, 10),
              temporalw = isTRUE(spec$temporal[1]),
              temporaltau = isTRUE(spec$temporal[1]),
              temporalz = isTRUE(spec$temporal[1]),
              ar2_w = isTRUE(spec$ar2[1]), ar2_tau = isTRUE(spec$ar2[1]),
              ar2_z = isTRUE(spec$ar2[1]), rho.upper = 15, nu.upper = 10)
  yhat_r <- forecast_block(fit = fit, seam = blk$seam, H = block_H,
                           x_block = x.block, s = s, ar2 = isTRUE(spec$ar2[1]))
  lam_r <- mean(fit$lambda); rm(fit); gc()
  sc_r <- four_scores(yhat_r, y.val); rm(yhat_r); gc()

  env <- new.env(parent = emptyenv())
  load(sprintf("results/blk1-%d-%d-%d.RData", st, m, d), envir = env)
  sc_h <- four_scores(env$res$yhat, env$res$y_val)
  lam_h <- env$res$chk$lambda; rm(env); gc()

  rows[[i]] <- data.frame(cell = sprintf("s%d d%d m%d", st, d, m),
                          lam_r = lam_r, lam_h = lam_h,
                          rbind(reflected = sc_r, healthy = sc_h),
                          which = c("reflected", "healthy"))
  cat(sprintf("%s | lam %+6.1f->%+5.1f | Brier %.4f->%.4f  CRPS %.3f->%.3f  Energy %.3f->%.3f  Vario %.4f->%.4f\n",
              rows[[i]]$cell[1], lam_r, lam_h,
              sc_r["brier"], sc_h["brier"], sc_r["crps"], sc_h["crps"],
              sc_r["energy"], sc_h["energy"], sc_r["vario"], sc_h["vario"]))
}
res <- do.call(rbind, rows)
save(res, file = "output/reflect_score_channels.RData")

# relative reflection damage per score (reflected vs healthy), % worse
cat("\n=== relative damage of reflection (reflected/healthy - 1), mean over 5 cells ===\n")
refl <- res[res$which == "reflected", ]
heal <- res[res$which == "healthy", ]
for (sc in c("brier", "crps", "energy", "vario")) {
  rel <- refl[[sc]] / heal[[sc]] - 1
  cat(sprintf("%-8s: %+6.0f%%  (per cell: %s)\n", sc, 100 * mean(rel),
              paste(sprintf("%+.0f%%", 100 * rel), collapse = " ")))
}
cat("\nwrote output/reflect_score_channels.RData\n")
