# ar2_twcrps_export.R
# ---------------------------------------------------------------------------
# Tail-band twCRPS for the AR(2) simulation spatial hold-out arm (settings
# 9-12). Confirms whether the interpolation null (Brier + plain CRPS both ~ns
# for AR(2) vs AR(1)) also holds for the tail-focused twCRPS. Mirrors
# ar2_crps_export.R but the per-cell score is the tail-band Brier integral
#   twCRPS = integral_[q.90,q.99] (F(z) - 1{y<=z})^2 dz
# on a fine z-grid (trapezoid), z from quantile(y[,,set,setting]) as in scores.R.
#
# Writes: output/results/ar2_twcrps_paired.csv
# ---------------------------------------------------------------------------

rm(list = ls())
.this <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), winslash = "/", mustWork = FALSE)) else "."
})
setwd(.this)
source("./helpers.R"); source("../../R/auxfunctions.R")
stopifnot(requireNamespace("scoringRules", quietly = TRUE))

data_path <- resolve_simstudy_data_path(NULL); load(data_path)   # y, ntest
obs <- c(rep(TRUE, dim(y)[1] - ntest), rep(FALSE, ntest))
results_dir <- derive_results_dir(data_path, "results")

settings <- c(9L, 10L, 11L, 12L)
setting_desc <- c(`9` = "K=1 strong", `10` = "K=1 weak", `11` = "K=5 strong", `12` = "K=5 weak")
datasets <- 1:10
band_probs <- c(0.90, 0.99); G <- 31L
meth_of <- function(st) if (st %in% c(9, 10)) c(noTS = 2L, AR2 = 7L, AR1 = 9L) else c(noTS = 4L, AR2 = 8L, AR1 = 10L)
contrasts <- list(c("AR2", "noTS", "AR(2) vs no-TS"),
                  c("AR2", "AR1", "AR(2) vs AR(1)"),
                  c("AR1", "noTS", "AR(1) vs no-TS"))

mean_twcrps <- function(yp, validate, zgrid) {
  d <- dim(yp); ypm <- matrix(yp, nrow = d[1]); yv <- as.vector(validate)
  acc <- numeric(length(yv)); prevBS <- NULL
  for (k in seq_along(zgrid)) {
    BS <- (colMeans(ypm > zgrid[k]) - as.numeric(yv > zgrid[k]))^2
    if (k > 1) acc <- acc + 0.5 * (BS + prevBS) * (zgrid[k] - zgrid[k - 1])
    prevBS <- BS
  }
  mean(acc, na.rm = TRUE)
}
paired <- function(v) {
  v <- v[is.finite(v)]
  tt <- t.test(v); w1 <- suppressWarnings(wilcox.test(v, alternative = "less", exact = FALSE))
  list(n = length(v), mean = mean(v), lo = tt$conf.int[1], hi = tt$conf.int[2],
       t_p2 = tt$p.value, w_p1 = w1$p.value)
}

rows <- list()
for (st in settings) {
  ms <- meth_of(st)
  tw <- matrix(NA_real_, length(datasets), length(ms), dimnames = list(datasets, names(ms)))
  for (di in seq_along(datasets)) {
    set <- datasets[di]
    validate <- y[!obs, , set, st]
    band <- quantile(y[, , set, st], probs = band_probs, na.rm = TRUE)
    zgrid <- seq(band[1], band[2], length.out = G)
    for (mj in seq_along(ms)) {
      f <- build_simstudy_result_file(results_dir = results_dir, setting_id = st,
                                      method_id = ms[mj], dataset_id = set, mrts_k = NA_integer_)
      if (!file.exists(f)) { cat("missing:", f, "\n"); next }
      env <- new.env(parent = emptyenv()); load(f, envir = env)
      if (!is.null(env$fit.1$yp)) tw[di, mj] <- mean_twcrps(env$fit.1$yp, validate, zgrid)
      rm(env); gc()
    }
  }
  for (ct in contrasts) {
    r <- paired(tw[, ct[1]] - tw[, ct[2]])
    rows[[length(rows) + 1]] <- data.frame(
      setting = st, desc = setting_desc[as.character(st)], contrast = ct[3],
      n = r$n, Delta = r$mean, CI_lo = r$lo, CI_hi = r$hi,
      t_p_2sided = r$t_p2, wilcox_p_1sided = r$w_p1, stringsAsFactors = FALSE)
  }
  cat(sprintf("setting %d done\n", st))
}
out <- do.call(rbind, rows)
write.csv(out, "output/results/ar2_twcrps_paired.csv", row.names = FALSE)
cat("Written: output/results/ar2_twcrps_paired.csv\n")
print(out[out$contrast == "AR(2) vs AR(1)", c("setting","contrast","Delta","CI_lo","CI_hi","wilcox_p_1sided")], row.names = FALSE)
