# ar2_crps_export.R
# ---------------------------------------------------------------------------
# CRPS companion to ar2_paired_tests.R for the AR(2) simulation spatial
# hold-out arm (settings 9-12). scores.R caches Brier/quant/energy but NOT
# CRPS, so we recompute per-dataset mean CRPS from the fits' predictive draws
# at the held-out test sites, then run a paired t (two-sided 95% CI) + Wilcoxon
# (one-sided, H1: temporal better) across the 10 datasets.
#
# Reuses scores.R's data plumbing (helpers.R, simdata `y`/`obs`, fit.1$yp).
# Writes: output/results/ar2_crps_paired.csv
# ---------------------------------------------------------------------------

rm(list = ls())
.this <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), winslash = "/", mustWork = FALSE)) else "."
})
setwd(.this)
source("./helpers.R")
source("../../R/auxfunctions.R")
stopifnot(requireNamespace("scoringRules", quietly = TRUE))

data_path <- resolve_simstudy_data_path(NULL)   # default simdata.RData
load(data_path)                                  # provides y[site,day,dataset,setting], ntest
obs <- c(rep(TRUE, dim(y)[1] - ntest), rep(FALSE, ntest))   # last ntest sites = test (as in scores.R:231)
results_dir <- derive_results_dir(data_path, "results")

settings <- c(9L, 10L, 11L, 12L)
setting_desc <- c(`9` = "K=1 strong", `10` = "K=1 weak", `11` = "K=5 strong", `12` = "K=5 weak")
datasets <- 1:10
mlab <- c(`2` = "no-TS", `4` = "no-TS", `7` = "AR2", `8` = "AR2", `9` = "AR1", `10` = "AR1")
# matched method sets per setting (K=1 uses 2/7/9; K=5 uses 4/8/10)
meth_of <- function(st) if (st %in% c(9, 10)) c(noTS = 2L, AR2 = 7L, AR1 = 9L) else c(noTS = 4L, AR2 = 8L, AR1 = 10L)
contrasts <- list(c("AR2", "noTS", "AR(2) vs no-TS"),
                  c("AR2", "AR1", "AR(2) vs AR(1)"),
                  c("AR1", "noTS", "AR(1) vs no-TS"))

# per-dataset mean CRPS over held-out test site-days
mean_crps <- function(yp, validate) {
  d <- dim(yp)                                    # iters x np x nt
  dat <- matrix(aperm(yp, c(2, 3, 1)), nrow = d[2] * d[3], ncol = d[1])  # (np*nt) x iters
  yv <- as.vector(validate)                       # site-fastest, matches dat rows
  mean(scoringRules::crps_sample(y = yv, dat = dat), na.rm = TRUE)
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
  crps <- matrix(NA_real_, length(datasets), length(ms), dimnames = list(datasets, names(ms)))
  for (di in seq_along(datasets)) {
    set <- datasets[di]
    validate <- y[!obs, , set, st]
    for (mj in seq_along(ms)) {
      f <- build_simstudy_result_file(results_dir = results_dir, setting_id = st,
                                      method_id = ms[mj], dataset_id = set, mrts_k = NA_integer_)
      if (!file.exists(f)) { cat("missing:", f, "\n"); next }
      env <- new.env(parent = emptyenv()); load(f, envir = env)
      fit <- env$fit.1
      if (!is.null(fit$yp)) crps[di, mj] <- mean_crps(fit$yp, validate)
      rm(env); gc()
    }
  }
  for (ct in contrasts) {
    r <- paired(crps[, ct[1]] - crps[, ct[2]])
    rows[[length(rows) + 1]] <- data.frame(
      setting = st, desc = setting_desc[as.character(st)], contrast = ct[3],
      n = r$n, Delta = r$mean, CI_lo = r$lo, CI_hi = r$hi,
      t_p_2sided = r$t_p2, wilcox_p_1sided = r$w_p1, stringsAsFactors = FALSE)
  }
  cat(sprintf("setting %d done\n", st))
}
out <- do.call(rbind, rows)
write.csv(out, "output/results/ar2_crps_paired.csv", row.names = FALSE)
cat("Written: output/results/ar2_crps_paired.csv\n")
print(out[, c("setting", "contrast", "Delta", "CI_lo", "CI_hi", "wilcox_p_1sided")], row.names = FALSE)
