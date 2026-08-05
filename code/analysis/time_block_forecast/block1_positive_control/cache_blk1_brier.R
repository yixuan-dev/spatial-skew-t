#########################################################################
# cache_blk1_brier.R -- one-time Brier-only extraction for the thesis.
#
# Scans results/blk1-<setting>-<method>-<dataset>.RData for settings 5
# (near-unit-root) and 7 (LM d=0.45), methods 1 = iid, 2 = AR(2),
# 4 = AR(1), and stores per-lead Brier scores at the full prob grid in a
# small cache so the thesis figure/table scripts never re-read the ~160MB
# fit files. Mirrors analyze_threeway.R's loading loop and
# analyze_blk1.R's Brier scoring (thresholds from the validation block,
# identical across methods within a dataset); no CRPS, so scoringRules
# is not needed.
#
# Output: output/tables/blk1_brier_cache.RData with
#   brier [lead 1:15, prob, dataset, method, setting]
#   pass  [dataset, method, setting]   guard flag per fit
#   settings, datasets, methods, probs, mlab
# Also prints the clean-n=10 per-lead q95 three-way table per setting, for
# diffing against tex/timeblock_expABC_legacy/timeblock_expABC_legacy.tex
# (renamed 2026-08-04 from timeblock_score_report; that report is now a
# frozen archive of this retired campaign).
#
#   & $R cache_blk1_brier.R
#########################################################################

rm(list = ls())
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1]),
    winslash = "/", mustWork = FALSE))
  if (dir.exists(script_dir)) setwd(script_dir)
}

settings <- c(4L, 5L, 7L)     # 4 = strong short-memory control (iid/AR2 only)
datasets <- 1:30              # scan wide; missing files are skipped
n_target <- 10L
methods <- c(1L, 2L, 4L)      # iid, AR(2), AR(1)
mlab <- c(`1` = "iid", `2` = "AR2", `4` = "AR1")
H <- 15L
probs <- c(0.90, 0.92, 0.94, 0.95, 0.96, 0.98, 0.99)
qi95 <- which(probs == 0.95)

brier <- array(NA_real_,
  dim = c(H, length(probs), length(datasets), length(methods), length(settings)),
  dimnames = list(NULL, sprintf("q%.2f", probs), datasets,
                  mlab[as.character(methods)], settings))
pass <- array(NA,
  dim = c(length(datasets), length(methods), length(settings)),
  dimnames = list(datasets, mlab[as.character(methods)], settings))

for (si in seq_along(settings)) {
  st <- settings[si]
  for (di in seq_along(datasets)) {
    for (mi in seq_along(methods)) {
      f <- sprintf("results/blk1-%d-%d-%d.RData", st, methods[mi], datasets[di])
      if (!file.exists(f)) next
      e <- new.env(parent = emptyenv()); load(f, envir = e)
      yhat <- e$res$yhat            # iters x sites x leads
      y_val <- e$res$y_val          # sites x leads
      thr <- quantile(y_val, probs = probs, na.rm = TRUE)
      for (h in seq_len(H)) {
        for (qi in seq_along(probs)) {
          phat <- colMeans(yhat[, , h] > thr[qi])
          brier[h, qi, di, mi, si] <-
            mean((as.numeric(y_val[, h] > thr[qi]) - phat)^2, na.rm = TRUE)
        }
      }
      pass[di, mi, si] <- isTRUE(e$res$chk$pass)
      rm(e, yhat, y_val); gc()
      cat(sprintf("scored: setting %d method %s dataset %d\n",
                  st, mlab[as.character(methods[mi])], datasets[di]))
    }
  }
}

out_dir <- "output/tables"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
save(brier, pass, settings, datasets, methods, probs, mlab, n_target,
     file = file.path(out_dir, "blk1_brier_cache.RData"))
cat(sprintf("\nWrote %s\n",
            normalizePath(file.path(out_dir, "blk1_brier_cache.RData"),
                          winslash = "/", mustWork = FALSE)))

# ---- verification: clean-n=10 per-lead q95 three-way means -------------
# Column order printed as iid / AR1 / AR2 to match analyze_threeway.R and
# the tables in timeblock_expABC_legacy.tex (tab:expB5brier, tab:expC7brier).
for (si in seq_along(settings)) {
  st <- settings[si]
  # clean = every method actually fit on the dataset passed the guard
  # (setting 4 has no AR(1) fits, so its NA column is ignored)
  clean_all <- which(apply(pass[, , si], 1,
                           function(r) any(!is.na(r)) && all(r[!is.na(r)])))
  clean <- head(clean_all, n_target)
  cat(sprintf("\n=================== setting %d ===================\n", st))
  cat(sprintf("three-way clean available: %s (n = %d); using first %d: %s\n",
              paste(datasets[clean_all], collapse = ","), length(clean_all),
              length(clean), paste(datasets[clean], collapse = ",")))
  lb <- apply(brier[, qi95, clean, , si, drop = FALSE], c(1, 4), mean)
  cat("\n-- Brier q95 by lead (iid / AR1 / AR2) --\n")
  cat(sprintf("%-4s %7s %7s %7s | %8s %8s %8s\n", "h", "iid", "AR1", "AR2",
              "AR1-iid", "AR2-iid", "AR2-AR1"))
  for (h in seq_len(H))
    cat(sprintf("%-4d %7.4f %7.4f %7.4f | %+8.4f %+8.4f %+8.4f\n", h,
                lb[h, 1], lb[h, 3], lb[h, 2],
                lb[h, 3] - lb[h, 1], lb[h, 2] - lb[h, 1], lb[h, 2] - lb[h, 3]))
  for (hs in list(1:5, 1:15)) {
    wm <- apply(brier[hs, qi95, clean, , si, drop = FALSE], 4, mean)
    cat(sprintf("leads %d-%d mean: iid %.4f  AR1 %.4f  AR2 %.4f  (AR2-iid %+.4f, AR2-AR1 %+.4f)\n",
                min(hs), max(hs), wm[1], wm[3], wm[2],
                wm[2] - wm[1], wm[2] - wm[3]))
  }
}
