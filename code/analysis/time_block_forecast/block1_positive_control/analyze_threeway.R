#########################################################################
# analyze_threeway.R -- iid / AR(1) / AR(2) three-way comparison on
# settings 5 (near-unit-root) and 7 (LM d=0.45), block 1.
#
# Methods: 1 = iid, 2 = AR(2), 4 = AR(1). Pairs on the datasets that pass
# the guard under ALL THREE methods (AR(1) is a fresh chain, so its flag
# set differs from AR(2)'s). Pre-registered comparisons (leads 1-5 and
# 1-15, one-sided paired t + Wilcoxon):
#   AR(1) - iid    : does one lag help?
#   AR(2) - iid    : does two-lag pooling help?
#   AR(2) - AR(1)  : the VALUE OF THE SECOND LAG (core quantity)
# Prints per-lead three-way tables and the paired tests.
#
#   & $R analyze_threeway.R
#########################################################################

rm(list = ls())
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1]),
    winslash = "/", mustWork = FALSE))
  if (dir.exists(script_dir)) setwd(script_dir)
}
stopifnot(requireNamespace("scoringRules", quietly = TRUE))

settings <- c(5, 7)
datasets <- 1:30              # scan wide; missing files are skipped
n_target <- 10L              # take the first n_target three-way-clean datasets
methods <- c(1, 2, 4)                       # iid, AR(2), AR(1)
mlab <- c(`1` = "iid", `2` = "AR2", `4` = "AR1")
H <- 15L
probs <- c(0.90, 0.92, 0.94, 0.95, 0.96, 0.98, 0.99)
qi95 <- which(probs == 0.95)
lead_short <- 1:5

score_case <- function(yhat, y_val) {
  Hh <- dim(yhat)[3]
  crps_h <- vapply(seq_len(Hh), function(h)
    mean(scoringRules::crps_sample(y = y_val[, h], dat = t(yhat[, , h])),
         na.rm = TRUE), numeric(1))
  thr <- quantile(y_val, probs = probs, na.rm = TRUE)
  b95_h <- vapply(seq_len(Hh), function(h) {
    phat <- colMeans(yhat[, , h] > thr[qi95])
    mean((as.numeric(y_val[, h] > thr[qi95]) - phat)^2, na.rm = TRUE)
  }, numeric(1))
  list(crps = crps_h, b95 = b95_h)
}

one_sided <- function(dvec) {           # H1: mean < 0 (the second method better)
  dvec <- dvec[!is.na(dvec)]
  if (length(dvec) < 3) return(c(n = length(dvec), mean = mean(dvec), t_p = NA, w_p = NA))
  c(n = length(dvec), mean = mean(dvec),
    t_p = t.test(dvec, alternative = "less")$p.value,
    w_p = wilcox.test(dvec, alternative = "less", exact = FALSE)$p.value)
}

# storage: crps/brier per [lead, dataset, method] and a pass flag per (dataset, method)
for (st in settings) {
  crps <- array(NA_real_, c(H, length(datasets), length(methods)))
  brier <- array(NA_real_, c(H, length(datasets), length(methods)))
  pass <- matrix(NA, length(datasets), length(methods),
                 dimnames = list(datasets, mlab[as.character(methods)]))
  for (di in seq_along(datasets)) {
    for (mi in seq_along(methods)) {
      f <- sprintf("results/blk1-%d-%d-%d.RData", st, methods[mi], datasets[di])
      if (!file.exists(f)) { cat("missing:", f, "\n"); next }
      e <- new.env(parent = emptyenv()); load(f, envir = e)
      sc <- score_case(e$res$yhat, e$res$y_val)
      crps[, di, mi] <- sc$crps
      brier[, di, mi] <- sc$b95
      pass[di, mi] <- isTRUE(e$res$chk$pass)
      rm(e); gc()
    }
  }

  clean_all <- which(apply(pass, 1, all)) # datasets clean under all three methods
  clean <- head(clean_all, n_target)      # first n_target by dataset index
  cat(sprintf("\n=================== setting %d ===================\n", st))
  cat("guard pass by (dataset x method):\n")
  print(pass[apply(!is.na(pass), 1, any), , drop = FALSE])
  cat(sprintf("three-way clean available: %s  (n = %d); using first %d: %s\n",
              paste(datasets[clean_all], collapse = ","), length(clean_all),
              length(clean), paste(datasets[clean], collapse = ",")))
  if (length(clean) < 3) { cat("too few clean datasets for tests\n"); next }

  # per-lead three-way means (over clean datasets)
  lc <- apply(crps[, clean, , drop = FALSE], c(1, 3), mean)   # H x method
  lb <- apply(brier[, clean, , drop = FALSE], c(1, 3), mean)
  cat("\n-- CRPS by lead (iid / AR1 / AR2), gaps vs iid and AR2-AR1 --\n")
  cat(sprintf("%-4s %7s %7s %7s | %8s %8s %8s\n", "h", "iid", "AR1", "AR2",
              "AR1-iid", "AR2-iid", "AR2-AR1"))
  for (h in 1:H)
    cat(sprintf("%-4d %7.3f %7.3f %7.3f | %+8.3f %+8.3f %+8.3f\n", h,
                lc[h, 1], lc[h, 3], lc[h, 2],
                lc[h, 3] - lc[h, 1], lc[h, 2] - lc[h, 1], lc[h, 2] - lc[h, 3]))
  cat("\n-- Brier q95 by lead (iid / AR1 / AR2) --\n")
  cat(sprintf("%-4s %7s %7s %7s | %8s %8s %8s\n", "h", "iid", "AR1", "AR2",
              "AR1-iid", "AR2-iid", "AR2-AR1"))
  for (h in 1:H)
    cat(sprintf("%-4d %7.4f %7.4f %7.4f | %+8.4f %+8.4f %+8.4f\n", h,
                lb[h, 1], lb[h, 3], lb[h, 2],
                lb[h, 3] - lb[h, 1], lb[h, 2] - lb[h, 1], lb[h, 2] - lb[h, 3]))

  # paired tests on per-dataset lead-window means
  win_mean <- function(arr, hs, mi) apply(arr[hs, clean, mi, drop = FALSE], 2, mean)
  cat("\n-- paired one-sided tests (H1: 2nd method better), n =", length(clean), "--\n")
  cat(sprintf("%-22s %9s %7s %7s\n", "contrast (leads)", "mean gap", "t p", "W p"))
  do_contrast <- function(arr, hs, a, b, label) {
    d <- win_mean(arr, hs, b) - win_mean(arr, hs, a)   # b - a; H1 b<a
    st <- one_sided(d)
    cat(sprintf("%-22s %+9.4f %7.3f %7.3f\n", label, st["mean"], st["t_p"], st["w_p"]))
  }
  for (sc_name in c("CRPS", "Brier")) {
    arr <- if (sc_name == "CRPS") crps else brier
    cat(sprintf("[%s]\n", sc_name))
    do_contrast(arr, lead_short, 1, 3, "AR1-iid  (1-5)")
    do_contrast(arr, lead_short, 1, 2, "AR2-iid  (1-5)")
    do_contrast(arr, lead_short, 3, 2, "AR2-AR1  (1-5)*")
    do_contrast(arr, 1:H, 1, 3, "AR1-iid  (1-15)")
    do_contrast(arr, 1:H, 1, 2, "AR2-iid  (1-15)")
    do_contrast(arr, 1:H, 3, 2, "AR2-AR1  (1-15)*")
  }
  cat("(* = core quantity: value of the second lag)\n")
}
