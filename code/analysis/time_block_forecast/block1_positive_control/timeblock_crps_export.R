# timeblock_crps_export.R
# ---------------------------------------------------------------------------
# CRPS companion to timeblock_paired_export.R. Recomputes per-dataset CRPS from
# the raw fits (there is no CRPS cache -- cache_blk1_brier.R is Brier-only) for
# the SAME clean-n=10 datasets, then attaches a paired t-test (two-sided 95% CI)
# and a one-sided Wilcoxon (H1: temporal model better) per setting x lead-window
# x contrast. Slow (loads the ~160MB fits); run in the background.
#
# Uses the cache only to pick the identical clean datasets, so the CRPS and
# Brier tables score exactly the same fits.
#
# Writes: output/tables/blk1_paired_ci_crps.csv
# ---------------------------------------------------------------------------

rm(list = ls())
.this <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), winslash = "/", mustWork = FALSE)) else "."
})
setwd(.this)
stopifnot(requireNamespace("scoringRules", quietly = TRUE))

e <- new.env(); load("output/tables/blk1_brier_cache.RData", envir = e)
pass <- e$pass; settings <- e$settings; n_target <- e$n_target
methods <- c(1L, 2L, 4L); mnm <- c("iid", "AR2", "AR1")   # ids + dim order
setting_desc <- c(`5` = "near-unit-root", `7` = "ARFIMA d=0.45")
H <- 15L; windows <- list(`1-5` = 1:5, `1-15` = 1:15)

contrasts <- list(c("AR1","iid","AR(1) - iid"),
                  c("AR2","iid","AR(2) - iid"),
                  c("AR2","AR1","AR(2) - AR(1)"))

crps_by_lead <- function(yhat, y_val) {
  vapply(seq_len(dim(yhat)[3]), function(h)
    mean(scoringRules::crps_sample(y = y_val[, h], dat = t(yhat[, , h])), na.rm = TRUE),
    numeric(1))
}
paired <- function(d) {
  d <- d[is.finite(d)]
  tt <- t.test(d)
  w1 <- suppressWarnings(wilcox.test(d, alternative = "less", exact = FALSE))
  list(n = length(d), mean = mean(d), lo = tt$conf.int[1], hi = tt$conf.int[2],
       t_p2 = tt$p.value, w_p1 = w1$p.value)
}

rows <- list()
for (si in seq_along(settings)) {
  st <- settings[si]
  clean <- head(which(apply(pass[, , si], 1, all)), n_target)
  # per-dataset CRPS-by-lead for each method
  crps <- array(NA_real_, c(H, length(clean), length(methods)))
  for (di in seq_along(clean)) for (mi in seq_along(methods)) {
    f <- sprintf("results/blk1-%d-%d-%d.RData", st, methods[mi], clean[di])
    if (!file.exists(f)) { cat("missing:", f, "\n"); next }
    le <- new.env(); load(f, envir = le)
    crps[, di, mi] <- crps_by_lead(le$res$yhat, le$res$y_val)
    rm(le); gc()
  }
  colnames_m <- mnm
  for (wn in names(windows)) {
    hs <- windows[[wn]]
    wm <- sapply(seq_along(methods), function(mi) apply(crps[hs, , mi, drop = FALSE], 2, mean))
    colnames(wm) <- colnames_m
    for (ct in contrasts) {
      r <- paired(wm[, ct[1]] - wm[, ct[2]])
      rows[[length(rows) + 1]] <- data.frame(
        setting = st, desc = setting_desc[as.character(st)], leads = wn,
        contrast = ct[3], n = r$n, Delta = r$mean, CI_lo = r$lo, CI_hi = r$hi,
        t_p_2sided = r$t_p2, wilcox_p_1sided = r$w_p1, stringsAsFactors = FALSE)
    }
  }
  cat(sprintf("setting %d done (clean: %s)\n", st, paste(clean, collapse = ",")))
}
out <- do.call(rbind, rows)
write.csv(out, "output/tables/blk1_paired_ci_crps.csv", row.names = FALSE)
cat("Written: output/tables/blk1_paired_ci_crps.csv\n")
print(out[, c("setting","leads","contrast","Delta","CI_lo","CI_hi","wilcox_p_1sided")], row.names = FALSE)
