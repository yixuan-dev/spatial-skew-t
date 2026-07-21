# timeblock_paired_export.R
# ---------------------------------------------------------------------------
# Attach uncertainty to the time-block forecast Brier differences (thesis
# Table 4 / persistence experiment). The 10 clean data sets are the
# independent replicates; for each setting x lead-window x contrast we take the
# per-dataset difference of window-mean q95 Brier and run a paired t-test
# (two-sided 95% CI) plus a two-sided Wilcoxon signed-rank test.
#
# Reads the Brier-only cache (no slow CRPS re-scoring):
#   output/tables/blk1_brier_cache.RData
#     brier[lead 1:15, prob, dataset, method(iid,AR2,AR1), setting(5,7)]
#     pass [dataset, method, setting]   guard flag
# Clean-n selection = first n_target datasets passing under ALL THREE methods
# (identical to cache_blk1_brier.R / blk1_thesis_table.csv).
#
# Writes: output/tables/blk1_paired_ci.csv
# ---------------------------------------------------------------------------

rm(list = ls())
.this <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), winslash = "/", mustWork = FALSE)) else "."
})
setwd(.this)

e <- new.env(); load("output/tables/blk1_brier_cache.RData", envir = e)
brier <- e$brier; pass <- e$pass; probs <- e$probs
settings <- e$settings; n_target <- e$n_target
mlab <- c(iid = "iid", AR2 = "AR2", AR1 = "AR1")   # dim-4 order in the cache
qi <- which.min(abs(probs - 0.95))                 # q95 index
windows <- list(`1-5` = 1:5, `1-15` = 1:15)
setting_desc <- c(`4` = "AR(2) (0.80,-0.35)", `5` = "near-unit-root", `7` = "ARFIMA d=0.45")

# contrast c(A, B): Delta = BS_A - BS_B, negative => A better (temporal wins)
contrasts_full <- list(
  c("AR1", "iid", "AR(1) - iid"),
  c("AR2", "iid", "AR(2) - iid"),
  c("AR2", "AR1", "AR(2) - AR(1)")
)
contrasts_iid_only <- list(c("AR2", "iid", "AR(2) - iid"))  # setting 4: no AR(1) fits

paired <- function(d) {
  d <- d[is.finite(d)]; n <- length(d)
  tt <- t.test(d)                                   # two-sided CI
  w2 <- suppressWarnings(wilcox.test(d, exact = FALSE))  # two-sided
  list(n = n, mean = mean(d), lo = tt$conf.int[1], hi = tt$conf.int[2],
       t_p2 = tt$p.value, w_p2 = w2$p.value)
}

rows <- list()
for (si in seq_along(settings)) {
  st <- settings[si]
  # clean = every method actually fit on the dataset passed (setting 4 lacks AR(1))
  clean <- head(which(apply(pass[, , si], 1,
                            function(r) any(!is.na(r)) && all(r[!is.na(r)]))), n_target)
  contrasts <- if (st == 4) contrasts_iid_only else contrasts_full
  for (wn in names(windows)) {
    hs <- windows[[wn]]
    # per-dataset window-mean Brier at q95 for each method (over clean datasets)
    wm <- sapply(names(mlab), function(m)
      apply(brier[hs, qi, clean, which(names(mlab) == m), si, drop = FALSE], 3, mean))
    # wm: clean-datasets x 3 methods
    for (ct in contrasts) {
      d <- wm[, ct[1]] - wm[, ct[2]]
      r <- paired(d)
      rows[[length(rows) + 1]] <- data.frame(
        setting = st, desc = setting_desc[as.character(st)], leads = wn,
        contrast = ct[3], n = r$n,
        BS_A = mean(wm[, ct[1]]), BS_B = mean(wm[, ct[2]]),
        Delta = r$mean, CI_lo = r$lo, CI_hi = r$hi,
        t_p_2sided = r$t_p2, wilcox_p_2sided = r$w_p2, stringsAsFactors = FALSE)
    }
  }
}
out <- do.call(rbind, rows)

fmt <- out
for (c0 in c("Delta", "CI_lo", "CI_hi")) fmt[[c0]] <- sprintf("%+.5f", out[[c0]])
fmt$t_p_2sided <- sprintf("%.3f", out$t_p_2sided)
fmt$wilcox_p_2sided <- sprintf("%.3f", out$wilcox_p_2sided)
fmt$sig2 <- ifelse(out$wilcox_p_2sided < 0.05, "*", "")   # two-sided Wilcoxon
cat("=== Time-block forecast: paired tests over clean n=10 datasets (q95 Brier) ===\n")
cat("=== Delta = BS_A - BS_B (negative => temporal better); CI = two-sided paired-t ===\n")
cat("=== wilcox_p = TWO-sided; '*' two-sided p<0.05 ===\n\n")
print(fmt[, c("setting", "desc", "leads", "contrast", "n", "Delta",
              "CI_lo", "CI_hi", "t_p_2sided", "wilcox_p_2sided", "sig2")],
      row.names = FALSE)

write.csv(out, "output/tables/blk1_paired_ci.csv", row.names = FALSE)
cat("\nWritten: output/tables/blk1_paired_ci.csv\n")
