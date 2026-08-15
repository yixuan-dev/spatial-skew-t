# ar2_paired_tests.R
# ---------------------------------------------------------------------------
# Replication-level (method-2) uncertainty for the SPATIAL HOLD-OUT arm of the
# AR(2) simulation experiment (thesis Figure 1, settings 9-12).
#
# The figure shows Brier-vs-quantile curves with no uncertainty band, and the
# text claims "the temporal curves sit essentially on top of their non-temporal
# counterparts at every quantile". This script attaches a paired test to that
# claim: the 10 simulated data sets are the independent replicates, so for each
# matched contrast we take the per-dataset Brier difference and run a paired
# t-test (95% CI) and a Wilcoxon signed-rank test across the 10 datasets.
#
# Reads the per-dataset Brier already cached by scores.R:
#   output/results/scores<setting>.RData  ->  brier.score[prob, dataset, method, 1]
#
# Method ids (helpers.R):  2 = Skew-t K=1 (no TS)   4 = Skew-t K=5 (no TS)
#                          7 = Skew-t K=1 AR(2)     8 = Skew-t K=5 AR(2)
#                          9 = Skew-t K=1 AR(1)    10 = Skew-t K=5 AR(1)
# Contrast A vs B reports Delta = BS_A - BS_B  (negative => A better).
# ---------------------------------------------------------------------------

rm(list = ls())
.this <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), winslash = "/", mustWork = FALSE)) else "."
})
setwd(.this)

settings     <- c(9, 10, 11, 12, 16, 19)
setting_desc <- c(`9` = "K=1 strong", `10` = "K=1 weak", `11` = "K=5 strong", `12` = "K=5 weak",
                  `16` = "K=1 (0.15,0.80)", `19` = "K=1 ARFIMA d=0.45")
target_probs <- c(0.90, 0.95, 0.98)

# matched contrasts: c(A_id, B_id, label)
contrasts_k1 <- list(
  c(7, 2, "AR(2) vs no-TS  (K=1)"),
  c(7, 9, "AR(2) vs AR(1)  (K=1)"),
  c(9, 2, "AR(1) vs no-TS  (K=1)")
)
contrasts_k5 <- list(
  c(8, 4, "AR(2) vs no-TS  (K=5)"),
  c(8, 10, "AR(2) vs AR(1)  (K=5)"),
  c(10, 4, "AR(1) vs no-TS  (K=5)")
)

paired <- function(d) {                       # d = per-dataset Brier difference (A - B)
  d <- d[is.finite(d)]
  n <- length(d)
  if (n < 3) return(list(n = n, mean = if (n) mean(d) else NA, lo = NA, hi = NA, t_p = NA, w_p = NA))
  tt <- t.test(d)                             # two-sided; H0: mean diff = 0
  # exact Wilcoxon at n = 10; this is the test sec:sim-uq describes when it
  # states the smallest attainable two-sided p-value
  wt <- suppressWarnings(wilcox.test(d))
  list(n = n, mean = mean(d), lo = tt$conf.int[1], hi = tt$conf.int[2],
       t_p = tt$p.value, w_p = wt$p.value)
}

rows <- list()
for (st in settings) {
  f <- sprintf("output/results/scores%d.RData", st)
  if (!file.exists(f)) { cat("missing cache:", f, "\n"); next }
  e <- new.env(); load(f, envir = e)
  bs <- e$brier.score                         # prob x dataset x method x mrts_k
  probs <- e$probs
  methods <- e$methods                        # integer ids in dim-3 order
  pos <- function(id) match(id, methods)
  # contrast set follows the GENERATING knot count: 9, 10, 16, 19 are K=1
  # designs, 11 and 12 are K=5. Setting 19 joins the K=1 list now that its
  # nine-model catalog is complete (before 2026-08-10 it had only the
  # Gaussian and the K=1 AR(1)/AR(2) fits, so only 7 vs 9 was available).
  kk <- if (st %in% c(9, 10, 16, 19)) contrasts_k1 else contrasts_k5

  for (pr in target_probs) {
    pi <- which.min(abs(probs - pr))
    for (ct in kk) {
      A <- as.integer(ct[1]); B <- as.integer(ct[2])
      if (is.na(pos(A)) || is.na(pos(B))) next
      d <- bs[pi, , pos(A), 1] - bs[pi, , pos(B), 1]   # length = 10 datasets
      r <- paired(d)
      rows[[length(rows) + 1]] <- data.frame(
        setting = st, desc = setting_desc[as.character(st)], level = probs[pi],
        contrast = ct[3], n = r$n, Delta = r$mean, CI_lo = r$lo, CI_hi = r$hi,
        t_p = r$t_p, w_p = r$w_p, stringsAsFactors = FALSE)
    }
  }
}
out <- do.call(rbind, rows)

# ---- report ----
fmt <- out
for (c0 in c("Delta", "CI_lo", "CI_hi")) fmt[[c0]] <- sprintf("%+.5f", out[[c0]])
fmt$t_p <- sprintf("%.3f", out$t_p); fmt$w_p <- sprintf("%.3f", out$w_p)
fmt$sig <- ifelse(out$CI_lo > 0 | out$CI_hi < 0, "*", "")
cat("=== AR(2) simulation spatial hold-out: paired tests over 10 datasets ===\n")
cat("=== Delta = BS_A - BS_B (negative => A better); 95% paired-t CI; '*' CI excludes 0 ===\n\n")
print(fmt[, c("setting", "desc", "level", "contrast", "n", "Delta", "CI_lo", "CI_hi", "w_p", "sig")],
      row.names = FALSE)

dir.create("output/results", recursive = TRUE, showWarnings = FALSE)
write.csv(out, "output/results/ar2_paired_tests.csv", row.names = FALSE)
cat("\nWritten: output/results/ar2_paired_tests.csv\n")
