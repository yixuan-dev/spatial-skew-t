#########################################################################
# compare_brier_prepost.R -- before/after comparison of predictive scores
# across the z.init fix, setting 4 (phi = (0.80, -0.35)).
#
#   before = output/results/scores4_BUGGY.RData  (scored 2026-05-24, from
#            the z.init = 0.01 runs, now archived in results_BUGGY_backup/)
#   after  = output/results/scores4.RData        (rescored from the fixed
#            reruns in results/, 2026-07-15)
#
# Both sides must score the SAME y_val under the SAME threshold RULE. The
# rule is recorded in the cache as `brier_threshold_basis` since
# 2026-07-31; caches written before that carry no such object and are
# implicitly the per-block rule. stopifnot(identical(probs)) can NOT see a
# rule change -- probs is the same grid under both rules -- so the basis
# is asserted separately below. The fixed reruns reuse the pipeline seed,
# so every (lead, quantile, dataset, block) cell pairs exactly.
# Aggregation is done two ways --
# mean AND median over the 25 (dataset, block) pairs -- because the
# residual ~16% exploding blocks (the phi.z near-unit-root issue) inflate
# the mean on the "after" side.
#
#   & $R compare_brier_prepost.R
#########################################################################

rm(list = ls())

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0L) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]),
    winslash = "/", mustWork = FALSE)
  if (dir.exists(dirname(script_path))) setwd(dirname(script_path))
}

load_scores <- function(path) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  env
}
bef <- load_scores("output/results/scores4_BUGGY.RData")
aft <- load_scores("output/results/scores4.RData")

stopifnot(identical(bef$probs, aft$probs))

# The threshold rule must match on both sides; probs alone cannot see it.
threshold_basis <- function(env) {
  if (exists("brier_threshold_basis", envir = env, inherits = FALSE)) {
    env$brier_threshold_basis
  } else {
    "validation_block"   # the implicit rule of every pre-2026-07-31 cache
  }
}
if (!identical(threshold_basis(bef), threshold_basis(aft))) {
  stop(sprintf(paste("threshold rule differs: before = '%s', after = '%s'",
                     "-- the before/after cells are not comparable"),
               threshold_basis(bef), threshold_basis(aft)), call. = FALSE)
}

# both files may carry the full dataset axis (1..50) with NA where no
# result file existed; restrict to the datasets scored on BOTH sides.
scored <- function(env) {
  ds <- dimnames(env$brier.lead)$dataset
  ds[apply(env$brier.lead, 3, function(x) !all(is.na(x)))]
}
use_sets <- intersect(scored(bef), scored(aft))
stopifnot(length(use_sets) > 0)
cat("paired datasets:", paste(use_sets, collapse = ","), "\n\n")

# When both sides carry the thresholds, require them identical on the
# datasets actually compared -- the assertion the header claims.
if (!is.null(bef$brier.thresholds) && !is.null(aft$brier.thresholds)) {
  stopifnot(isTRUE(all.equal(bef$brier.thresholds[, use_sets, drop = FALSE],
                             aft$brier.thresholds[, use_sets, drop = FALSE],
                             tolerance = 0)))
}

probs <- bef$probs
H <- dim(bef$brier.lead)[1]
nsets <- length(use_sets)
nblk <- dim(bef$brier.lead)[5]
M_IID <- 1L; M_AR2 <- 2L   # get_tbf_method_catalog(): 1 = iid, 2 = AR(2)

# after-side fit diagnostics: which (dataset, block) failed the spread
# assertion (the residual phi.z near-unit-root explosions). Chains match
# because both reruns reuse get_tbf_seed().
fail_after <- matrix(FALSE, max(as.integer(use_sets)), nblk)
if (file.exists("output/diag/fit_diag_4.RData")) {
  denv <- new.env(parent = emptyenv())
  load("output/diag/fit_diag_4.RData", envir = denv)
  d <- denv$diag
  for (r in seq_len(nrow(d))) {
    if (!isTRUE(d$C_spread[r])) fail_after[d$dataset[r], d$block[r]] <- TRUE
  }
}

# ---- block-level aggregates -------------------------------------------
# collapse leads first -> one number per (dataset, block); then summarise
# over the 25 pairs with mean and median.
blk_brier <- function(env, m, qi) {
  apply(env$brier.lead[, qi, use_sets, m, , drop = FALSE], c(3, 5), mean)
}
blk_crps <- function(env, m) {
  apply(env$crps.lead[, use_sets, m, , drop = FALSE], c(2, 4), mean)
}

fmt <- function(x, d = 4) formatC(x, format = "f", digits = d)

cat("=========================================================\n")
cat(" Brier score, setting 4: before (z.init=0.01) vs after\n")
cat(" cell = per-(dataset,block) mean over leads 1..", H, "\n", sep = "")
cat(" summarised over ", nsets * nblk, " (dataset,block) pairs\n", sep = "")
cat("=========================================================\n")
cat(sprintf("%-5s | %21s | %21s | %21s\n", "", "AR(2) mean (median)",
            "iid mean (median)", "AR2-iid paired mean"))
cat(sprintf("%-5s | %10s %10s | %10s %10s | %10s %10s\n",
            "q", "before", "after", "before", "after", "before", "after"))
sum_tab <- NULL
for (qi in seq_along(probs)) {
  b2 <- blk_brier(bef, M_AR2, qi); a2 <- blk_brier(aft, M_AR2, qi)
  b1 <- blk_brier(bef, M_IID, qi); a1 <- blk_brier(aft, M_IID, qi)
  cat(sprintf("%-5s | %10s %10s | %10s %10s | %+10.4f %+10.4f\n",
              paste0("q", probs[qi] * 100),
              sprintf("%s(%s)", fmt(mean(b2), 3), fmt(median(b2), 3)),
              sprintf("%s(%s)", fmt(mean(a2), 3), fmt(median(a2), 3)),
              sprintf("%s(%s)", fmt(mean(b1), 3), fmt(median(b1), 3)),
              sprintf("%s(%s)", fmt(mean(a1), 3), fmt(median(a1), 3)),
              mean(b2 - b1), mean(a2 - a1)))
  sum_tab <- rbind(sum_tab, data.frame(
    prob = probs[qi],
    ar2_before_mean = mean(b2), ar2_after_mean = mean(a2),
    ar2_before_med = median(b2), ar2_after_med = median(a2),
    iid_before_mean = mean(b1), iid_after_mean = mean(a1),
    iid_before_med = median(b1), iid_after_med = median(a1),
    gap_before_mean = mean(b2 - b1), gap_after_mean = mean(a2 - a1),
    gap_before_med = median(b2 - b1), gap_after_med = median(a2 - a1)))
}

cat("\n(gap = AR2 - iid on paired (dataset,block) cells; negative = AR2 better)\n")

# ---- CRPS same format ---------------------------------------------------
cb2 <- blk_crps(bef, M_AR2); ca2 <- blk_crps(aft, M_AR2)
cb1 <- blk_crps(bef, M_IID); ca1 <- blk_crps(aft, M_IID)
cat("\n================== CRPS, same protocol ==================\n")
cat(sprintf("AR(2)  mean: %s -> %s   median: %s -> %s\n",
            fmt(mean(cb2), 3), fmt(mean(ca2), 3), fmt(median(cb2), 3), fmt(median(ca2), 3)))
cat(sprintf("iid    mean: %s -> %s   median: %s -> %s\n",
            fmt(mean(cb1), 3), fmt(mean(ca1), 3), fmt(median(cb1), 3), fmt(median(ca1), 3)))
cat(sprintf("AR2-iid paired mean gap: %+0.4f -> %+0.4f   median gap: %+0.4f -> %+0.4f\n",
            mean(cb2 - cb1), mean(ca2 - ca1), median(cb2 - cb1), median(ca2 - ca1)))

# ---- culprit-block detail (q95 Brier + CRPS per (dataset, block)) -------
qi95 <- which(probs == 0.95)
b95 <- blk_brier(bef, M_AR2, qi95); a95 <- blk_brier(aft, M_AR2, qi95)
i95a <- blk_brier(aft, M_IID, qi95)
cat("\n========== per-(dataset, block) AR(2) detail, Brier q95 ==========\n")
cat(sprintf("%-4s %-4s | %9s %9s %9s | %9s %9s | %s\n", "dset", "blk",
            "Br before", "Br after", "iid after", "CRPS bef", "CRPS aft", "flag"))
for (d in seq_len(nsets)) {
  dsid <- as.integer(use_sets[d])
  for (b in seq_len(nblk)) {
    flag <- if (fail_after[dsid, b]) "<<< after-side FAIL C (residual explosion)" else ""
    cat(sprintf("%-4d %-4d | %9s %9s %9s | %9s %9s | %s\n", dsid, b,
                fmt(b95[d, b]), fmt(a95[d, b]), fmt(i95a[d, b]),
                fmt(cb2[d, b], 3), fmt(ca2[d, b], 3), flag))
  }
}

# ---- money plot: lead curves, Brier q95 + CRPS ---------------------------
lead_brier <- function(env, m, qi) {
  apply(env$brier.lead[, qi, use_sets, m, , drop = FALSE], 1, mean)
}
lead_crps <- function(env, m) {
  apply(env$crps.lead[, use_sets, m, , drop = FALSE], 1, mean)
}

dir.create("output/plots", recursive = TRUE, showWarnings = FALSE)
pdf("output/plots/brier_prepost-set4.pdf", width = 10, height = 4.5)
par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
cols <- c(ar2_bef = "#d95f02", ar2_aft = "#1b9e77",
          iid_bef = "#fdbf6f", iid_aft = "#66c2a5")
plot_four <- function(y_ab, y_aa, y_ib, y_ia, ylab, main) {
  ylim <- range(y_ab, y_aa, y_ib, y_ia)
  plot(1:H, y_ab, type = "l", lwd = 2, col = cols["ar2_bef"], lty = 2,
       xlab = "lead h", ylab = ylab, main = main, ylim = ylim)
  lines(1:H, y_aa, lwd = 2.5, col = cols["ar2_aft"])
  lines(1:H, y_ib, lwd = 2, col = cols["iid_bef"], lty = 2)
  lines(1:H, y_ia, lwd = 2.5, col = cols["iid_aft"])
  legend("topleft", bty = "n", lwd = 2, cex = 0.8,
         col = cols, lty = c(2, 1, 2, 1),
         legend = c("AR(2) before (buggy z.init)", "AR(2) after (fixed)",
                    "iid before", "iid after"))
}
plot_four(lead_brier(bef, M_AR2, qi95), lead_brier(aft, M_AR2, qi95),
          lead_brier(bef, M_IID, qi95), lead_brier(aft, M_IID, qi95),
          "Brier score", "Brier q95, setting 4 (mean over datasets & blocks)")
plot_four(lead_crps(bef, M_AR2), lead_crps(aft, M_AR2),
          lead_crps(bef, M_IID), lead_crps(aft, M_IID),
          "CRPS", "CRPS, setting 4")
dev.off()

write.csv(sum_tab, "output/tables/brier_prepost_summary-set4.csv", row.names = FALSE)
cat("\nwrote output/plots/brier_prepost-set4.pdf\n")
cat("wrote output/tables/brier_prepost_summary-set4.csv\n")
