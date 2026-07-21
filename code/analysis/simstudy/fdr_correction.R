# fdr_correction.R
# ---------------------------------------------------------------------------
# Replace the hand-waved "does not survive multiplicity" with a formal
# Benjamini-Hochberg (BH) FDR correction over the family of paired model-
# comparison p-values reported across the whole UQ analysis (ozone site
# bootstrap + simulation + time-block). Reports how many are nominally
# significant, how many survive BH at FDR 0.05, and the adjusted p of the
# few cells that were ever nominally significant.
#
# Writes (relative to repo code/analysis): ../fdr_correction_summary.csv here,
#   actually -> simstudy/output/results/fdr_correction.csv
# ---------------------------------------------------------------------------

rm(list = ls())
.this <- local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), winslash = "/", mustWork = FALSE)) else "." })
setwd(.this)
root <- normalizePath(file.path(.this, ".."), winslash = "/", mustWork = FALSE)  # code/analysis
rd <- function(...) { f <- file.path(root, ...); if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else NULL }

gather <- function() {
  rows <- list()
  add <- function(arm, score, label, p) if (length(p) && is.finite(p)) rows[[length(rows) + 1]] <<-
    data.frame(arm = arm, score = score, comparison = label, raw_p = p, stringsAsFactors = FALSE)

  OZ <- "ozone/US-all-auto/output/us-all-auto/tables"
  for (split in c("200-200", "300-100")) {
    br <- rd(OZ, sprintf("brier_bootstrap_pairs_%s.csv", split))
    if (!is.null(br)) for (i in seq_len(nrow(br)))
      add(sprintf("ozone-%s", split), "brier", sprintf("%s @%.3f", br$pair[i], br$level[i]), br$boot_p[i])
    for (scn in c("crps", "twcrps")) {
      d <- rd(OZ, sprintf("%s_bootstrap_pairs_%s.csv", scn, split))
      if (!is.null(d)) for (i in seq_len(nrow(d))) add(sprintf("ozone-%s", split), scn, d$pair[i], d$boot_p[i])
    }
  }
  # simulation (Wilcoxon one-sided); Brier file uses w_p
  sb <- rd("simstudy/output/results/ar2_paired_tests.csv")
  if (!is.null(sb)) for (i in seq_len(nrow(sb)))
    add("sim", "brier", sprintf("set%d %s @%.2f", sb$setting[i], sb$contrast[i], sb$level[i]), sb$w_p[i])
  for (scn in c("crps", "twcrps")) {
    d <- rd(sprintf("simstudy/output/results/ar2_%s_paired.csv", scn))
    if (!is.null(d)) for (i in seq_len(nrow(d)))
      add("sim", scn, sprintf("set%d %s", d$setting[i], d$contrast[i]), d$wilcox_p_1sided[i])
  }
  # time-block (Wilcoxon one-sided)
  TB <- "time_block_forecast/block1_positive_control/output/tables"
  for (tag in c("", "_crps", "_twcrps")) {
    d <- rd(TB, sprintf("blk1_paired_ci%s.csv", tag))
    scn <- if (tag == "") "brier" else sub("_", "", tag)
    if (!is.null(d)) for (i in seq_len(nrow(d)))
      add("timeblock", scn, sprintf("set%d %s %s", d$setting[i], d$leads[i], d$contrast[i]), d$wilcox_p_1sided[i])
  }
  do.call(rbind, rows)
}

df <- gather()
df$p_BH <- p.adjust(df$raw_p, method = "BH")
df$p_BY <- p.adjust(df$raw_p, method = "BY")   # Benjamini-Yekutieli, valid under dependence

m <- nrow(df)
cat(sprintf("Family size m = %d paired comparisons\n", m))
cat(sprintf("Nominally significant (raw p < 0.05): %d\n", sum(df$raw_p < 0.05)))
cat(sprintf("Survive BH at FDR 0.05: %d\n", sum(df$p_BH < 0.05)))
cat(sprintf("Survive BY at FDR 0.05: %d\n", sum(df$p_BY < 0.05)))
cat(sprintf("Smallest raw p = %.4f  ->  BH = %.3f,  BY = %.3f\n",
            min(df$raw_p), min(df$p_BH), min(df$p_BY)))

cat("\n-- cells that were nominally significant (raw p < 0.05) --\n")
sig <- df[df$raw_p < 0.05, ]; sig <- sig[order(sig$raw_p), ]
print(within(sig, { raw_p <- round(raw_p, 4); p_BH <- round(p_BH, 3); p_BY <- round(p_BY, 3) }), row.names = FALSE)

dir.create("output/results", recursive = TRUE, showWarnings = FALSE)
write.csv(df[order(df$raw_p), ], "output/results/fdr_correction.csv", row.names = FALSE)
cat("\nWritten: output/results/fdr_correction.csv\n")
