#########################################################################
# compare3-set4-K20.R
#
# Three-way mean Brier / Quantile-score comparison at setting 4:
#   (1) original             : ../simstudy/scores4_0mrts.RData
#   (2) mrts (K = 20)        : ../simstudy/scores4_20mrts.RData
#   (3) prop  (K = 20)       : output/results/scores4-prop.RData [,,,2]
#
# Methods 1..5; raw mean scores (not relative).
#
# Outputs:
#   output/tables/compare3-set4-K20.csv
#   output/plots/compare3-set4-K20-bs-method<m>.pdf
#   output/plots/compare3-set4-K20-qs-method<m>.pdf
#   output/plots/compare3-set4-K20-bs-panel.pdf
#   output/plots/compare3-set4-K20-qs-panel.pdf
#########################################################################

rm(list = ls())

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]),
                               winslash = "/", mustWork = FALSE)
  if (dir.exists(dirname(script_path))) setwd(dirname(script_path))
}

tables_dir <- "output/tables"
plots_dir  <- "output/plots"
for (d in c(tables_dir, plots_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ---- load three score sources ---------------------------------------
load_scores <- function(path) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  env
}
e_orig <- load_scores("../simstudy/scores4_0mrts.RData")
e_mrts <- load_scores("../simstudy/scores4_20mrts.RData")
e_prop <- load_scores("output/results/scores4-prop.RData")

probs <- e_prop$probs
n_methods <- 5L
methods <- 1:n_methods
prop_k <- 20L
ki_prop <- match(prop_k, e_prop$prop_ks)
stopifnot(!is.na(ki_prop))

# average over datasets (dim 2) -> [probs, method]
mean_over_sets <- function(arr3) apply(arr3, c(1, 3), mean, na.rm = TRUE)

bs_orig <- mean_over_sets(e_orig$brier.score)
bs_mrts <- mean_over_sets(e_mrts$brier.score)
bs_prop <- mean_over_sets(e_prop$brier.score[, , , ki_prop])
qs_orig <- mean_over_sets(e_orig$quant.score)
qs_mrts <- mean_over_sets(e_mrts$quant.score)
qs_prop <- mean_over_sets(e_prop$quant.score[, , , ki_prop])

# ---- long table ------------------------------------------------------
long_rows <- list()
for (m in methods) {
  for (i in seq_along(probs)) {
    long_rows[[length(long_rows) + 1L]] <- data.frame(
      method = m, quantile = probs[i],
      brier_original = bs_orig[i, m],
      brier_mrts20   = bs_mrts[i, m],
      brier_prop20   = bs_prop[i, m],
      quant_original = qs_orig[i, m],
      quant_mrts20   = qs_mrts[i, m],
      quant_prop20   = qs_prop[i, m]
    )
  }
}
out_table <- do.call(rbind, long_rows)
write.csv(out_table, file.path(tables_dir, "compare3-set4-K20.csv"),
          row.names = FALSE)

# ---- plotting helpers -----------------------------------------------
method_label <- c("Gaussian",
                  "skew-t (K=1)",
                  "t (K=1, T=q0.80)",
                  "skew-t (K=5)",
                  "t (K=5, T=q0.80)")

# styles for the three model variants
src_lab <- c("original", "mrts (K=20)", "prop (K=20)")
src_lty <- c(1, 1, 1)
src_pch <- c(21, 22, 24)
src_col <- c("gray30", "firebrick4", "dodgerblue4")
src_bg  <- c("gray70", "firebrick1", "dodgerblue1")

draw_three <- function(orig_vec, mrts_vec, prop_vec, ylab, main_txt) {
  ymin <- min(c(orig_vec, mrts_vec, prop_vec), na.rm = TRUE)
  ymax <- max(c(orig_vec, mrts_vec, prop_vec), na.rm = TRUE)
  plot(probs, orig_vec, type = "o",
       lty = src_lty[1], pch = src_pch[1], col = src_col[1], bg = src_bg[1],
       ylim = c(ymin, ymax),
       xlab = "Threshold quantile", ylab = ylab, main = main_txt)
  lines(probs, mrts_vec, lty = src_lty[2], col = src_col[2])
  points(probs, mrts_vec, pch = src_pch[2], col = src_col[2], bg = src_bg[2])
  lines(probs, prop_vec, lty = src_lty[3], col = src_col[3])
  points(probs, prop_vec, pch = src_pch[3], col = src_col[3], bg = src_bg[3])
  legend("topright", legend = src_lab,
         lty = src_lty, pch = src_pch, col = src_col, pt.bg = src_bg,
         bty = "n", cex = 0.8)
}

# ---- per-method individual PDFs --------------------------------------
for (m in methods) {
  pdf(file.path(plots_dir,
                sprintf("compare3-set4-K20-bs-method%d.pdf", m)),
      width = 6, height = 5)
  draw_three(bs_orig[, m], bs_mrts[, m], bs_prop[, m],
             "Mean Brier score",
             sprintf("Setting 4 — method %d: %s", m, method_label[m]))
  dev.off()

  pdf(file.path(plots_dir,
                sprintf("compare3-set4-K20-qs-method%d.pdf", m)),
      width = 6, height = 5)
  draw_three(qs_orig[, m], qs_mrts[, m], qs_prop[, m],
             "Mean Quantile score",
             sprintf("Setting 4 — method %d: %s", m, method_label[m]))
  dev.off()
}

# ---- combined panel figure (5 methods + legend) ---------------------
panel_layout <- function(orig_mat, mrts_mat, prop_mat, ylab) {
  par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3.0, 1.0))
  for (m in methods) {
    draw_three(orig_mat[, m], mrts_mat[, m], prop_mat[, m],
               ylab,
               sprintf("method %d: %s", m, method_label[m]))
  }
  plot.new()
  legend("center", legend = src_lab,
         lty = src_lty, pch = src_pch, col = src_col, pt.bg = src_bg,
         bty = "n", cex = 1.2)
}

pdf(file.path(plots_dir, "compare3-set4-K20-bs-panel.pdf"),
    width = 12, height = 8)
panel_layout(bs_orig, bs_mrts, bs_prop, "Mean Brier score")
dev.off()

pdf(file.path(plots_dir, "compare3-set4-K20-qs-panel.pdf"),
    width = 12, height = 8)
panel_layout(qs_orig, qs_mrts, qs_prop, "Mean Quantile score")
dev.off()

cat("Wrote table to ", file.path(tables_dir, "compare3-set4-K20.csv"), "\n")
cat("Wrote plots to ", plots_dir, "\n")
