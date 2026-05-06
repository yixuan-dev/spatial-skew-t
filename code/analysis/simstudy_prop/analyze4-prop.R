#########################################################################
# analyze4-prop.R
#
# Post-processing for output/results/scores4-prop.RData.
# Mirrors the analysis done in ../simstudy/results.R (mean/median scores,
# relative-to-Gaussian Brier / Quantile, Friedman tests, lambda coverage,
# diagnostic plots), with an extra axis for the prop basis rank.
#
# Score arrays are indexed: [quantile, dataset, method, prop_k].
#
# Outputs:
#   output/results/simresults4-prop.RData
#   output/tables/score_mean4-prop.csv
#   output/tables/score_median4-prop.csv
#   output/tables/score_rel_gauss4-prop.csv
#   output/tables/best_method_per_K4-prop.csv
#   output/tables/lambda_coverage4-prop.csv
#   output/plots/bs_rel_gauss_by_quantile-K<k>.pdf  (one per prop_k)
#   output/plots/qs_rel_gauss_by_quantile-K<k>.pdf
#   output/plots/bs_mean_vs_K-q<qq>.pdf             (one per selected quantile)
#   output/plots/qs_mean_vs_K-q<qq>.pdf
#   output/plots/lambda_ci_vs_dataset-method<m>-K<k>.pdf
#########################################################################

rm(list = ls())

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]),
                               winslash = "/", mustWork = FALSE)
  if (dir.exists(dirname(script_path))) setwd(dirname(script_path))
}

results_dir <- "output/results"
tables_dir  <- "output/tables"
plots_dir   <- "output/plots"
for (d in c(results_dir, tables_dir, plots_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

load(file.path(results_dir, "scores4-prop.RData"))
# Provides: quant.score, brier.score, beta.0/1/2, tau.alpha, tau.beta,
#           rho, nu, gamma, lambda, elapsed_sec, probs, intervals,
#           prop_ks, datasets, methods, setting

n_probs   <- length(probs)
n_sets    <- length(datasets)
n_methods <- length(methods)
n_ks      <- length(prop_ks)

method_label <- c("1: Gaussian",
                  "2: skew-t (K=1)",
                  "3: t (K=1, T=q0.80)",
                  "4: skew-t (K=5)",
                  "5: t (K=5, T=q0.80)")

# ---- mean / median scores -------------------------------------------
mean_apply <- function(a) apply(a, c(1, 3, 4), mean,   na.rm = TRUE)
med_apply  <- function(a) apply(a, c(1, 3, 4), median, na.rm = TRUE)

bs_mean <- mean_apply(brier.score)   # [probs, method, prop_k]
qs_mean <- mean_apply(quant.score)
bs_med  <- med_apply(brier.score)
qs_med  <- med_apply(quant.score)

dimnames(bs_mean) <- dimnames(qs_mean) <-
  dimnames(bs_med) <- dimnames(qs_med) <-
  list(quantile = as.character(probs),
       method   = as.character(methods),
       prop_k   = as.character(prop_ks))

# ---- long table of mean / median scores -----------------------------
score_long <- function(arr_mean, arr_med, score) {
  out <- expand.grid(quantile = probs, method = methods, prop_k = prop_ks,
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  out$score <- score
  out$mean  <- as.vector(arr_mean)
  out$median <- as.vector(arr_med)
  out
}
score_table <- rbind(
  score_long(bs_mean, bs_med, "brier"),
  score_long(qs_mean, qs_med, "quant")
)
write.csv(score_table[, c("score", "method", "prop_k", "quantile",
                          "mean", "median")],
          file.path(tables_dir, "score_mean4-prop.csv"), row.names = FALSE)

# ---- relative scores vs Gaussian (method 1) at the same prop_k ------
rel_to_gauss <- function(arr) {
  # arr[probs, method, prop_k]
  g <- arr[, 1, , drop = FALSE]                    # gauss only
  arr_div <- sweep(arr, c(1, 3), array(g, dim = c(dim(arr)[1], 1, dim(arr)[3])), "/")
  arr_div  # method 1 column will be exactly 1
}
bs_rel_mean <- rel_to_gauss(bs_mean)
qs_rel_mean <- rel_to_gauss(qs_mean)
bs_rel_med  <- rel_to_gauss(bs_med)
qs_rel_med  <- rel_to_gauss(qs_med)

rel_long <- function(arr_mean, arr_med, score) {
  out <- expand.grid(quantile = probs, method = methods, prop_k = prop_ks,
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  out$score   <- score
  out$rel_mean   <- as.vector(arr_mean)
  out$rel_median <- as.vector(arr_med)
  out
}
rel_table <- rbind(
  rel_long(bs_rel_mean, bs_rel_med, "brier"),
  rel_long(qs_rel_mean, qs_rel_med, "quant")
)
write.csv(rel_table[, c("score", "method", "prop_k", "quantile",
                        "rel_mean", "rel_median")],
          file.path(tables_dir, "score_rel_gauss4-prop.csv"),
          row.names = FALSE)

# ---- best method-K combo per quantile (lowest mean) -----------------
best_combo <- function(arr_mean, score) {
  out <- data.frame(score = score, quantile = probs,
                    best_method = NA_integer_,
                    best_prop_k = NA_integer_,
                    best_value  = NA_real_)
  for (i in seq_len(n_probs)) {
    mat <- arr_mean[i, , ]                  # method x prop_k
    idx <- which(mat == min(mat, na.rm = TRUE), arr.ind = TRUE)[1, ]
    out$best_method[i] <- methods[idx[1]]
    out$best_prop_k[i] <- prop_ks[idx[2]]
    out$best_value[i]  <- mat[idx[1], idx[2]]
  }
  out
}
best_table <- rbind(best_combo(bs_mean, "brier"),
                    best_combo(qs_mean, "quant"))
write.csv(best_table, file.path(tables_dir, "best_method_per_K4-prop.csv"),
          row.names = FALSE)

# ---- lambda coverage for skew-t methods (true lambda = 3) -----------
# intervals: c(0.01, 0.025, 0.05, 0.1, 0.9, 0.95, 0.975, 0.99) -> idx 2 & 7 = 95% CI
ci_lo_idx <- which(abs(intervals - 0.025) < 1e-12)
ci_hi_idx <- which(abs(intervals - 0.975) < 1e-12)
lambda_true <- 3
skew_methods <- c(2L, 4L)
cov_rows <- list()
for (m in skew_methods) {
  mi <- match(m, methods)
  for (ki in seq_along(prop_ks)) {
    lo <- lambda[ci_lo_idx, , mi, ki]
    hi <- lambda[ci_hi_idx, , mi, ki]
    cover <- mean(lo <= lambda_true & hi >= lambda_true, na.rm = TRUE)
    width <- mean(hi - lo, na.rm = TRUE)
    cov_rows[[length(cov_rows) + 1L]] <- data.frame(
      method = m, prop_k = prop_ks[ki],
      coverage_95 = cover, mean_width_95 = width,
      n = sum(!is.na(lo) & !is.na(hi))
    )
  }
}
cov_table <- do.call(rbind, cov_rows)
write.csv(cov_table, file.path(tables_dir, "lambda_coverage4-prop.csv"),
          row.names = FALSE)

# ---- save analysis objects ------------------------------------------
save(bs_mean, qs_mean, bs_med, qs_med,
     bs_rel_mean, qs_rel_mean, bs_rel_med, qs_rel_med,
     score_table, rel_table, best_table, cov_table,
     probs, prop_ks, methods, datasets, intervals, setting,
     file = file.path(results_dir, "simresults4-prop.RData"))

# ---- plots ----------------------------------------------------------
mlty <- c(1, 1, 3, 3, 5)
mpch <- c(21, 22, 23, 24, 25)
mcol <- c("gray30", "firebrick4", "dodgerblue4", "firebrick1", "dodgerblue1")
mbg  <- c("gray70", "firebrick2", "dodgerblue2", "firebrick1", "dodgerblue1")

# (1) relative Brier / Quant vs quantile, lines per method, one PDF per prop_k
plot_rel_vs_quantile <- function(arr_rel, score_lab, file_prefix) {
  for (ki in seq_along(prop_ks)) {
    pdf_file <- file.path(plots_dir,
      sprintf("%s_rel_gauss_by_quantile-K%d.pdf", file_prefix, prop_ks[ki]))
    pdf(pdf_file, width = 7, height = 5)
    on.exit(if (dev.cur() != 1) dev.off(), add = TRUE)
    mat <- arr_rel[, , ki]                  # probs x method
    ymin <- min(mat, 1, na.rm = TRUE)
    ymax <- max(mat, 1, na.rm = TRUE)
    plot(probs, mat[, 1], type = "o", ylim = c(ymin, ymax),
         pch = mpch[1], lty = mlty[1], col = mcol[1], bg = mbg[1],
         xlab = "Threshold quantile",
         ylab = sprintf("Relative %s score (vs. Gaussian)", score_lab),
         main = sprintf("Setting %d, prop_k = %d", setting, prop_ks[ki]))
    abline(h = 1, lty = 2, col = "gray60")
    for (j in 2:n_methods) {
      lines(probs, mat[, j], lty = mlty[j], col = mcol[j])
      points(probs, mat[, j], pch = mpch[j], col = mcol[j], bg = mbg[j])
    }
    legend("topleft", legend = method_label, lty = mlty, pch = mpch,
           col = mcol, pt.bg = mbg, cex = 0.7, bty = "n")
    dev.off()
  }
}
plot_rel_vs_quantile(bs_rel_mean, "Brier",    "bs")
plot_rel_vs_quantile(qs_rel_mean, "Quantile", "qs")

# (2) mean Brier / Quant vs prop_k for selected quantiles, lines per method
selected_q_idx <- c(1, 6, 9, 10)             # 0.90, 0.95, 0.98, 0.99
plot_mean_vs_K <- function(arr_mean, score_lab, file_prefix) {
  for (qi in selected_q_idx) {
    pdf_file <- file.path(plots_dir,
      sprintf("%s_mean_vs_K-q%03d.pdf", file_prefix, round(probs[qi] * 1000)))
    pdf(pdf_file, width = 7, height = 5)
    mat <- arr_mean[qi, , ]                  # method x prop_k
    ymin <- min(mat, na.rm = TRUE)
    ymax <- max(mat, na.rm = TRUE)
    plot(prop_ks, mat[1, ], type = "o", ylim = c(ymin, ymax),
         pch = mpch[1], lty = mlty[1], col = mcol[1], bg = mbg[1],
         xlab = "prop basis rank K",
         ylab = sprintf("Mean %s score", score_lab),
         main = sprintf("Setting %d, q = %.3f", setting, probs[qi]))
    for (j in 2:n_methods) {
      lines(prop_ks, mat[j, ], lty = mlty[j], col = mcol[j])
      points(prop_ks, mat[j, ], pch = mpch[j], col = mcol[j], bg = mbg[j])
    }
    legend("topright", legend = method_label, lty = mlty, pch = mpch,
           col = mcol, pt.bg = mbg, cex = 0.7, bty = "n")
    dev.off()
  }
}
plot_mean_vs_K(bs_mean, "Brier",    "bs")
plot_mean_vs_K(qs_mean, "Quantile", "qs")

# (3) lambda 95% CI vs dataset, one PDF per (skew method, prop_k)
for (m in skew_methods) {
  mi <- match(m, methods)
  for (ki in seq_along(prop_ks)) {
    pdf_file <- file.path(plots_dir,
      sprintf("lambda_ci_vs_dataset-method%d-K%d.pdf", m, prop_ks[ki]))
    pdf(pdf_file, width = 8, height = 5)
    lo <- lambda[ci_lo_idx, , mi, ki]
    hi <- lambda[ci_hi_idx, , mi, ki]
    plot(NA, xlim = c(1, n_sets), ylim = range(c(lo, hi, lambda_true), na.rm = TRUE),
         xlab = "dataset", ylab = expression(lambda),
         main = sprintf("95%% CI for lambda — method %d, prop_k = %d", m, prop_ks[ki]))
    abline(h = lambda_true, col = "red", lty = 2)
    segments(seq_len(n_sets), lo, seq_len(n_sets), hi,
             col = ifelse(lo <= lambda_true & hi >= lambda_true, "black", "orange"))
    points(seq_len(n_sets), (lo + hi) / 2, pch = 19, cex = 0.4)
    dev.off()
  }
}

cat("\nWrote tables to", tables_dir, "\n")
cat("Wrote plots  to", plots_dir, "\n")
cat("Saved analysis objects to", file.path(results_dir, "simresults4-prop.RData"), "\n")
