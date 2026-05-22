#########################################################################
# plots-prop.R - Stage 3 of the simstudy_prop post-fit pipeline.
#
# Reads output/results/simresults<setting>-prop<suffix>.RData (produced
# by tables-prop.R) and renders the figures under output/plots/.
#
# Split out of tables-prop.R so the table and figure stages re-run
# independently; mirrors code/analysis/simstudy/plots.R.
#
# Usage:
#   Rscript plots-prop.R --setting=<id>
#                        [--data=<path>]
#
# Examples:
#   Rscript plots-prop.R --setting=4
#   Rscript plots-prop.R --setting=1 --data=simdata_def.RData
#
# Outputs (suffix = "" for simdata.RData, "_def" for simdata_def.RData, ...):
#   output/plots/{bs,qs}_rel_gauss_by_quantile-set<setting><suffix>-K<k>.pdf
#   output/plots/{bs,qs}_mean_vs_K-set<setting><suffix>-q<qq>.pdf
#   output/plots/lambda_ci_vs_dataset-set<setting><suffix>-method<m>-K<k>.pdf
#########################################################################

rm(list = ls())

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0L) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]),
                               winslash = "/", mustWork = FALSE)
  if (dir.exists(dirname(script_path))) setwd(dirname(script_path))
}

source("./prop_simstudy_helpers.R")

# ---- CLI parsing -----------------------------------------------------
cli_args <- commandArgs(trailingOnly = TRUE)
parsed <- extract_leading_flags(cli_args, c("data", "setting"))
flags <- parsed$values

if (is.null(flags$setting) || !nzchar(flags$setting)) {
  stop("plots-prop.R: --setting=<id> is required.", call. = FALSE)
}
setting_id <- as.integer(parse_index_expr(flags$setting, "setting"))
if (length(setting_id) != 1L || setting_id < 1L) {
  stop("plots-prop.R: --setting must be a single positive integer.", call. = FALSE)
}

data_suffix <- if (!is.null(flags$data) && nzchar(flags$data)) {
  derive_data_suffix(flags$data)
} else {
  ""
}

results_dir <- "output/results"
plots_dir   <- "output/plots"
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
}

# ---- load the Stage-2 artifact ---------------------------------------
simresults_file <- file.path(
  results_dir,
  sprintf("simresults%d-prop%s.RData", setting_id, data_suffix)
)
if (!file.exists(simresults_file)) {
  stop(sprintf(
    "Aggregated results not found: %s\n  Run tables-prop.R --setting=%d%s first.",
    simresults_file, setting_id,
    if (nzchar(data_suffix)) sprintf(" --data=%s", flags$data) else ""
  ), call. = FALSE)
}
load(simresults_file)
# Provides: bs_mean/qs_mean/bs_med/qs_med, bs_rel_mean/qs_rel_mean/...,
#           score_long_table, score_table, rel_table, best_table,
#           cov_table, lambda, probs, prop_ks, methods, datasets,
#           intervals, setting, data_suffix

cat(sprintf(
  "plots-prop: setting=%d cache=%s suffix='%s'\n",
  setting_id, simresults_file,
  if (nzchar(data_suffix)) data_suffix else "<none>"
))

n_probs   <- length(probs)
n_sets    <- length(datasets)
n_methods <- length(methods)
n_ks      <- length(prop_ks)

method_catalog <- get_prop_method_catalog()
method_label_for <- function(m) {
  row <- method_catalog[method_catalog$method_id == m, , drop = FALSE]
  if (nrow(row) == 1L) sprintf("%d: %s", m, row$label[1]) else as.character(m)
}
method_label <- vapply(methods, method_label_for, character(1))

# lambda CI bookkeeping (mirrors tables-prop.R's lambda_coverage block)
ci_lo_idx <- which(abs(intervals - 0.025) < 1e-12)
ci_hi_idx <- which(abs(intervals - 0.975) < 1e-12)
lambda_true  <- 3
skew_methods <- intersect(c(2L, 4L), methods)

# ---- figures ---------------------------------------------------------
set_tag <- sprintf("set%d%s", setting_id, data_suffix)
mlty <- c(1, 1, 3, 3, 5)
mpch <- c(21, 22, 23, 24, 25)
mcol <- c("gray30", "firebrick4", "dodgerblue4", "firebrick1", "dodgerblue1")
mbg  <- c("gray70", "firebrick2", "dodgerblue2", "firebrick1", "dodgerblue1")

# (1) relative score vs quantile, lines per method, one PDF per prop_k
plot_rel_vs_quantile <- function(arr_rel, score_lab, file_prefix) {
  for (ki in seq_along(prop_ks)) {
    pdf_file <- file.path(plots_dir,
      sprintf("%s_rel_gauss_by_quantile-%s-K%d.pdf",
              file_prefix, set_tag, prop_ks[ki]))
    pdf(pdf_file, width = 7, height = 5)
    mat <- arr_rel[, , ki]
    ymin <- min(mat, 1, na.rm = TRUE)
    ymax <- max(mat, 1, na.rm = TRUE)
    plot(probs, mat[, 1], type = "o", ylim = c(ymin, ymax),
         pch = mpch[1], lty = mlty[1], col = mcol[1], bg = mbg[1],
         xlab = "Threshold quantile",
         ylab = sprintf("Relative %s score (vs. Gaussian)", score_lab),
         main = sprintf("Setting %d%s, prop_k = %d",
                        setting_id, data_suffix, prop_ks[ki]))
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

# (2) mean score vs prop_k for selected quantiles, lines per method.
# Only meaningful when there is more than one prop_k; skip otherwise.
if (n_ks >= 2) {
  selected_q_idx <- c(1, 6, 9, 10)             # 0.90, 0.95, 0.98, 0.99
  selected_q_idx <- selected_q_idx[selected_q_idx <= n_probs]
  plot_mean_vs_K <- function(arr_mean, score_lab, file_prefix) {
    for (qi in selected_q_idx) {
      pdf_file <- file.path(plots_dir,
        sprintf("%s_mean_vs_K-%s-q%03d.pdf",
                file_prefix, set_tag, round(probs[qi] * 1000)))
      pdf(pdf_file, width = 7, height = 5)
      mat <- matrix(arr_mean[qi, , , drop = FALSE],
                    nrow = n_methods, ncol = n_ks)
      ymin <- min(mat, na.rm = TRUE)
      ymax <- max(mat, na.rm = TRUE)
      plot(prop_ks, mat[1, ], type = "o", ylim = c(ymin, ymax),
           pch = mpch[1], lty = mlty[1], col = mcol[1], bg = mbg[1],
           xlab = "prop basis rank K",
           ylab = sprintf("Mean %s score", score_lab),
           main = sprintf("Setting %d%s, q = %.3f",
                          setting_id, data_suffix, probs[qi]))
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
} else {
  cat("  (skipping mean-vs-K plots: only ", n_ks, " prop_k value)\n", sep = "")
}

# (3) lambda 95% CI vs dataset, one PDF per (skew method, prop_k)
if (length(ci_lo_idx) == 1L && length(ci_hi_idx) == 1L) {
  for (m in skew_methods) {
    mi <- match(m, methods)
    for (ki in seq_along(prop_ks)) {
      pdf_file <- file.path(plots_dir,
        sprintf("lambda_ci_vs_dataset-%s-method%d-K%d.pdf",
                set_tag, m, prop_ks[ki]))
      pdf(pdf_file, width = 8, height = 5)
      lo <- lambda[ci_lo_idx, , mi, ki]
      hi <- lambda[ci_hi_idx, , mi, ki]
      if (all(is.na(c(lo, hi)))) { dev.off(); file.remove(pdf_file); next }
      plot(NA, xlim = c(1, n_sets),
           ylim = range(c(lo, hi, lambda_true), na.rm = TRUE),
           xlab = "dataset", ylab = expression(lambda),
           main = sprintf("95%% CI for lambda - method %d, prop_k = %d",
                          m, prop_ks[ki]))
      abline(h = lambda_true, col = "red", lty = 2)
      segments(seq_len(n_sets), lo, seq_len(n_sets), hi,
               col = ifelse(lo <= lambda_true & hi >= lambda_true,
                            "black", "orange"))
      points(seq_len(n_sets), (lo + hi) / 2, pch = 19, cex = 0.4)
      dev.off()
    }
  }
}

cat("\nWrote plots to ", plots_dir, "\n", sep = "")
