#########################################################################
# tables-prop.R - Stage 2 of the simstudy_prop post-fit pipeline.
#
# Reads output/results/scores<setting>-prop<suffix>.RData (produced by
# scores-prop.R) and emits CSV tables + an aggregated .RData artifact
# under output/tables/ and output/results/.
#
# This script supersedes analyze4-prop.R and the table-emission halves of
# analyze_mrts.R / results-prop.R.
#
# Usage:
#   Rscript tables-prop.R --setting=<id>
#                         [--data=<path>]
#                         [--no-plots]
#
# Examples:
#   Rscript tables-prop.R --setting=4
#   Rscript tables-prop.R --setting=1 --data=simdata_def.RData
#   Rscript tables-prop.R --setting=4 --no-plots
#
# Outputs (suffix = "" for simdata.RData, "_def" for simdata_def.RData, ...):
#   output/tables/score_long<setting>-prop<suffix>.csv
#   output/tables/score_mean<setting>-prop<suffix>.csv
#   output/tables/score_rel_gauss<setting>-prop<suffix>.csv
#   output/tables/best_method_per_K<setting>-prop<suffix>.csv
#   output/tables/lambda_coverage<setting>-prop<suffix>.csv
#   output/results/simresults<setting>-prop<suffix>.RData
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
no_plots <- any(parsed$args == "--no-plots") ||
            any(cli_args == "--no-plots")

if (is.null(flags$setting) || !nzchar(flags$setting)) {
  stop("tables-prop.R: --setting=<id> is required.", call. = FALSE)
}
setting_id <- as.integer(parse_index_expr(flags$setting, "setting"))
if (length(setting_id) != 1L || setting_id < 1L) {
  stop("tables-prop.R: --setting must be a single positive integer.", call. = FALSE)
}

# ---- locate stage-1 cache --------------------------------------------
data_suffix <- if (!is.null(flags$data) && nzchar(flags$data)) {
  derive_data_suffix(flags$data)
} else {
  ""
}

results_dir <- "output/results"
tables_dir  <- "output/tables"
plots_dir   <- "output/plots"
for (d in c(results_dir, tables_dir, plots_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

scores_file <- file.path(
  results_dir,
  sprintf("scores%d-prop%s.RData", setting_id, data_suffix)
)
if (!file.exists(scores_file)) {
  stop(sprintf(
    "Score cache not found: %s\n  Run scores-prop.R --setting=%d%s first.",
    scores_file, setting_id,
    if (nzchar(data_suffix)) sprintf(" --data=%s", flags$data) else ""
  ), call. = FALSE)
}
load(scores_file)
# Provides: quant.score, brier.score, beta.0/1/2, tau.alpha, tau.beta,
#           rho, nu, gamma, lambda, elapsed_sec, probs, intervals,
#           prop_ks, datasets, methods, setting,
#           data_path, data_suffix, fits_dir   (provenance)

cat(sprintf(
  "tables-prop: setting=%d cache=%s suffix='%s' plots=%s\n",
  setting_id, scores_file,
  if (nzchar(data_suffix)) data_suffix else "<none>",
  if (no_plots) "off" else "on"
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

# ---- mean / median scores ------------------------------------------
mean_apply <- function(a) apply(a, c(1, 3, 4), mean,   na.rm = TRUE)
med_apply  <- function(a) apply(a, c(1, 3, 4), median, na.rm = TRUE)

bs_mean <- mean_apply(brier.score)   # [probs, method, prop_k]
qs_mean <- mean_apply(quant.score)
bs_med  <- med_apply(brier.score)
qs_med  <- med_apply(quant.score)

dn_mean <- list(quantile = as.character(probs),
                method   = as.character(methods),
                prop_k   = as.character(prop_ks))
dimnames(bs_mean) <- dimnames(qs_mean) <- dn_mean
dimnames(bs_med)  <- dimnames(qs_med)  <- dn_mean

# ---- score_long: per-(score, method, prop_k, dataset, quantile) ------
# Replaces analyze_mrts.R's mrts_brier_long.csv. One row per cell of the
# 4-D score arrays so per-dataset slicing remains possible.
flatten_long <- function(arr, score) {
  grid <- expand.grid(
    quantile = probs,
    dataset  = datasets,
    method   = methods,
    prop_k   = prop_ks,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$score <- score
  grid$value <- as.vector(arr)
  grid[, c("score", "method", "prop_k", "dataset", "quantile", "value")]
}
score_long_table <- rbind(
  flatten_long(brier.score, "brier"),
  flatten_long(quant.score, "quant")
)
write.csv(
  score_long_table,
  file.path(tables_dir, sprintf("score_long%d-prop%s.csv", setting_id, data_suffix)),
  row.names = FALSE
)

# ---- score_mean: aggregated per (score, method, prop_k, quantile) ----
score_long <- function(arr_mean, arr_med, score) {
  out <- expand.grid(quantile = probs, method = methods, prop_k = prop_ks,
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  out$score  <- score
  out$mean   <- as.vector(arr_mean)
  out$median <- as.vector(arr_med)
  out
}
score_table <- rbind(
  score_long(bs_mean, bs_med, "brier"),
  score_long(qs_mean, qs_med, "quant")
)
write.csv(
  score_table[, c("score", "method", "prop_k", "quantile", "mean", "median")],
  file.path(tables_dir, sprintf("score_mean%d-prop%s.csv", setting_id, data_suffix)),
  row.names = FALSE
)

# ---- score_rel_gauss: relative scores vs method 1 (Gaussian) ---------
rel_to_gauss <- function(arr) {
  if (!"1" %in% dimnames(arr)$method) {
    return(array(NA_real_, dim = dim(arr), dimnames = dimnames(arr)))
  }
  g <- arr[, "1", , drop = FALSE]
  sweep(arr, c(1, 3),
        array(g, dim = c(dim(arr)[1], 1, dim(arr)[3])),
        "/")
}
bs_rel_mean <- rel_to_gauss(bs_mean)
qs_rel_mean <- rel_to_gauss(qs_mean)
bs_rel_med  <- rel_to_gauss(bs_med)
qs_rel_med  <- rel_to_gauss(qs_med)

rel_long <- function(arr_mean, arr_med, score) {
  out <- expand.grid(quantile = probs, method = methods, prop_k = prop_ks,
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  out$score      <- score
  out$rel_mean   <- as.vector(arr_mean)
  out$rel_median <- as.vector(arr_med)
  out
}
rel_table <- rbind(
  rel_long(bs_rel_mean, bs_rel_med, "brier"),
  rel_long(qs_rel_mean, qs_rel_med, "quant")
)
write.csv(
  rel_table[, c("score", "method", "prop_k", "quantile", "rel_mean", "rel_median")],
  file.path(tables_dir, sprintf("score_rel_gauss%d-prop%s.csv", setting_id, data_suffix)),
  row.names = FALSE
)

# ---- best_method_per_K: lowest mean per quantile ---------------------
best_combo <- function(arr_mean, score) {
  out <- data.frame(score = score, quantile = probs,
                    best_method = NA_integer_,
                    best_prop_k = NA_integer_,
                    best_value  = NA_real_)
  for (i in seq_len(n_probs)) {
    # Force 2-D shape so the n_ks == 1 (or n_methods == 1) case still works.
    mat <- matrix(arr_mean[i, , , drop = FALSE],
                  nrow = n_methods, ncol = n_ks)
    if (all(is.na(mat))) next
    idx <- which(mat == min(mat, na.rm = TRUE), arr.ind = TRUE)[1, ]
    out$best_method[i] <- methods[idx[1]]
    out$best_prop_k[i] <- prop_ks[idx[2]]
    out$best_value[i]  <- mat[idx[1], idx[2]]
  }
  out
}
best_table <- rbind(best_combo(bs_mean, "brier"),
                    best_combo(qs_mean, "quant"))
write.csv(
  best_table,
  file.path(tables_dir, sprintf("best_method_per_K%d-prop%s.csv", setting_id, data_suffix)),
  row.names = FALSE
)

# ---- lambda_coverage: 95% interval coverage for skew methods ---------
ci_lo_idx <- which(abs(intervals - 0.025) < 1e-12)
ci_hi_idx <- which(abs(intervals - 0.975) < 1e-12)
lambda_true  <- 3
skew_methods <- intersect(c(2L, 4L), methods)
cov_rows <- list()
if (length(ci_lo_idx) == 1L && length(ci_hi_idx) == 1L) {
  for (m in skew_methods) {
    mi <- match(m, methods)
    for (ki in seq_along(prop_ks)) {
      lo <- lambda[ci_lo_idx, , mi, ki]
      hi <- lambda[ci_hi_idx, , mi, ki]
      cover <- mean(lo <= lambda_true & hi >= lambda_true, na.rm = TRUE)
      width <- mean(hi - lo, na.rm = TRUE)
      cov_rows[[length(cov_rows) + 1L]] <- data.frame(
        method        = m,
        prop_k        = prop_ks[ki],
        coverage_95   = cover,
        mean_width_95 = width,
        n             = sum(!is.na(lo) & !is.na(hi))
      )
    }
  }
}
cov_table <- if (length(cov_rows) > 0L) {
  do.call(rbind, cov_rows)
} else {
  data.frame(
    method        = integer(0),
    prop_k        = integer(0),
    coverage_95   = numeric(0),
    mean_width_95 = numeric(0),
    n             = integer(0)
  )
}
write.csv(
  cov_table,
  file.path(tables_dir, sprintf("lambda_coverage%d-prop%s.csv", setting_id, data_suffix)),
  row.names = FALSE
)

# ---- save aggregated .RData artifact --------------------------------
simresults_file <- file.path(
  results_dir,
  sprintf("simresults%d-prop%s.RData", setting_id, data_suffix)
)
save(
  bs_mean, qs_mean, bs_med, qs_med,
  bs_rel_mean, qs_rel_mean, bs_rel_med, qs_rel_med,
  score_long_table, score_table, rel_table, best_table, cov_table,
  probs, prop_ks, methods, datasets, intervals, setting,
  data_suffix,
  file = simresults_file
)

# ---- plots (default-on; --no-plots to skip) -------------------------
if (!no_plots) {
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
}

cat("\nWrote tables to ", tables_dir, "\n", sep = "")
if (!no_plots) cat("Wrote plots  to ", plots_dir,  "\n", sep = "")
cat("Saved analysis objects to ", simresults_file, "\n", sep = "")
