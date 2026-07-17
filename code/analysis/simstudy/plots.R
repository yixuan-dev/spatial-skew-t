#########################################################################
# plots.R - Stage 3 of the simstudy post-fit pipeline.
#
# Reads output/results/simresults<setting><suffix>.RData (produced by
# tables.R) and renders the figures under output/plots/.
#
# Split out of tables.R so the table and figure stages re-run
# independently; mirrors code/analysis/simstudy_prop/plots-prop.R.
#
# Usage:
#   Rscript plots.R --setting=<id>
#                   [--data=<path>]
#
# Examples:
#   Rscript plots.R --setting=4
#   Rscript plots.R --setting=1 --data=simdata_def.RData
#
# Outputs: each figure lives in output/plots/<type>/, and <type> is also the
#  filename prefix, so file and folder share one descriptive name. <type> is one
#  of brier_score, quantile_score, energy_score, variogram_score,
#  predictive_rmse, recovery_rmse, lambda_ci.
#  (suffix = "" for simdata.RData, "_def" for simdata_def.RData, ...):
#   <type>/<type>_rel_gauss_by_quantile-set<setting><suffix>-K<k>.pdf   (brier_score, quantile_score)
#   <type>/<type>_mean_vs_K-set<setting><suffix>-q<qq>.pdf              (brier_score, quantile_score)
#   <type>/<type>_mean_vs_K-set<setting><suffix>.pdf                   (energy_score, variogram_score, predictive_rmse [omits methods 3 & 5], recovery_rmse)
#   lambda_ci/lambda_ci_vs_dataset-set<setting><suffix>-method<m>-K<k>.pdf
#########################################################################

rm(list = ls())

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0L) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]),
    winslash = "/", mustWork = FALSE
  )
  if (dir.exists(dirname(script_path))) setwd(dirname(script_path))
}

source("./helpers.R")

# ---- CLI parsing -----------------------------------------------------
cli_args <- commandArgs(trailingOnly = TRUE)
parsed <- extract_leading_flags(cli_args, c("data", "setting"))
flags <- parsed$values

if (is.null(flags$setting) || !nzchar(flags$setting)) {
  stop("plots.R: --setting=<id> is required.", call. = FALSE)
}
setting_id <- as.integer(parse_index_expr(flags$setting, "setting"))
if (length(setting_id) != 1L || setting_id < 1L) {
  stop("plots.R: --setting must be a single positive integer.", call. = FALSE)
}

data_suffix <- if (!is.null(flags$data) && nzchar(flags$data)) {
  derive_data_suffix(flags$data)
} else {
  ""
}

results_dir <- "output/results"
plots_dir <- "output/plots"
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
}

# Route each figure into a descriptive output/plots/<type>/ folder. The plot
# type (e.g. "brier_score") is used for BOTH the folder and the filename
# prefix, so files and their folder share the same descriptive name.
plot_pdf_path <- function(type, filename) {
  sub_dir <- file.path(plots_dir, type)
  if (!dir.exists(sub_dir)) {
    dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)
  }
  file.path(sub_dir, filename)
}

# ---- load the Stage-2 artifact ---------------------------------------
simresults_file <- file.path(
  results_dir,
  sprintf("simresults%d%s.RData", setting_id, data_suffix)
)
if (!file.exists(simresults_file)) {
  stop(sprintf(
    "Aggregated results not found: %s\n  Run tables.R --setting=%d%s first.",
    simresults_file, setting_id,
    if (nzchar(data_suffix)) sprintf(" --data=%s", flags$data) else ""
  ), call. = FALSE)
}
load(simresults_file)
# Provides: bs_mean/qs_mean/bs_med/qs_med, bs_rel_mean/qs_rel_mean/...,
#           score_summary_table, cov_table, lambda, probs, mrts_ks, methods,
#           datasets, intervals, setting, data_suffix; plus
#           energy.score / vario.score when the Stage-1 cache had them.

cat(sprintf(
  "plots: setting=%d cache=%s suffix='%s'\n",
  setting_id, simresults_file,
  if (nzchar(data_suffix)) data_suffix else "<none>"
))

n_probs <- length(probs)
n_sets <- length(datasets)
n_methods <- length(methods)
n_ks <- length(mrts_ks)

method_catalog <- get_simstudy_method_catalog(include_maxstable = TRUE)
method_label_for <- function(m) {
  row <- method_catalog[method_catalog$method_id == m, , drop = FALSE]
  if (nrow(row) == 1L) sprintf("%d: %s", m, row$label[1]) else as.character(m)
}
method_label <- vapply(methods, method_label_for, character(1))

has_multivar <- exists("energy.score") && exists("vario.score")
has_pred_rmse <- exists("pred.rmse")
has_recovery <- exists("recovery.rmse")

# lambda CI bookkeeping (mirrors tables.R's lambda_coverage block)
ci_lo_idx <- which(abs(intervals - 0.025) < 1e-12)
ci_hi_idx <- which(abs(intervals - 0.975) < 1e-12)
lambda_true <- 3
skew_methods <- intersect(c(2L, 4L, 7L, 8L, 9L, 10L), methods)

# ---- figures ---------------------------------------------------------
set_tag <- sprintf("set%d%s", setting_id, data_suffix)

# Data-generating process label per setting id (for figure titles).
# Stored as bquote() language objects so math symbols render through
# plot() main / mtext().
dgp_label_map <- list(
  "1" = bquote("Gaussian"),
  "2" = bquote(italic(t) * " (K=1)"),
  "3" = bquote(italic(t) * " (K=5)"),
  "4" = bquote("Skew-" * italic(t) * " (K=1, " * lambda * "=3)"),
  "5" = bquote("Skew-" * italic(t) * " (K=5, " * lambda * "=3)"),
  "6" = bquote("Max-stable, Reich and Shaby"),
  "7" = bquote("Transformed Skew-" * italic(t) * ", T=q(0.80)"),
  "8" = bquote("Max-stable, Brown-Resnick"),
  "9" = bquote("Skew-" * italic(t) * " (K=1, " * lambda * "=3), AR(2): " *
    phi[1] * "=0.8, " * phi[2] * "=-0.35"),
  "10" = bquote("Skew-" * italic(t) * " (K=1, " * lambda * "=3), AR(2): " *
    phi[1] * "=0.12, " * phi[2] * "=-0.05"),
  "11" = bquote("Skew-" * italic(t) * " (K=5, " * lambda * "=3), AR(2): " *
    phi[1] * "=0.8, " * phi[2] * "=-0.35"),
  "12" = bquote("Skew-" * italic(t) * " (K=5, " * lambda * "=3), AR(2): " *
    phi[1] * "=0.12, " * phi[2] * "=-0.05"),
  "13" = bquote("Skew-" * italic(t) * " (K=1, " * lambda * "=3), AR(2) on " *
    italic(z) * " only: " * phi[1] * "=0.8, " * phi[2] * "=-0.35"),
  "14" = bquote("Skew-" * italic(t) * " (K=1, " * lambda * "=3), AR(2) on " *
    tau * " only: " * phi[1] * "=0.8, " * phi[2] * "=-0.35"),
  "15" = bquote("Skew-" * italic(t) * " (K=5, " * lambda * "=3), AR(2) on " *
    italic(w) * " only: " * phi[1] * "=0.8, " * phi[2] * "=-0.35")
)
# The non-stationary-mean sim (simdata_nonsta.RData -> suffix "_nonsta") has
# its own DGP per setting; without this its settings would be mislabelled
# with the simdata catalog (e.g. setting 1 -> "Gaussian").
#
# Settings 1-2 put the non-stationarity in the MEAN (MRTS can recover it),
# 3-8 in the DEPENDENCE (it structurally cannot), 9-10 in the higher moments.
# Settings 1 and 6 are the headline pair: the SAME cosine field f1 drives the
# mean in 1 and the correlation range in 6. See setup_nonsta.R.
dgp_label_map_nonsta <- list(
  "1"  = bquote("Skew-" * italic(t) * " (K=1, " * lambda * "=3), non-stationary MEAN: fixed cosine-bump surface"),
  "2"  = bquote("Skew-" * italic(t) * " (K=1, " * lambda * "=3), non-stationary MEAN: time-varying cosine-bump surface"),
  "3"  = bquote("Skew-" * italic(t) * " (K=1, " * lambda * "=3), non-stationary DEP.: low-rank cosine random effect"),
  "4"  = bquote("Skew-" * italic(t) * " (K=1, " * lambda * "=3), non-stationary DEP.: Paciorek-Schervish local anisotropy (in " * italic(C) * ")"),
  "5"  = bquote("Skew-" * italic(t) * " (K=1, " * lambda * "=3), non-stationary DEP.: Paciorek-Schervish field (additive)"),
  "6"  = bquote("Skew-" * italic(t) * " (K=1, " * lambda * "=3), non-stationary DEP.: covariate-driven range " * rho * "(" * italic(s) * ")" ~ "(in " * italic(C) * ")"),
  "7"  = bquote("Skew-" * italic(t) * " (K=1, " * lambda * "=3), non-stationary DEP.: covariate-driven range (additive)"),
  "8"  = bquote("Skew-" * italic(t) * " (K=1, " * lambda * "=3), non-stationary DEP.: Fuentes 4-regime GP mixture (in " * italic(C) * ")"),
  "9"  = bquote("Skew-" * italic(t) * " (K=1, " * lambda * "=3), non-stationary TAILS: Tukey " * italic(g) * "-and-" * italic(h) * " field"),
  "10" = bquote("Skew-" * italic(t) * " (K=1), non-stationary SKEWNESS: " * lambda * "(" * italic(s) * ") = 3 + 3" * italic(f)[2] * "(" * italic(s) * ")"),
  "11" = bquote("Skew-" * italic(t) * " (K=1, " * lambda * "=3), non-stationary MEAN, HIGH SNR (3" * symbol("\264") * " setting 1): Brier positive control"),
  "12" = bquote("Skew-" * italic(t) * " (K=1, " * lambda * "=3), RICH MEAN (6 bumps, sd 11): MRTS large-K showcase")
)
# The deformed-covariance sim (simdata_def.RData -> suffix "_def") is all
# Skew-t (K=1, lambda=3) with a deformed covariance per setting (see
# setup_def.R); without this its settings inherit the simdata catalog
# (1 -> "Gaussian"), which is wrong.
dgp_label_map_def <- list(
  "1" = bquote("Skew-" * italic(t) * " (K=1, " * lambda * "=3), isotropic (stationary baseline)"),
  "2" = bquote("Skew-" * italic(t) * " (K=1, " * lambda * "=3), geometric anisotropy (" * theta * "=" * pi * "/4, r=0.5)"),
  "3" = bquote("Skew-" * italic(t) * " (K=1, " * lambda * "=3), geometric anisotropy (" * theta * "=" * pi * "/6, r=0.25)"),
  "4" = bquote("Skew-" * italic(t) * " (K=1, " * lambda * "=3), non-stationary cov. (axial deformation)"),
  "5" = bquote("Skew-" * italic(t) * " (K=1, " * lambda * "=3), non-stationary cov. (sine deformation)"),
  "6" = bquote("Skew-" * italic(t) * " (K=1, " * lambda * "=3), non-stationary cov. (composed deformation)")
)

# dgp_title(): expression built from setting_id + data_suffix,
# optionally suffixed by an extra `bquote(...)` snippet for context.
dgp_title <- function(suffix_expr = NULL) {
  if (identical(data_suffix, "_nonsta")) {
    base <- dgp_label_map_nonsta[[as.character(setting_id)]]
    if (is.null(base)) base <- bquote(paste("nonsta setting ", .(setting_id)))
  } else if (identical(data_suffix, "_def")) {
    base <- dgp_label_map_def[[as.character(setting_id)]]
    if (is.null(base)) base <- bquote(paste("deformed setting ", .(setting_id)))
  } else {
    base <- dgp_label_map[[as.character(setting_id)]]
    if (is.null(base)) base <- bquote(paste("setting ", .(setting_id)))
    if (nzchar(data_suffix)) {
      base <- bquote(.(base) ~ .(paste0("[", data_suffix, "]")))
    }
  }
  if (is.null(suffix_expr)) base else bquote(.(base) * ", " * .(suffix_expr))
}

# One entry per Morris method id 1..10 (palette[methods[j]], not palette[j]).
# Systematic encoding: colour = family (Skew-t red / Sym-t blue),
# lty = knots (K=1 solid / K=5 dotted), pch = base vs temporal (baseline
# filled square; AR(2) "+"; AR(1) "x" -- open marks so the baseline square
# underneath is not obscured). methods 1 (Gaussian) / 6 (max-stable) are
# outside the family x K scheme; 6 uses a filled diamond so "x" is free for
# AR(1).
mlty <- c(1, 1, 1, 3, 3, 6, 1, 3, 1, 3)
mpch <- c(21, 22, 22, 22, 22, 23, 3, 3, 4, 4)
mcol <- c(
  "gray30", "firebrick4", "dodgerblue4", "firebrick4", "dodgerblue4",
  "darkgreen", "firebrick4", "firebrick4", "firebrick4", "firebrick4"
)
mbg <- c(
  "gray70", "firebrick2", "dodgerblue2", "firebrick2", "dodgerblue2",
  "lightgreen", "firebrick2", "firebrick2", "firebrick2", "firebrick2"
)
# Style vectors are indexed positionally by method_id (mlty[methods[j]]). Guard
# against an out-of-range id silently falling through to NA (default styling).
stopifnot(all(methods >= 1L & methods <= length(mlty)))

# Shared drawing helpers for the per-method line plots (blocks 1-3 below).
# add_method_series(): draw the first series with plot() then overlay the rest
# with lines()+points(); style each by its method_id via the m* vectors. `ymat`
# is oriented one column per method (callers with method-as-row pass t(mat)).
# `idx` selects which method positions to draw; `...` is forwarded to plot().
add_method_series <- function(x, ymat, idx, ...) {
  m1 <- methods[idx[1]]
  plot(x, ymat[, idx[1]],
    type = "o",
    pch = mpch[m1], lty = mlty[m1], col = mcol[m1], bg = mbg[m1], ...
  )
  for (j in idx[-1]) {
    mj <- methods[j]
    lines(x, ymat[, j], lty = mlty[mj], col = mcol[mj])
    points(x, ymat[, j], pch = mpch[mj], col = mcol[mj], bg = mbg[mj])
  }
}
add_method_legend <- function(pos, idx = seq_len(n_methods)) {
  legend(pos,
    legend = method_label[idx],
    lty = mlty[methods[idx]], pch = mpch[methods[idx]],
    col = mcol[methods[idx]], pt.bg = mbg[methods[idx]],
    cex = 0.7, bty = "n"
  )
}

# (1) relative score vs quantile, lines per method, one PDF per mrts_k
plot_rel_vs_quantile <- function(arr_rel, score_lab, file_prefix) {
  for (ki in seq_along(mrts_ks)) {
    # Force 2-D so the n_methods == 1 slice does not drop to a vector.
    mat <- matrix(arr_rel[, , ki], nrow = n_probs, ncol = n_methods)
    if (all(is.na(mat))) next # e.g. no Gaussian (method 1) -> rel undefined
    pdf_file <- plot_pdf_path(
      file_prefix,
      sprintf(
        "%s_rel_gauss_by_quantile-%s-K%d.pdf",
        file_prefix, set_tag, mrts_ks[ki]
      )
    )
    pdf(pdf_file, width = 7, height = 5)
    ymin <- min(mat, 1, na.rm = TRUE)
    ymax <- max(mat, 1, na.rm = TRUE)
    add_method_series(probs, mat, seq_len(n_methods),
      ylim = c(ymin, ymax),
      xlab = "Threshold quantile",
      ylab = sprintf("Relative %s score", score_lab),
      main = dgp_title()
    )
    abline(h = 1, lty = 2, col = "gray60")
    add_method_legend("topleft")
    dev.off()
  }
}
plot_rel_vs_quantile(bs_rel_mean, "Brier", "brier_score")
plot_rel_vs_quantile(qs_rel_mean, "Quantile", "quantile_score")

# (2) mean score vs mrts_k for selected quantiles, lines per method.
# Only meaningful when there is more than one mrts_k; skip otherwise.
if (n_ks >= 2) {
  # Select target quantiles by value (not position) so the plot stays correct
  # if `probs` changes; tolerance match mirrors the lambda CI index lookup above.
  target_q <- c(0.90, 0.95, 0.98, 0.99)
  selected_q_idx <- which(vapply(
    probs, function(p) any(abs(p - target_q) < 1e-9), logical(1)
  ))
  plot_mean_vs_K <- function(arr_mean, score_lab, file_prefix) {
    for (qi in selected_q_idx) {
      pdf_file <- plot_pdf_path(
        file_prefix,
        sprintf(
          "%s_mean_vs_K-%s-q%03d.pdf",
          file_prefix, set_tag, round(probs[qi] * 1000)
        )
      )
      pdf(pdf_file, width = 7, height = 5)
      mat <- matrix(arr_mean[qi, , , drop = FALSE],
        nrow = n_methods, ncol = n_ks
      )
      ymin <- min(mat, na.rm = TRUE)
      ymax <- max(mat, na.rm = TRUE)
      add_method_series(mrts_ks, t(mat), seq_len(n_methods),
        ylim = c(ymin, ymax),
        xlab = "Number of MRTS basis functions",
        ylab = bquote(.(sprintf("%s score", score_lab)) * "," ~
          q == .(sprintf("%.3f", probs[qi]))),
        main = dgp_title()
      )
      add_method_legend("topright")
      dev.off()
    }
  }
  plot_mean_vs_K(bs_mean, "Brier", "brier_score")
  plot_mean_vs_K(qs_mean, "Quantile", "quantile_score")
} else {
  cat("  (skipping mean-vs-K plots: only ", n_ks, " mrts_k value)\n", sep = "")
}

# (3) energy / variogram score vs mrts_k, lines per method.
# One number per (method, mrts_k), so this only varies along K.
if ((has_multivar || has_pred_rmse || has_recovery) && n_ks >= 2) {
  plot_multivar_vs_K <- function(arr, score_lab, file_prefix,
                                 ylab = sprintf("Mean %s score", score_lab),
                                 drop_methods = integer(0)) {
    mat <- apply(arr, c(2, 3), mean, na.rm = TRUE) # method x mrts_k
    keep <- which(!(methods %in% drop_methods))
    if (length(keep) == 0L || all(is.na(mat[keep, , drop = FALSE]))) {
      return(invisible(NULL))
    }
    pdf_file <- plot_pdf_path(
      file_prefix,
      sprintf("%s_mean_vs_K-%s.pdf", file_prefix, set_tag)
    )
    pdf(pdf_file, width = 7, height = 5)
    ymin <- min(mat[keep, ], na.rm = TRUE)
    ymax <- max(mat[keep, ], na.rm = TRUE)
    add_method_series(mrts_ks, t(mat), keep,
      ylim = c(ymin, ymax),
      xlab = "Number of MRTS basis functions",
      ylab = ylab,
      main = dgp_title()
    )
    add_method_legend("topright", keep)
    dev.off()
  }
  if (has_multivar) {
    plot_multivar_vs_K(energy.score, "energy", "energy_score")
    plot_multivar_vs_K(vario.score, "variogram", "variogram_score")
  }
  if (has_pred_rmse) {
    # point-prediction RMSE (vs observed y); omit the symmetric-t outliers
    # (methods 3 & 5) that otherwise dominate the y-axis scale
    plot_multivar_vs_K(pred.rmse, "predictive RMSE", "predictive_rmse",
      ylab = "Predictive RMSE",
      drop_methods = c(3, 5)
    )
  }
  if (has_recovery) {
    # mean-surface recovery RMSE (vs true mean); the mean-estimation channel
    plot_multivar_vs_K(recovery.rmse, "recovery RMSE", "recovery_rmse",
      ylab = "Mean-surface recovery RMSE  (test sites)"
    )
  }
} else if (!has_multivar && !has_pred_rmse && !has_recovery) {
  cat("  (skipping energy/variogram/pred-rmse/recovery plots: no such scores)\n")
} else {
  # scores exist but only one mrts_k -> nothing varies along K
  cat("  (skipping energy/variogram/pred-rmse/recovery plots: only ",
    n_ks, " mrts_k value)\n",
    sep = ""
  )
}

# (4) lambda 95% CI vs dataset, one PDF per (skew method, mrts_k)
if (length(ci_lo_idx) == 1L && length(ci_hi_idx) == 1L) {
  for (m in skew_methods) {
    mi <- match(m, methods)
    for (ki in seq_along(mrts_ks)) {
      pdf_file <- plot_pdf_path(
        "lambda_ci",
        sprintf(
          "lambda_ci_vs_dataset-%s-method%d-K%d.pdf",
          set_tag, m, mrts_ks[ki]
        )
      )
      lo <- lambda[ci_lo_idx, , mi, ki]
      hi <- lambda[ci_hi_idx, , mi, ki]
      # Skip before opening the device so no empty PDF is written then removed.
      if (all(is.na(c(lo, hi)))) next
      pdf(pdf_file, width = 8, height = 5)
      plot(NA,
        xlim = c(1, n_sets),
        ylim = range(c(lo, hi, lambda_true), na.rm = TRUE),
        xlab = "dataset", ylab = expression(lambda),
        main = sprintf(
          "95%% CI for lambda - method %d, mrts_k = %d",
          m, mrts_ks[ki]
        )
      )
      abline(h = lambda_true, col = "red", lty = 2)
      segments(seq_len(n_sets), lo, seq_len(n_sets), hi,
        col = ifelse(lo <= lambda_true & hi >= lambda_true,
          "black", "orange"
        )
      )
      points(seq_len(n_sets), (lo + hi) / 2, pch = 19, cex = 0.4)
      dev.off()
    }
  }
}

cat("\nWrote plots to ", plots_dir, "\n", sep = "")
