#########################################################################
# plots_for_paper.R - paper-ready combined figures from simresult artifacts.
#
# Merges per-setting brier_score plots (normally one PDF per setting from
# plots.R) into single multi-panel figures for the thesis. Reads
# output/results/simresults<set><suffix>.RData (produced by tables.R) and writes
# the combined PDFs straight into the thesis figure dir (myLatex/pdf).
#
# Shared paper-format conventions (all figures): 3 rows x 2 cols layout; four
# (a)-(d) panels in regions 1-4; region 5 (row 3, col 1) left blank; a single
# shared legend in region 6 (row 3, col 2) as one vertical column with bare
# method labels; enlarged fonts; no per-panel legends.
#
# Figure 1 - relative Brier score vs threshold quantile (AR(2) settings 9-12).
#   Relative to the Gaussian model, so Gaussian (method 1) is the dashed y = 1
#   reference and is omitted from the curves/legend.
#   (a) set9  Skew-t K=1, lambda=3, AR(2) strong phi=(0.80,-0.35)
#   (b) set10 Skew-t K=1, lambda=3, AR(2) weak   phi=(0.12,-0.05)
#   (c) set11 Skew-t K=5, lambda=3, AR(2) strong phi=(0.80,-0.35)
#   (d) set12 Skew-t K=5, lambda=3, AR(2) weak   phi=(0.12,-0.05)
#
# Figure 2 - mean Brier score vs number of MRTS basis functions (MRTS probe).
#   Absolute score, so Gaussian (method 1) is kept as a real curve.
#   (a) set5         q=0.95   Skew-t K=5, lambda=3
#   (b) set5         q=0.98
#   (c) set3_nonsta  q=0.95   Skew-t K=1, lambda=3, additive non-stationary spatial random effect
#   (d) set3_nonsta  q=0.98
#
# Usage (from code/analysis/simstudy/, via PowerShell Rscript.exe):
#   Rscript.exe plots_for_paper.R
#
# Outputs (myLatex/pdf):
#   bs_rel_gauss_by_quantile-set9_to_12-K0.pdf
#   brier_score_mean_vs_K.pdf
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

results_dir <- "output/results"
out_dir <- "../../../myLatex/pdf"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# ---- shared method styling + drawing helpers (mirror plots.R) ---------
# One entry per Morris method id 1..10 (vectors indexed by method_id).
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

# Bare method labels (no "id:" prefix) for the paper legend.
method_catalog <- get_simstudy_method_catalog(include_maxstable = TRUE)
label_for_id <- function(m) {
  row <- method_catalog[method_catalog$method_id == m, , drop = FALSE]
  if (nrow(row) == 1L) row$label[1] else as.character(m)
}

# draw_series(): plot the first series, overlay the rest. Columns of `ymat`
# align positionally with `method_ids`; each series is styled by its method id.
draw_series <- function(x, ymat, method_ids, ...) {
  m1 <- method_ids[1]
  plot(x, ymat[, 1],
    type = "o",
    pch = mpch[m1], lty = mlty[m1], col = mcol[m1], bg = mbg[m1], ...
  )
  for (j in seq_along(method_ids)[-1]) {
    m <- method_ids[j]
    lines(x, ymat[, j], lty = mlty[m], col = mcol[m])
    points(x, ymat[, j], pch = mpch[m], col = mcol[m], bg = mbg[m])
  }
}

# draw_legend_region(): region 5 blank, then the shared vertical legend in
# region 6. Call once, after the four panels, with the union of drawn ids.
draw_legend_region <- function(method_ids) {
  par(mar = c(0, 0, 0, 0))
  plot.new() # region 5 (row 3, col 1): intentionally blank
  plot.new() # region 6 (row 3, col 2): shared legend
  legend("center",
    legend = vapply(method_ids, label_for_id, character(1)),
    lty = mlty[method_ids], pch = mpch[method_ids],
    col = mcol[method_ids], pt.bg = mbg[method_ids],
    ncol = 1, cex = 1.25, bty = "n"
  )
}

open_panel_device <- function(out_file) {
  pdf(out_file, width = 9, height = 11)
  layout(matrix(1:6, nrow = 3, byrow = TRUE), heights = c(1, 1, 0.8))
  par(mar = c(4.6, 4.8, 2.4, 1.0), cex.axis = 1.1, cex.lab = 1.35, cex.main = 1.5)
}

load_sim <- function(setting_id, suffix = "") {
  f <- file.path(results_dir, sprintf("simresults%d%s.RData", setting_id, suffix))
  if (!file.exists(f)) {
    stop(sprintf("Aggregated results not found: %s (run tables.R first)", f),
      call. = FALSE
    )
  }
  e <- new.env()
  load(f, envir = e)
  e
}

# ---- Figure 1: relative Brier vs quantile, AR(2) settings 9-12 --------
make_fig_rel_by_quantile <- function() {
  sets <- c(9L, 10L, 11L, 12L)
  tags <- c("(a)", "(b)", "(c)", "(d)")
  k0 <- 0L
  envs <- lapply(sets, load_sim)
  probs <- envs[[1]]$probs
  methods <- envs[[1]]$methods
  for (e in envs) {
    stopifnot(
      identical(e$methods, methods),
      isTRUE(all.equal(e$probs, probs)),
      k0 %in% e$mrts_ks
    )
  }
  # Gaussian (method 1) is the relative reference (== 1, shown as dashed line).
  draw_cols <- which(methods != 1L)
  draw_ids <- methods[draw_cols]

  out_file <- file.path(out_dir, "bs_rel_gauss_by_quantile-set9_to_12-K0.pdf")
  open_panel_device(out_file)
  for (i in seq_along(sets)) {
    e <- envs[[i]]
    ki <- which(e$mrts_ks == k0)
    mat <- matrix(e$bs_rel_mean[, , ki], nrow = length(probs), ncol = length(methods))
    ymat <- mat[, draw_cols, drop = FALSE]
    # Range over the drawn (non-Gaussian) series; keep 1 in view for the ref line.
    draw_series(probs, ymat, draw_ids,
      ylim = c(min(ymat, 1, na.rm = TRUE), max(ymat, 1, na.rm = TRUE)),
      xlab = "Threshold quantile",
      ylab = "Relative Brier score",
      main = tags[i]
    )
    abline(h = 1, lty = 2, col = "gray60")
  }
  draw_legend_region(draw_ids)
  dev.off()
  cat(sprintf("Wrote %s\n", normalizePath(out_file, winslash = "/", mustWork = FALSE)))
}

# ---- Figure 2: mean Brier vs MRTS K, set5 + set3_nonsta ---------------
make_fig_mean_vs_K <- function() {
  panels <- list(
    list(setting = 5L, suffix = "", q = 0.95, tag = "(a)"),
    list(setting = 5L, suffix = "", q = 0.98, tag = "(b)"),
    list(setting = 3L, suffix = "_nonsta", q = 0.95, tag = "(c)"),
    list(setting = 3L, suffix = "_nonsta", q = 0.98, tag = "(d)")
  )
  out_file <- file.path(out_dir, "brier_score_mean_vs_K.pdf")
  open_panel_device(out_file)
  legend_ids <- integer(0)
  for (p in panels) {
    e <- load_sim(p$setting, p$suffix)
    probs <- e$probs
    methods <- e$methods
    mrts_ks <- e$mrts_ks
    qi <- which(abs(probs - p$q) < 1e-9)
    stopifnot(length(qi) == 1L, length(mrts_ks) >= 2L)
    # method x mrts_k, then transpose so columns align with method_ids.
    mat <- matrix(e$bs_mean[qi, , ], nrow = length(methods), ncol = length(mrts_ks))
    draw_series(mrts_ks, t(mat), methods,
      ylim = range(mat, na.rm = TRUE),
      xlab = "Number of MRTS basis functions",
      ylab = bquote("Brier score," ~ q == .(sprintf("%.3f", p$q))),
      main = p$tag
    )
    legend_ids <- union(legend_ids, methods)
  }
  draw_legend_region(sort(legend_ids))
  dev.off()
  cat(sprintf("Wrote %s\n", normalizePath(out_file, winslash = "/", mustWork = FALSE)))
}

make_fig_rel_by_quantile()
make_fig_mean_vs_K()
