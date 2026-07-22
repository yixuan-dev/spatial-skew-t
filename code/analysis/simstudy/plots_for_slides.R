#########################################################################
# plots_for_slides.R - beamer-format variants of two thesis figures.
#
# The thesis figures from plots_for_paper.R are portrait multi-panel
# pages sized for A4 floats; on the 16:9 defense slides their panels
# shrink beyond legibility.  This script redraws only the panels the
# talk narrates, in a wide single-row layout with enlarged fonts, and
# writes them next to the thesis figures in myLatex/pdf.  The full
# grids remain on the backup slides via the plots_for_paper.R PDFs.
# Reads the same aggregated artifacts
# output/results/simresults<set><suffix>.RData; styling mirrors
# plots_for_paper.R.
#
# Slide figure A - relative Brier score vs threshold quantile, spatial
#   hold-out evaluation, three settings in one row (thesis Figure 1
#   panels (a), (c), (e)):
#   set9   Skew-t K=1, lambda=3, AR(2) strong  phi=(0.80,-0.35)
#   set11  Skew-t K=5, lambda=3, AR(2) strong  phi=(0.80,-0.35)
#   set16  Skew-t K=1, lambda=3, AR(2)         phi=(0.15, 0.80)
# Slide figure B - mean Brier score vs number of MRTS basis functions
#   at q=0.95, the two nonlinear-mean DGPs in one row (thesis Figure 4
#   panels (a), (c)):
#   set16_nonsta  tanh arc front
#   set20_nonsta  log/exp/ratio transform mix
#
# Usage (from code/analysis/simstudy/, via PowerShell Rscript.exe):
#   Rscript.exe plots_for_slides.R
#
# Outputs (myLatex/pdf):
#   slides-bs_rel_gauss-set9_11_16.pdf
#   slides-brier_vs_K-set16_20_nonsta-q95.pdf
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

# ---- method styling (mirrors plots_for_paper.R) -----------------------
mlty <- c(1, 1, 1, 3, 3, 6, 1, 3, 1, 3)
mpch <- c(21, 22, 22, 22, 22, 23, 3, 3, 4, 4)
ar2_ids <- c(7L, 8L)
style_legacy <- list(
  col = c(
    "gray30", "firebrick4", "dodgerblue4", "firebrick4", "dodgerblue4",
    "darkgreen", "firebrick4", "firebrick4", "firebrick4", "firebrick4"
  ),
  bg = c(
    "gray70", "firebrick2", "dodgerblue2", "firebrick2", "dodgerblue2",
    "lightgreen", "firebrick2", "firebrick2", "firebrick2", "firebrick2"
  ),
  lwd = rep(1, 10)
)
style_ar2 <- list(
  col = c(
    "gray55", "gray25", "slategray", "gray25", "slategray",
    "gray70", "firebrick3", "firebrick3", "dodgerblue3", "dodgerblue3"
  ),
  bg = c(
    "gray80", "gray50", "lightsteelblue1", "gray50", "lightsteelblue1",
    "gray85", "firebrick1", "firebrick1", "dodgerblue1", "dodgerblue1"
  ),
  lwd = ifelse(seq_len(10) %in% ar2_ids, 2.5, 1)
)

method_catalog <- get_simstudy_method_catalog(include_maxstable = TRUE)
label_for_id <- function(m) {
  row <- method_catalog[method_catalog$method_id == m, , drop = FALSE]
  if (nrow(row) == 1L) row$label[1] else as.character(m)
}

# Slide variant of draw_series(): same series order (AR(2) last, on top),
# with point size and line width scaled up for projection.
pt_cex <- 1.3
lwd_mult <- 1.6
draw_series_slide <- function(x, ymat, method_ids, style, ...) {
  ord <- order(method_ids %in% ar2_ids)
  ymat <- ymat[, ord, drop = FALSE]
  method_ids <- method_ids[ord]
  m1 <- method_ids[1]
  plot(x, ymat[, 1],
    type = "o", cex = pt_cex,
    pch = mpch[m1], lty = mlty[m1], col = style$col[m1], bg = style$bg[m1],
    lwd = style$lwd[m1] * lwd_mult, ...
  )
  for (j in seq_along(method_ids)[-1]) {
    m <- method_ids[j]
    lines(x, ymat[, j],
      lty = mlty[m], col = style$col[m], lwd = style$lwd[m] * lwd_mult
    )
    points(x, ymat[, j],
      pch = mpch[m], col = style$col[m], bg = style$bg[m],
      lwd = style$lwd[m] * lwd_mult, cex = pt_cex
    )
  }
}

draw_legend_row <- function(labels, ids, style, ncol) {
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center",
    legend = labels,
    lty = mlty[ids], pch = mpch[ids],
    col = style$col[ids], pt.bg = style$bg[ids],
    lwd = style$lwd[ids] * lwd_mult, pt.cex = pt_cex,
    ncol = ncol, cex = 1.45, bty = "n"
  )
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

# 16:9 slide device: panel row on top, legend strip below.
open_slide_device <- function(out_file, n_panels) {
  pdf(out_file, width = 14, height = 4.9)
  layout(
    matrix(c(seq_len(n_panels), rep(n_panels + 1L, n_panels)),
      nrow = 2, byrow = TRUE
    ),
    heights = c(1, 0.24)
  )
  par(
    mar = c(4.4, 5.0, 2.8, 1.0),
    cex.axis = 1.45, cex.lab = 1.7, cex.main = 1.8
  )
}

# ---- Slide figure A: relative Brier vs quantile, 3 settings -----------
make_slide_rel_by_quantile <- function() {
  panels <- list(
    list(
      setting = 9L, ids = NULL,
      main = expression(paste(italic(K) == 1, ",   ", phi == "(0.80, -0.35)"))
    ),
    list(
      setting = 11L, ids = NULL,
      main = expression(paste(italic(K) == 5, ",   ", phi == "(0.80, -0.35)"))
    ),
    list(
      setting = 16L, ids = NULL,
      main = expression(paste(italic(K) == 1, ",   ", phi == "(0.15, 0.80)"))
    )
  )
  k0 <- 0L
  out_file <- file.path(out_dir, "slides-bs_rel_gauss-set9_11_16.pdf")
  open_slide_device(out_file, n_panels = length(panels))
  legend_ids <- integer(0)
  for (p in panels) {
    e <- load_sim(p$setting)
    stopifnot(k0 %in% e$mrts_ks)
    draw_ids <- if (is.null(p$ids)) e$methods[e$methods != 1L] else p$ids
    draw_cols <- match(draw_ids, e$methods)
    stopifnot(!anyNA(draw_cols))
    ki <- which(e$mrts_ks == k0)
    mat <- matrix(e$bs_rel_mean[, , ki],
      nrow = length(e$probs), ncol = length(e$methods)
    )
    ymat <- mat[, draw_cols, drop = FALSE]
    draw_series_slide(e$probs, ymat, draw_ids, style_ar2,
      ylim = c(min(ymat, 1, na.rm = TRUE), max(ymat, 1, na.rm = TRUE)),
      xlab = "Threshold quantile",
      ylab = "Relative Brier score",
      main = p$main
    )
    abline(h = 1, lty = 2, col = "gray60")
    legend_ids <- union(legend_ids, draw_ids)
  }
  legend_ids <- sort(legend_ids)
  draw_legend_row(vapply(legend_ids, label_for_id, character(1)),
    legend_ids, style_ar2,
    ncol = 4
  )
  dev.off()
  cat(sprintf("Wrote %s\n", normalizePath(out_file, winslash = "/", mustWork = FALSE)))
}

# ---- Slide figure B: Brier vs MRTS K at q=0.95, both nont DGPs --------
make_slide_mean_vs_K_nont <- function() {
  panels <- list(
    list(setting = 16L, main = "Arc front"),
    list(setting = 20L, main = "Transform mix")
  )
  q_show <- 0.95
  out_file <- file.path(out_dir, "slides-brier_vs_K-set16_20_nonsta-q95.pdf")
  open_slide_device(out_file, n_panels = length(panels))
  legend_ids <- integer(0)
  for (p in panels) {
    e <- load_sim(p$setting, "_nonsta")
    qi <- which(abs(e$probs - q_show) < 1e-9)
    stopifnot(length(qi) == 1L, length(e$mrts_ks) >= 2L)
    mat <- matrix(e$bs_mean[qi, , ],
      nrow = length(e$methods), ncol = length(e$mrts_ks)
    )
    draw_series_slide(e$mrts_ks, t(mat), e$methods, style_legacy,
      ylim = range(mat, na.rm = TRUE),
      xlab = "Number of MRTS basis functions",
      ylab = bquote("Brier score," ~ q == 0.95),
      main = p$main
    )
    legend_ids <- union(legend_ids, e$methods)
  }
  draw_legend_row(
    vapply(sort(legend_ids), label_for_id, character(1)),
    sort(legend_ids), style_legacy,
    ncol = 3
  )
  dev.off()
  cat(sprintf("Wrote %s\n", normalizePath(out_file, winslash = "/", mustWork = FALSE)))
}

make_slide_rel_by_quantile()
make_slide_mean_vs_K_nont()
