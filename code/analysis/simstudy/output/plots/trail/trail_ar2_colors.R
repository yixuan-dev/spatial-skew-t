#########################################################################
# trail_ar2_colors.R - DRAFT (trail) colour redesigns for the AR(2) figures.
#
# Goal: make the AR(2) methods (7, 8) stand out. New encoding principle:
#   colour = temporal class (AR(2) red, everything else NON-red),
#   lty     = knots (K=1 solid / K=5 dotted, unchanged from plots.R),
#   AR(2) drawn LAST so it sits on top of overlapping curves.
#
# Example figure: relative Brier vs threshold quantile, setting 11, K=0
# (the trail counterpart of
#  output/plots/brier_score/brier_score_rel_gauss_by_quantile-set11-K0.pdf).
#
# Variants (each in a bold [AR(2) lwd=2.5] and a _lwd1 [all lwd=1] version):
#   v1_gray         AR(2) red, ALL other methods greyscale
#   v2_coolcontrast AR(2) red, AR(1) blue, baselines grey, Sym-t slate
#   v3_legend_out   v2 styling + legend moved outside the plot region
#
# Usage (any cwd; the script locates itself):
#   Rscript.exe output/plots/trail/trail_ar2_colors.R
#
# DRAFT ONLY - plots.R / plots_for_paper.R are untouched until a variant
# is chosen.
#########################################################################

rm(list = ls())

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0L) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]),
    winslash = "/", mustWork = FALSE
  )
  if (dir.exists(dirname(script_path))) setwd(dirname(script_path))
}

# trail/ lives three levels below the simstudy dir: output/plots/trail
source("../../../helpers.R")
load("../../results/simresults11.RData")
# Provides (among others): bs_rel_mean, probs, methods, mrts_ks.

method_catalog <- get_simstudy_method_catalog(include_maxstable = TRUE)
method_label_for <- function(m) {
  row <- method_catalog[method_catalog$method_id == m, , drop = FALSE]
  if (nrow(row) == 1L) sprintf("%d: %s", m, row$label[1]) else as.character(m)
}
method_label <- vapply(methods, method_label_for, character(1))

n_probs <- length(probs)
n_methods <- length(methods)
stopifnot(all(methods >= 1L & methods <= 10L))

ki <- which(mrts_ks == 0L)
stopifnot(length(ki) == 1L)
mat <- matrix(bs_rel_mean[, , ki], nrow = n_probs, ncol = n_methods)

ar2_ids <- c(7L, 8L)
set_title <- bquote("Skew-" * italic(t) * " (K=5, " * lambda * "=3), AR(2): " *
  phi[1] * "=0.8, " * phi[2] * "=-0.35")

# ---- style tables (indexed positionally by method_id 1..10) -----------
# lty/pch keep the plots.R conventions (lty = K, pch = temporal class).
mlty <- c(1, 1, 1, 3, 3, 6, 1, 3, 1, 3)
mpch <- c(21, 22, 22, 22, 22, 23, 3, 3, 4, 4)

styles <- list(
  v1_gray = list(
    col = c(
      "gray60", "gray35", "gray55", "gray35", "gray55",
      "gray70", "firebrick3", "firebrick3", "gray45", "gray45"
    ),
    bg = c(
      "gray80", "gray55", "gray75", "gray55", "gray75",
      "gray85", "firebrick1", "firebrick1", "gray65", "gray65"
    )
  ),
  v2_coolcontrast = list(
    col = c(
      "gray55", "gray25", "slategray", "gray25", "slategray",
      "gray70", "firebrick3", "firebrick3", "dodgerblue3", "dodgerblue3"
    ),
    bg = c(
      "gray80", "gray50", "lightsteelblue1", "gray50", "lightsteelblue1",
      "gray85", "firebrick1", "firebrick1", "dodgerblue1", "dodgerblue1"
    )
  )
)

# ---- drawing ----------------------------------------------------------
# AR(2) columns are drawn last so the red curves sit on top.
draw_panel <- function(style, emphasize) {
  lwd_for <- function(m) if (emphasize && m %in% ar2_ids) 2.5 else 1
  draw_order <- c(
    which(!(methods %in% ar2_ids)),
    which(methods %in% ar2_ids)
  )
  plot(NA,
    xlim = range(probs),
    ylim = c(min(mat, 1, na.rm = TRUE), max(mat, 1, na.rm = TRUE)),
    xlab = "Threshold quantile",
    ylab = "Relative Brier score",
    main = set_title
  )
  abline(h = 1, lty = 2, col = "gray60")
  for (j in draw_order) {
    m <- methods[j]
    lines(probs, mat[, j], lty = mlty[m], col = style$col[m], lwd = lwd_for(m))
    points(probs, mat[, j],
      pch = mpch[m], col = style$col[m], bg = style$bg[m], lwd = lwd_for(m)
    )
  }
}

panel_legend <- function(style, emphasize, pos = "topleft", cex = 0.7) {
  lwds <- ifelse(emphasize & methods %in% ar2_ids, 2.5, 1)
  legend(pos,
    legend = method_label,
    lty = mlty[methods], pch = mpch[methods],
    col = style$col[methods], pt.bg = style$bg[methods], lwd = lwds,
    cex = cex, bty = "n"
  )
}

out_file <- function(tag) sprintf("brier_rel_set11-K0_trail_%s.pdf", tag)

# v1 / v2: legend inside (same placement as plots.R, for comparability)
for (variant in names(styles)) {
  style <- styles[[variant]]
  for (emphasize in c(TRUE, FALSE)) {
    tag <- paste0(sub("^v(\\d)_", "v\\1_", variant), if (emphasize) "" else "_lwd1")
    pdf(out_file(tag), width = 7, height = 5)
    draw_panel(style, emphasize)
    panel_legend(style, emphasize)
    dev.off()
    cat("Wrote", out_file(tag), "\n")
  }
}

# v3: v2 styling, legend outside the plot region (right-hand panel)
for (emphasize in c(TRUE, FALSE)) {
  tag <- paste0("v3_legend_out", if (emphasize) "" else "_lwd1")
  pdf(out_file(tag), width = 9, height = 5)
  layout(matrix(1:2, nrow = 1), widths = c(2.6, 1.15))
  par(mar = c(4.5, 4.5, 3, 0.6))
  draw_panel(styles$v2_coolcontrast, emphasize)
  par(mar = c(0, 0, 0, 0))
  plot.new()
  panel_legend(styles$v2_coolcontrast, emphasize, pos = "left", cex = 0.8)
  dev.off()
  cat("Wrote", out_file(tag), "\n")
}
