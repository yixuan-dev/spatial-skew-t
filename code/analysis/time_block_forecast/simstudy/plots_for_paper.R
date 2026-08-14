#########################################################################
# plots_for_paper.R - paper-ready figures for the TIME-BLOCK experiment.
#
# Companion to code/analysis/simstudy/plots_for_paper.R, which owns the
# spatial hold-out figures. The time-block arm keeps its own script so
# that a figure and the campaign it summarises live in the same
# directory: this file reads only this study's caches.
#
# Figure - paired Brier difference by forecast lead, Experiment A.
#   Reads output/results/scores{4,5,7}_hn.RData (lambda ~ HN(0, 20),
#   full-series thresholds, five expanding-window blocks, ten data sets,
#   methods 1 = i.i.d. / 4 = AR(1) / 2 = AR(2)). The five blocks are
#   averaged within a data set BEFORE the difference is formed, so the
#   quantity drawn is the quantity the paired tests use (n = 10).
#   Point estimates only: as in the spatial figures, the intervals live
#   in the table, not in the figure (thesis sec:sim-uq).
#   Panels are ordered by the persistence of the generating process:
#   (a) set4 (0.80, -0.35)   (b) set5 (0.15, 0.80)   (c) set7 ARFIMA d=0.45
#
# Usage (from code/analysis/time_block_forecast/simstudy/):
#   Rscript.exe plots_for_paper.R
#
# Outputs (myLatex/pdf) - the thesis figure:
#   bs_q95_leaddiff-set4_5_7.pdf
# Outputs (output/plots) - the study's own copy of the same figure:
#   leaddiff_brier_q95-set4_5_7.pdf
#
# Supersedes the block-1 figure of the retired campaign, whose frozen
# record is tex/timeblock_expABC_legacy.
#########################################################################

rm(list = ls())

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0L) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1]),
    winslash = "/", mustWork = FALSE
  ))
  if (dir.exists(script_dir)) setwd(script_dir)
}

results_dir <- "output/results"
plots_dir <- "output/plots"
paper_dir <- "../../../../myLatex/pdf"
for (d in c(plots_dir, paper_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Method labels by id, in the order scores.R stores them.
MLAB <- c("1" = "iid", "2" = "AR2", "4" = "AR1")

SETS <- list(
  list(id = 4L, tag = "(a)", lab = "AR(2) (0.80, -0.35)"),
  list(id = 5L, tag = "(b)", lab = "AR(2) (0.15, 0.80)"),
  list(id = 7L, tag = "(c)", lab = "ARFIMA(0, d, 0), d = 0.45")
)

# Colours match style_ar2 of code/analysis/simstudy/plots_for_paper.R
# (AR(2) firebrick, AR(1) dodgerblue), so the two arms of the thesis read
# as one visual system. AR(2) - AR(1) is drawn first and solid: it is the
# second-lag contrast the section turns on.
CONTRASTS <- list(
  list(a = "AR2", b = "AR1", lab = "AR(2) - AR(1)",
       col = "firebrick3", lty = 1, pch = 16, lwd = 2.6),
  list(a = "AR2", b = "iid", lab = "AR(2) - i.i.d.",
       col = "dodgerblue3", lty = 2, pch = 17, lwd = 2.2)
)

# Three panels in one row with the shared legend spanning the bottom. The
# canvas is 9 in wide like the spatial paper figures, so that at
# \includegraphics[width=\textwidth] every figure in the thesis renders
# its type at the same size.
open_panel_device_3 <- function(out_file) {
  pdf(out_file, width = 9, height = 4.4)
  layout(matrix(c(1, 2, 3, 4, 4, 4), nrow = 2, byrow = TRUE), heights = c(1, 0.24))
  par(mar = c(4.4, 4.6, 2.4, 0.8), cex.axis = 1.0, cex.lab = 1.25, cex.main = 1.3)
}

# lead_differences(): per-lead paired difference for one setting, in units
# of 1e-3. Gates the cache on the two provenance fields scores.R writes,
# so a figure can never be drawn from the wrong threshold rule or the
# wrong prior arm.
lead_differences <- function(setting_id, prob) {
  f <- file.path(results_dir, sprintf("scores%d_hn.RData", setting_id))
  if (!file.exists(f)) {
    stop(sprintf("Cache not found: %s (run the Experiment A driver first)", f),
      call. = FALSE
    )
  }
  e <- new.env()
  load(f, envir = e)
  stopifnot(identical(e$brier_threshold_basis, "full_series"), isTRUE(e$hn_prior))
  qi <- which.min(abs(e$probs - prob))
  # brier.lead[lead, prob, dataset, method, block] -> lead x dataset x method
  a <- apply(e$brier.lead[, qi, , , , drop = FALSE], c(1, 3, 4), mean)
  dimnames(a)[[3]] <- MLAB[as.character(e$methods)]
  vapply(CONTRASTS, function(ct) rowMeans(a[, , ct$a] - a[, , ct$b]) * 1000,
    numeric(dim(a)[1])
  )
}

make_fig_leaddiff <- function(prob, out_file) {
  per_setting <- lapply(SETS, function(s) lead_differences(s$id, prob))
  H <- nrow(per_setting[[1]])
  ylim <- range(unlist(per_setting))

  open_panel_device_3(out_file)
  for (i in seq_along(SETS)) {
    d <- per_setting[[i]]
    plot(NA,
      xlim = c(1, H), ylim = ylim,
      xlab = "Forecast lead (days)",
      # the units go in parentheses rather than behind a multiplication
      # sign: "%*%" is set from the Symbol font, which not every PDF
      # rasteriser carries
      ylab = if (i == 1) expression("Brier difference (" * 10^-3 * ")") else "",
      # panel tag only: the caption carries the data-generating design,
      # so the panel head stays as light as in the other paper figures
      main = SETS[[i]]$tag
    )
    abline(h = 0, col = "gray35", lwd = 1.4)
    for (j in seq_along(CONTRASTS)) {
      ct <- CONTRASTS[[j]]
      lines(seq_len(H), d[, j], col = ct$col, lty = ct$lty, lwd = ct$lwd)
      points(seq_len(H), d[, j], col = ct$col, pch = ct$pch, cex = 0.9)
    }
  }
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center",
    legend = vapply(CONTRASTS, function(x) x$lab, character(1)),
    col = vapply(CONTRASTS, function(x) x$col, character(1)),
    lty = vapply(CONTRASTS, function(x) x$lty, numeric(1)),
    pch = vapply(CONTRASTS, function(x) x$pch, numeric(1)),
    lwd = vapply(CONTRASTS, function(x) x$lwd, numeric(1)),
    # text.width forces the legend columns apart; with the default the two
    # entries pack against each other and read as one string
    ncol = length(CONTRASTS), cex = 1.25, bty = "n",
    text.width = 0.17, x.intersp = 0.9
  )
  dev.off()
  cat(sprintf("Wrote %s\n", normalizePath(out_file, winslash = "/", mustWork = FALSE)))
}

# the thesis figure
make_fig_leaddiff(0.95, file.path(paper_dir, "bs_q95_leaddiff-set4_5_7.pdf"))
# the study's own copy, for the record alongside the campaign it summarises
make_fig_leaddiff(0.95, file.path(plots_dir, "leaddiff_brier_q95-set4_5_7.pdf"))
