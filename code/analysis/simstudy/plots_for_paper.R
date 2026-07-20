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
# Figure 3 - relative Brier score vs threshold quantile, persistence DGPs,
#   spatial hold-out arm. Two panels; only the AR(1)/AR(2) K=1 skew-t
#   fits are drawn (Gaussian is the dashed y = 1 reference).
#   (a) set16 near-unit-root AR(2) phi=(0.15,0.80)
#   (b) set19 ARFIMA(0,d,0) d=0.45
#
# Figure 4 - mean Brier score at q95 by forecast lead, time-block arm
#   (block 1: fit days 1-50, forecast 51-65). Reads the cache written by
#   time_block_forecast/block1_positive_control/cache_blk1_brier.R and
#   averages over the first 10 datasets that pass the guard under all
#   three methods (iid / AR(1) / AR(2)). Also writes the lead-window
#   aggregate table used for the thesis table tab:timeblock-brier.
#   (a) tbf set5 near-unit-root   (b) tbf set7 ARFIMA d=0.45
#
# Usage (from code/analysis/simstudy/, via PowerShell Rscript.exe):
#   Rscript.exe plots_for_paper.R
#
# Outputs (myLatex/pdf):
#   bs_rel_gauss_by_quantile-set9_to_12-K0.pdf
#   brier_score_mean_vs_K.pdf
#   bs_rel_gauss_by_quantile-set16_19-K0.pdf
#   bs_q95_by_lead-blk1-set5_7.pdf
# Output (time_block_forecast/block1_positive_control/output/tables):
#   blk1_thesis_table.csv
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
# lty/pch are shared; colours come as two palettes, as in plots.R:
#  - style_ar2 (colour = temporal class; AR(2) red + thick, AR(1) blue,
#    rest muted) for Figure 1, the AR(2) experiment, so AR(2) stands out;
#  - style_legacy (colour = family; Skew-t red / Sym-t blue) for Figure 2,
#    the MRTS experiment, where skew vs symmetric is the story.
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

# Bare method labels (no "id:" prefix) for the paper legend.
method_catalog <- get_simstudy_method_catalog(include_maxstable = TRUE)
label_for_id <- function(m) {
  row <- method_catalog[method_catalog$method_id == m, , drop = FALSE]
  if (nrow(row) == 1L) row$label[1] else as.character(m)
}

# draw_series(): plot the first series, overlay the rest. Columns of `ymat`
# align positionally with `method_ids`; each series is styled by its method id
# via `style` (one of style_legacy / style_ar2). AR(2) is drawn last so its
# curves sit on top of overlapping ones.
draw_series <- function(x, ymat, method_ids, style, ...) {
  ord <- order(method_ids %in% ar2_ids) # stable: keeps order, AR(2) last
  ymat <- ymat[, ord, drop = FALSE]
  method_ids <- method_ids[ord]
  m1 <- method_ids[1]
  plot(x, ymat[, 1],
    type = "o",
    pch = mpch[m1], lty = mlty[m1], col = style$col[m1], bg = style$bg[m1],
    lwd = style$lwd[m1], ...
  )
  for (j in seq_along(method_ids)[-1]) {
    m <- method_ids[j]
    lines(x, ymat[, j], lty = mlty[m], col = style$col[m], lwd = style$lwd[m])
    points(x, ymat[, j],
      pch = mpch[m], col = style$col[m], bg = style$bg[m], lwd = style$lwd[m]
    )
  }
}

# draw_legend_region(): region 5 blank, then the shared vertical legend in
# region 6. Call once, after the four panels, with the union of drawn ids.
draw_legend_region <- function(method_ids, style) {
  par(mar = c(0, 0, 0, 0))
  plot.new() # region 5 (row 3, col 1): intentionally blank
  plot.new() # region 6 (row 3, col 2): shared legend
  legend("center",
    legend = vapply(method_ids, label_for_id, character(1)),
    lty = mlty[method_ids], pch = mpch[method_ids],
    col = style$col[method_ids], pt.bg = style$bg[method_ids],
    lwd = style$lwd[method_ids],
    ncol = 1, cex = 1.25, bty = "n"
  )
}

open_panel_device <- function(out_file) {
  pdf(out_file, width = 9, height = 11)
  layout(matrix(1:6, nrow = 3, byrow = TRUE), heights = c(1, 1, 0.8))
  par(mar = c(4.6, 4.8, 2.4, 1.0), cex.axis = 1.1, cex.lab = 1.35, cex.main = 1.5)
}

# Two-panel variant of the shared convention: panels (a)-(b) in row 1,
# one shared legend spanning both columns in row 2 (drawn as a single row).
open_panel_device_2 <- function(out_file) {
  pdf(out_file, width = 9, height = 5.6)
  layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(1, 0.32))
  par(mar = c(4.6, 4.8, 2.4, 1.0), cex.axis = 1.1, cex.lab = 1.35, cex.main = 1.5)
}

draw_legend_span <- function(labels, ids, style) {
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center",
    legend = labels,
    lty = mlty[ids], pch = mpch[ids],
    col = style$col[ids], pt.bg = style$bg[ids], lwd = style$lwd[ids],
    ncol = length(ids), cex = 1.25, bty = "n"
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
    draw_series(probs, ymat, draw_ids, style_ar2,
      ylim = c(min(ymat, 1, na.rm = TRUE), max(ymat, 1, na.rm = TRUE)),
      xlab = "Threshold quantile",
      ylab = "Relative Brier score",
      main = tags[i]
    )
    abline(h = 1, lty = 2, col = "gray60")
  }
  draw_legend_region(draw_ids, style_ar2)
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
    draw_series(mrts_ks, t(mat), methods, style_legacy,
      ylim = range(mat, na.rm = TRUE),
      xlab = "Number of MRTS basis functions",
      ylab = bquote("Brier score," ~ q == .(sprintf("%.3f", p$q))),
      main = p$tag
    )
    legend_ids <- union(legend_ids, methods)
  }
  draw_legend_region(sort(legend_ids), style_legacy)
  dev.off()
  cat(sprintf("Wrote %s\n", normalizePath(out_file, winslash = "/", mustWork = FALSE)))
}

# ---- Figure 3: relative Brier vs quantile, persistence DGPs 16 & 19 ---
# Method sets differ between the two settings (set16 has the full catalog,
# set19 only 1/7/9), so columns are matched per panel instead of asserting
# identical method vectors as make_fig_rel_by_quantile() does.
make_fig_rel_by_quantile_persist <- function() {
  sets <- c(16L, 19L)
  tags <- c("(a)", "(b)")
  draw_ids <- c(7L, 9L) # skew-t K=1 AR(2), AR(1); Gaussian is the y=1 reference
  k0 <- 0L
  out_file <- file.path(out_dir, "bs_rel_gauss_by_quantile-set16_19-K0.pdf")
  open_panel_device_2(out_file)
  for (i in seq_along(sets)) {
    e <- load_sim(sets[i])
    stopifnot(k0 %in% e$mrts_ks)
    draw_cols <- match(draw_ids, e$methods)
    stopifnot(!anyNA(draw_cols))
    ki <- which(e$mrts_ks == k0)
    mat <- matrix(e$bs_rel_mean[, , ki],
      nrow = length(e$probs), ncol = length(e$methods)
    )
    ymat <- mat[, draw_cols, drop = FALSE]
    draw_series(e$probs, ymat, draw_ids, style_ar2,
      ylim = c(min(ymat, 1, na.rm = TRUE), max(ymat, 1, na.rm = TRUE)),
      xlab = "Threshold quantile",
      ylab = "Relative Brier score",
      main = tags[i]
    )
    abline(h = 1, lty = 2, col = "gray60")
  }
  draw_legend_span(vapply(draw_ids, label_for_id, character(1)),
    draw_ids, style_ar2
  )
  dev.off()
  cat(sprintf("Wrote %s\n", normalizePath(out_file, winslash = "/", mustWork = FALSE)))
}

# ---- Figure 4: Brier q95 by forecast lead, time-block block 1 ----------
# Cache layout (cache_blk1_brier.R): brier[lead, prob, dataset, method,
# setting] with methods (iid, AR2, AR1) and settings (5, 7); pass[dataset,
# method, setting]. Curves are styled by borrowing Morris method ids:
# iid -> 2 (gray), AR(2) -> 7 (firebrick, thick), AR(1) -> 9 (dodgerblue).
make_fig_brier_by_lead_blk1 <- function() {
  blk1_dir <- "../time_block_forecast/block1_positive_control/output/tables"
  cache_file <- file.path(blk1_dir, "blk1_brier_cache.RData")
  if (!file.exists(cache_file)) {
    stop(sprintf("Cache not found: %s (run cache_blk1_brier.R first)", cache_file),
      call. = FALSE
    )
  }
  e <- new.env()
  load(cache_file, envir = e)
  qi95 <- which(e$probs == 0.95)
  H <- dim(e$brier)[1]
  tags <- c("(a)", "(b)")
  # cache method order (iid, AR2, AR1) -> style ids and legend labels
  style_ids <- c(2L, 7L, 9L)
  labels <- c(
    "Skew-t, K=1 (i.i.d. in time)", "Skew-t, K=1, AR(2)",
    "Skew-t, K=1, AR(1)"
  )

  out_file <- file.path(out_dir, "bs_q95_by_lead-blk1-set5_7.pdf")
  open_panel_device_2(out_file)
  rows <- list()
  for (si in seq_along(e$settings)) {
    clean <- head(which(apply(e$pass[, , si], 1, all)), e$n_target)
    stopifnot(length(clean) == e$n_target)
    lb <- apply(e$brier[, qi95, clean, , si, drop = FALSE], c(1, 4), mean)
    draw_series(seq_len(H), lb, style_ids, style_ar2,
      ylim = range(lb, na.rm = TRUE),
      xlab = "Forecast lead (days)",
      ylab = bquote("Brier score," ~ q == 0.95),
      main = tags[si]
    )
    abline(v = 5.5, lty = 3, col = "gray60")
    for (hs in list(1:5, 1:15)) {
      wm <- colMeans(lb[hs, , drop = FALSE])
      rows[[length(rows) + 1L]] <- data.frame(
        setting = e$settings[si],
        leads = sprintf("1-%d", max(hs)),
        n = length(clean),
        iid = wm[1], ar1 = wm[3], ar2 = wm[2],
        ar1_minus_iid = wm[3] - wm[1],
        ar2_minus_iid = wm[2] - wm[1],
        ar2_minus_ar1 = wm[2] - wm[3]
      )
    }
  }
  # legend in reading order: iid, AR(1), AR(2)
  draw_legend_span(labels[c(1, 3, 2)], style_ids[c(1, 3, 2)], style_ar2)
  dev.off()
  cat(sprintf("Wrote %s\n", normalizePath(out_file, winslash = "/", mustWork = FALSE)))

  tab <- do.call(rbind, rows)
  tab_file <- file.path(blk1_dir, "blk1_thesis_table.csv")
  write.csv(tab, tab_file, row.names = FALSE)
  cat(sprintf("Wrote %s\n", normalizePath(tab_file, winslash = "/", mustWork = FALSE)))
  print(tab, digits = 4)
}

make_fig_rel_by_quantile()
make_fig_mean_vs_K()
make_fig_rel_by_quantile_persist()
make_fig_brier_by_lead_blk1()
