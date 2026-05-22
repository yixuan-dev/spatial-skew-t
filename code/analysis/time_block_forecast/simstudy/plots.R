#########################################################################
# plots.R - Stage 3 of the time-block forecast post-fit pipeline.
#
# Reads output/results/simresults<setting><suffix>.RData (from tables.R)
# and renders the lead-time curve figures of Section 4.3 /
# Definition 7 ("Lead-time score curve") of ar2_rethink.tex.
#
# Each figure plots S_bar(h) for every method on common axes:
#   - a shaded +/- 1 SE band, the cross-block dispersion of
#     Definition 9 (drawn as a band rather than whiskers so the AR(2)
#     and i.i.d. uncertainty regions are directly comparable);
#   - a dashed vertical reference at the memory horizon h*, where
#     Proposition (horizon) predicts the curves converge;
#   - the crossing lead -- the first h at which AR(2) stops beating the
#     i.i.d. baseline, the headline quantity of the redesign -- marked
#     with a vertical rule and label.
#
# Outputs (suffix "" for simdata.RData):
#   output/plots/leadcurve_crps-set<setting><suffix>.pdf
#   output/plots/leadcurve_brier-set<setting><suffix>.pdf
#
# Usage:
#   Rscript plots.R --setting=<id> [--data=<path>]
#########################################################################

rm(list = ls())

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0L) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]),
    winslash = "/", mustWork = FALSE)
  if (dir.exists(dirname(script_path))) setwd(dirname(script_path))
}

source("./time_block_helpers.R")

# ---- CLI parsing ------------------------------------------------------
cli_args <- commandArgs(trailingOnly = TRUE)
parsed <- extract_leading_flags(cli_args, c("data", "setting"))
flags <- parsed$values

if (is.null(flags$setting) || !nzchar(flags$setting)) {
  stop("plots.R: --setting=<id> is required.", call. = FALSE)
}
setting_id <- as.integer(parse_index_expr(flags$setting, "setting"))

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

# ---- load the Stage 2 artifact ---------------------------------------
simresults_file <- file.path(results_dir,
  sprintf("simresults%d%s.RData", setting_id, data_suffix))
if (!file.exists(simresults_file)) {
  stop(sprintf("Aggregated results not found: %s\n  Run tables.R --setting=%d first.",
    simresults_file, setting_id), call. = FALSE)
}
load(simresults_file)
# Provides: curve_table, rel_table, crossing_table, joint_table,
#           h_star, probs, methods, datasets, setting, block_H,
#           block_seams, data_suffix

catalog <- get_tbf_method_catalog()
method_label <- vapply(methods, function(m) {
  row <- catalog[catalog$method_id == m, , drop = FALSE]
  if (nrow(row) == 1L) sprintf("%d: %s", m, row$label[1]) else as.character(m)
}, character(1))
H <- block_H
set_tag <- sprintf("set%d%s", setting_id, data_suffix)

cat(sprintf("plots: setting=%d cache=%s -> %s\n",
  setting_id, simresults_file, plots_dir))

# ---- style ------------------------------------------------------------
mcol <- c("gray30", "firebrick3", "dodgerblue3")
mpch <- c(21, 22, 24)
method_col <- function(mi) mcol[(mi - 1) %% length(mcol) + 1]
method_pch <- function(mi) mpch[(mi - 1) %% length(mpch) + 1]

#########################################################################
# plot_curve() - one lead-time curve figure (Definition 7 (i)).
#
#   score_name - "crps" or "brier"; selects rows of curve_table.
#   ylab       - y-axis label.
#
# Layering: SE bands first (so the mean curves and markers sit on top),
# then the h* reference, the mean curves, and finally the crossing-lead
# rules. Returns the path of the PDF written.
#########################################################################
plot_curve <- function(score_name, ylab) {
  ct <- curve_table[curve_table$score == score_name, ]
  if (nrow(ct) == 0L) {
    cat("  no rows for score '", score_name, "', skipped\n", sep = "")
    return(invisible(NULL))
  }
  pdf_file <- file.path(plots_dir,
    sprintf("leadcurve_%s-%s.pdf", score_name, set_tag))
  pdf(pdf_file, width = 7, height = 5)
  on.exit(dev.off())

  # y-range must span the full +/- 1 SE band, not just the means.
  ylim <- range(c(ct$mean - ct$se, ct$mean + ct$se, ct$mean), na.rm = TRUE)
  plot(NA, xlim = c(1, H), ylim = ylim,
    xlab = "Lead time h", ylab = ylab,
    main = sprintf("Lead-time curve - setting %d%s%s", setting_id, data_suffix,
      if (is.finite(h_star)) sprintf("  (h* = %d)", h_star) else ""))

  # memory-horizon reference: where AR(2) and i.i.d. should converge.
  if (is.finite(h_star) && h_star >= 1 && h_star <= H) {
    abline(v = h_star, lty = 3, col = "gray50")
  }

  # ---- shaded +/- 1 SE bands (drawn first, behind the curves) --------
  for (mi in seq_along(methods)) {
    sub <- ct[ct$method == methods[mi], ]
    sub <- sub[order(sub$lead), ]
    band <- sub[is.finite(sub$se) & is.finite(sub$mean), ]
    if (nrow(band) >= 2L) {
      cc <- method_col(mi)
      polygon(c(band$lead, rev(band$lead)),
        c(band$mean + band$se, rev(band$mean - band$se)),
        col = adjustcolor(cc, alpha.f = 0.18), border = NA)
    }
  }

  # ---- mean curves + per-lead markers --------------------------------
  for (mi in seq_along(methods)) {
    sub <- ct[ct$method == methods[mi], ]
    sub <- sub[order(sub$lead), ]
    cc <- method_col(mi)
    lines(sub$lead, sub$mean, col = cc, lwd = 2)
    points(sub$lead, sub$mean, pch = method_pch(mi), col = cc, bg = cc)
  }

  # ---- crossing-lead markers -----------------------------------------
  # crossing_table holds, per non-baseline method, the first lead at
  # which the AR(2) / i.i.d. ratio reaches 1 (the headline quantity).
  cross <- crossing_table[crossing_table$score == score_name, ]
  ytext <- ylim[1] + 0.95 * diff(ylim)
  for (j in seq_len(nrow(cross))) {
    cl <- cross$crossing_lead[j]
    if (!is.finite(cl) || cl < 1 || cl > H) next
    mi <- match(cross$method[j], methods)
    cc <- if (is.na(mi)) "black" else method_col(mi)
    abline(v = cl, lty = 2, lwd = 1.5, col = cc)
    text(cl, ytext, sprintf(" crossing (m%d): h=%d", cross$method[j], cl),
      pos = 4, col = cc, cex = 0.7, xpd = NA)
  }

  legend("topleft", legend = method_label,
    col = vapply(seq_along(methods), method_col, character(1)),
    pch = vapply(seq_along(methods), method_pch, numeric(1)),
    lwd = 2, bty = "n", cex = 0.8)
  invisible(pdf_file)
}

# ---- render -----------------------------------------------------------
plot_curve("crps", "Mean CRPS")
plot_curve("brier", "Mean Brier score")

cat("Wrote plots to ", plots_dir, "\n", sep = "")
