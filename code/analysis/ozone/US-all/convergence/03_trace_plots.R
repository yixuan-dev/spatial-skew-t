# convergence/03_trace_plots.R
# ---------------------------------------------------------------------------
# Trace plots of the saved post-burn-in draws (iterations 25001-30000,
# thin=1) for every diagnosed scalar chain, one multi-panel PNG per
# setting x fold. Panel titles carry the split R-hat and ESS from
# 02_diagnostics_table.R; marginal panels are titled orange, failing ones
# red (thresholds in the README).
#
# Requires: caches from 01 and output/us-all/tables/convergence_diagnostics.csv.
# Output:   output/us-all/plots/convergence/trace-us-all-<N>-fold<d>.png
# ---------------------------------------------------------------------------

rm(list = ls())
.this <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), winslash = "/", mustWork = FALSE)) else "."
})
setwd(file.path(.this, ".."))
source("convergence/00_conv_lib.R")

paths <- conv_paths()

diag_csv <- file.path(paths$tables, "convergence_diagnostics.csv")
if (!file.exists(diag_csv)) stop("run 02_diagnostics_table.R first", call. = FALSE)
diags <- read.csv(diag_csv, stringsAsFactors = FALSE)

title_col <- function(flag) {
  switch(flag, fail = "red", degenerate = "red", marginal = "darkorange", "black")
}

for (i in seq_len(nrow(conv_settings))) {
  s <- conv_settings$setting[i]
  slim <- load_setting_slim(s, paths)
  for (d in 1:2) {
    chains <- slim$folds[[d]]
    P <- ncol(chains)
    nd <- nrow(chains)
    iter_idx <- seq_len(nd) + if (nd == 5000L) 25000L else 0L
    dg <- diags[diags$setting == s & diags$fold == d, ]

    nc <- ceiling(sqrt(P))
    nr <- ceiling(P / nc)
    out <- file.path(paths$plots, sprintf("trace-us-all-%d-fold%d.png", s, d))
    png(out, width = 480 * nc, height = 340 * nr, res = 110)
    par(mfrow = c(nr, nc), mar = c(3.2, 3.2, 2.6, 0.8), mgp = c(2.1, 0.7, 0),
        oma = c(0, 0, 2.2, 0))
    for (pn in colnames(chains)) {
      row <- dg[dg$param == pn, ]
      subtitle <- if (nrow(row) == 1 && is.finite(row$rhat)) {
        sprintf("%s  Rhat=%.3f  ESS=%.0f", pn, row$rhat, row$ess)
      } else if (nrow(row) == 1 && row$degenerate) {
        sprintf("%s  [degenerate]", pn)
      } else {
        pn
      }
      cl <- if (nrow(row) == 1) title_col(row$flag) else "black"
      plot(iter_idx, chains[, pn], type = "l", lwd = 0.4,
           xlab = "iteration", ylab = pn, main = "")
      title(main = subtitle, col.main = cl, cex.main = 1.0)
    }
    mtext(sprintf("us-all setting %d, fold %d (%s backend, %d post-burn draws)",
                  s, d, slim$backend_detected[d], nd),
          side = 3, line = 0.6, outer = TRUE, cex = 1.0)
    dev.off()
    cat(sprintf("wrote %s (%d panels)\n", out, P))
  }
}
