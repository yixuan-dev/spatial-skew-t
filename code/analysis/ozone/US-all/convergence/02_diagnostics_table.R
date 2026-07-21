# convergence/02_diagnostics_table.R
# ---------------------------------------------------------------------------
# Single-chain convergence diagnostics for the nine headline settings, from
# the caches written by 01_extract_chains.R. Per (setting, fold, parameter):
# rank-normalized split R-hat (the 5000 post-burn draws split into 4
# contiguous segments as pseudo-chains; Vehtari et al. 2021), coda ESS, and
# the coda Geweke z-score.
#
# The ar2 backend stores tau.alpha/tau.beta doubled relative to the legacy
# backend, so the tables carry a backend column and no cross-backend numeric
# comparison should be made on those rows.
#
# Output: output/us-all/tables/convergence_diagnostics.csv  (long format)
#         output/us-all/tables/convergence_summary.csv      (setting x fold)
# ---------------------------------------------------------------------------

rm(list = ls())
.this <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), winslash = "/", mustWork = FALSE)) else "."
})
setwd(file.path(.this, ".."))
source("convergence/00_conv_lib.R")

paths <- conv_paths()

long <- list()
for (i in seq_len(nrow(conv_settings))) {
  s <- conv_settings$setting[i]
  slim <- load_setting_slim(s, paths)
  for (d in 1:2) {
    diag_d <- diagnose_single_chain(slim$folds[[d]])
    diag_d <- cbind(
      setting = s,
      fold = d,
      backend = slim$backend_detected[d],
      diag_d,
      stringsAsFactors = FALSE
    )
    long[[length(long) + 1L]] <- diag_d
  }
  cat(sprintf("diagnosed setting %d (%d params x 2 folds)\n",
              s, ncol(slim$folds[[1]])))
}
long <- do.call(rbind, long)
long$flag <- classify_flag(long$rhat, long$ess, long$degenerate)

num_cols <- c("mean", "sd", "rhat_bulk", "rhat_tail", "rhat", "ess", "geweke_z")
long[num_cols] <- lapply(long[num_cols], function(x) round(x, 4))
write.csv(long, file.path(paths$tables, "convergence_diagnostics.csv"),
          row.names = FALSE)

by_fold <- split(long, interaction(long$setting, long$fold, drop = TRUE))
summary_df <- do.call(rbind, lapply(by_fold, function(g) {
  ok <- is.finite(g$rhat)
  oke <- is.finite(g$ess)
  data.frame(
    setting = g$setting[1],
    fold = g$fold[1],
    backend = g$backend[1],
    n_params = nrow(g),
    max_rhat = if (any(ok)) max(g$rhat[ok]) else NA_real_,
    max_rhat_param = if (any(ok)) g$param[ok][which.max(g$rhat[ok])] else NA_character_,
    min_ess = if (any(oke)) min(g$ess[oke]) else NA_real_,
    min_ess_param = if (any(oke)) g$param[oke][which.min(g$ess[oke])] else NA_character_,
    n_rhat_gt_1.01 = sum(g$rhat > 1.01, na.rm = TRUE),
    n_rhat_gt_1.05 = sum(g$rhat > 1.05, na.rm = TRUE),
    n_rhat_gt_1.10 = sum(g$rhat > 1.10, na.rm = TRUE),
    n_ess_lt_400 = sum(g$ess < 400, na.rm = TRUE),
    n_ess_lt_100 = sum(g$ess < 100, na.rm = TRUE),
    n_geweke_gt_2 = sum(abs(g$geweke_z) > 2, na.rm = TRUE),
    n_geweke_gt_3 = sum(abs(g$geweke_z) > 3, na.rm = TRUE),
    n_degenerate = sum(g$degenerate),
    stringsAsFactors = FALSE
  )
}))
summary_df <- summary_df[order(summary_df$setting, summary_df$fold), ]
summary_df$max_rhat <- round(summary_df$max_rhat, 4)
summary_df$min_ess <- round(summary_df$min_ess, 1)
write.csv(summary_df, file.path(paths$tables, "convergence_summary.csv"),
          row.names = FALSE)

cat("\n== convergence_summary.csv ==\n")
print(summary_df, row.names = FALSE)
cat(sprintf("\nwrote %s and %s\n",
            file.path(paths$tables, "convergence_diagnostics.csv"),
            file.path(paths$tables, "convergence_summary.csv")))
