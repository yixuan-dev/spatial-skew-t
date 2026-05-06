## ---------------------------------------------------------------------------
## mrts_helpfulness.R
##
## Stage 2: load the per-(setting, K) score bundles produced by
## mrts_scores.R, pair MRTS variants against the corresponding K=0
## baseline, and write the comparison CSVs.
##
## Files consumed:
##   scores{setting}_{K}mrts.RData  (K = 0 is baseline)
##
## Convention: delta = MRTS - baseline; lower score is better, so
##             delta < 0 means MRTS is better.
##
## Run from this directory (after mrts_scores.R):
##   Rscript mrts_helpfulness.R
##
## Outputs (comparison_mrts/helpfulness/):
##   mrts_helpfulness_paired.csv   per-dataset paired deltas
##   mrts_helpfulness_summary.csv  per (setting, family, K, prob) summary
## ---------------------------------------------------------------------------

rm(list = ls())

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]),
                               winslash = "/", mustWork = FALSE)
  script_dir <- dirname(script_path)
  if (dir.exists(script_dir)) setwd(script_dir)
}

source("./mrts_cov_helpers.R")

output_dir <- trimws(Sys.getenv("SIMSTUDY_MRTS_OUTPUT_DIR",
                                unset = "comparison_mrts/helpfulness"))
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

## (setting, K) pairs to compare -- exactly what was simulated.
experiments <- list(
  list(setting = 4L, k_values = c(10L, 15L, 20L)),
  list(setting = 5L, k_values = c(15L))
)

method_catalog <- get_simstudy_method_catalog(include_maxstable = FALSE)

load_scores <- function(setting, mrts_k) {
  f <- sprintf("scores%d_%dmrts.RData", setting, mrts_k)
  if (!file.exists(f)) {
    stop(sprintf("Missing score bundle '%s'. Run mrts_scores.R first.", f),
         call. = FALSE)
  }
  e <- new.env()
  load(f, envir = e)
  list(
    quant   = e$quant.score,
    brier   = e$brier.score,
    probs   = e$probs,
    methods = e$methods,
    datasets = e$datasets
  )
}

## ---------------------------------------------------------------------------
## per-dataset paired deltas
## ---------------------------------------------------------------------------
paired_rows <- list()
row_id <- 1L

for (exp in experiments) {
  setting_id <- exp$setting
  base <- load_scores(setting_id, 0L)
  probs    <- base$probs
  methods  <- base$methods
  datasets <- base$datasets

  for (k in exp$k_values) {
    mrts <- load_scores(setting_id, k)
    stopifnot(identical(mrts$probs, probs),
              identical(mrts$methods, methods),
              identical(mrts$datasets, datasets))

    for (m_idx in seq_along(methods)) {
      method_id <- methods[m_idx]
      family <- method_catalog$family[method_catalog$method_id == method_id]
      label  <- method_catalog$label [method_catalog$method_id == method_id]

      for (d_idx in seq_along(datasets)) {
        dataset_id <- datasets[d_idx]
        for (q_idx in seq_along(probs)) {
          base_b <- base$brier[q_idx, d_idx, m_idx]
          mrts_b <- mrts$brier[q_idx, d_idx, m_idx]
          base_q <- base$quant[q_idx, d_idx, m_idx]
          mrts_q <- mrts$quant[q_idx, d_idx, m_idx]

          paired_rows[[row_id]] <- data.frame(
            setting        = setting_id,
            dataset        = dataset_id,
            method_id      = method_id,
            family         = family,
            method_label   = label,
            mrts_k         = k,
            quantile       = probs[q_idx],
            baseline_brier = base_b,
            mrts_brier     = mrts_b,
            delta_brier    = mrts_b - base_b,
            baseline_quant = base_q,
            mrts_quant     = mrts_q,
            delta_quant    = mrts_q - base_q,
            stringsAsFactors = FALSE
          )
          row_id <- row_id + 1L
        }
      }
    }
  }
}
paired <- do.call(rbind, paired_rows)
rownames(paired) <- NULL
write.csv(paired, file.path(output_dir, "mrts_helpfulness_paired.csv"),
          row.names = FALSE)

## ---------------------------------------------------------------------------
## per (setting, family, K, prob) summary
## ---------------------------------------------------------------------------
summary_grid <- unique(paired[, c("setting", "family", "method_id",
                                  "method_label", "mrts_k", "quantile")])
summary_grid <- summary_grid[order(summary_grid$setting,
                                   summary_grid$method_id,
                                   summary_grid$mrts_k,
                                   summary_grid$quantile), ]

summary_rows <- vector("list", nrow(summary_grid))
for (i in seq_len(nrow(summary_grid))) {
  g <- summary_grid[i, ]
  rows <- paired[paired$setting   == g$setting   &
                 paired$method_id == g$method_id &
                 paired$mrts_k    == g$mrts_k    &
                 paired$quantile  == g$quantile, ]

  db <- rows$delta_brier
  dq <- rows$delta_quant
  n_b <- sum(is.finite(db))
  n_q <- sum(is.finite(dq))

  summary_rows[[i]] <- data.frame(
    setting            = g$setting,
    family             = g$family,
    method_id          = g$method_id,
    method_label       = g$method_label,
    mrts_k             = g$mrts_k,
    quantile           = g$quantile,
    n_pairs            = nrow(rows),
    n_valid_brier      = n_b,
    n_valid_quant      = n_q,
    mean_delta_brier   = if (n_b > 0) mean(db,   na.rm = TRUE) else NA_real_,
    median_delta_brier = if (n_b > 0) median(db, na.rm = TRUE) else NA_real_,
    win_rate_brier     = if (n_b > 0) mean(db < 0, na.rm = TRUE) else NA_real_,
    mean_delta_quant   = if (n_q > 0) mean(dq,   na.rm = TRUE) else NA_real_,
    median_delta_quant = if (n_q > 0) median(dq, na.rm = TRUE) else NA_real_,
    win_rate_quant     = if (n_q > 0) mean(dq < 0, na.rm = TRUE) else NA_real_,
    stringsAsFactors = FALSE
  )
}
summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(output_dir, "mrts_helpfulness_summary.csv"),
          row.names = FALSE)

cat("\nWritten:\n",
    " -", file.path(output_dir, "mrts_helpfulness_paired.csv"), "\n",
    " -", file.path(output_dir, "mrts_helpfulness_summary.csv"), "\n")
