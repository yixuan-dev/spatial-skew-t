#########################################################################
# tables.R - Stage 2 of the simstudy post-fit pipeline.
#
# Reads output/results/scores<setting><suffix>.RData (produced by
# scores.R) and emits CSV tables + an aggregated .RData artifact under
# output/tables/ and output/results/.
#
# Mirrors code/analysis/simstudy_prop/tables-prop.R, with prop_k -> mrts_k.
# Supersedes the analyze-side bits of results.R / mrts_helpfulness.R.
#
# Plotting moved to Stage 3 (plots.R), which consumes the simresults
# artifact written here. This script produces tables only and has no
# graphics output.
#
# Usage:
#   Rscript tables.R --setting=<id>
#                    [--data=<path>]
#
# Examples:
#   Rscript tables.R --setting=4
#   Rscript tables.R --setting=1 --data=simdata_def.RData
#
# Outputs (suffix = "" for simdata.RData, "_def" for simdata_def.RData, ...):
#   output/tables/score_long<setting><suffix>.csv
#   output/tables/score_mean<setting><suffix>.csv
#   output/tables/score_rel_gauss<setting><suffix>.csv
#   output/tables/multivar_score<setting><suffix>.csv
#   output/tables/best_method_per_K<setting><suffix>.csv
#   output/tables/lambda_coverage<setting><suffix>.csv
#   output/results/simresults<setting><suffix>.RData
#########################################################################

rm(list = ls())

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0L) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]),
                               winslash = "/", mustWork = FALSE)
  if (dir.exists(dirname(script_path))) setwd(dirname(script_path))
}

source("./helpers.R")

# ---- CLI parsing -----------------------------------------------------
cli_args <- commandArgs(trailingOnly = TRUE)
parsed <- extract_leading_flags(cli_args, c("data", "setting"))
flags <- parsed$values

if (is.null(flags$setting) || !nzchar(flags$setting)) {
  stop("tables.R: --setting=<id> is required.", call. = FALSE)
}
setting_id <- as.integer(parse_index_expr(flags$setting, "setting"))
if (length(setting_id) != 1L || setting_id < 1L) {
  stop("tables.R: --setting must be a single positive integer.", call. = FALSE)
}

# ---- locate stage-1 cache --------------------------------------------
data_suffix <- if (!is.null(flags$data) && nzchar(flags$data)) {
  derive_data_suffix(flags$data)
} else {
  ""
}

results_dir <- "output/results"
tables_dir  <- "output/tables"
for (d in c(results_dir, tables_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

scores_file <- file.path(
  results_dir,
  sprintf("scores%d%s.RData", setting_id, data_suffix)
)
if (!file.exists(scores_file)) {
  stop(sprintf(
    "Score cache not found: %s\n  Run scores.R --setting=%d%s first.",
    scores_file, setting_id,
    if (nzchar(data_suffix)) sprintf(" --data=%s", flags$data) else ""
  ), call. = FALSE)
}
load(scores_file)
# Provides: quant.score, brier.score, energy.score, vario.score,
#           beta.0/1/2, tau.alpha, tau.beta, rho, nu, gamma, lambda,
#           elapsed_sec, probs, intervals, vs_p,
#           mrts_ks, datasets, methods, setting,
#           data_path, data_suffix, results_dir   (provenance)
# energy.score / vario.score / vs_p are absent in caches written by an
# older scores.R; the multivariate-score block below is guarded for that.

cat(sprintf(
  "tables: setting=%d cache=%s suffix='%s'\n",
  setting_id, scores_file,
  if (nzchar(data_suffix)) data_suffix else "<none>"
))

n_probs   <- length(probs)
n_methods <- length(methods)
n_ks      <- length(mrts_ks)

# ---- mean / median scores ------------------------------------------
mean_apply <- function(a) apply(a, c(1, 3, 4), mean,   na.rm = TRUE)
med_apply  <- function(a) apply(a, c(1, 3, 4), median, na.rm = TRUE)

bs_mean <- mean_apply(brier.score)   # [probs, method, mrts_k]
qs_mean <- mean_apply(quant.score)
bs_med  <- med_apply(brier.score)
qs_med  <- med_apply(quant.score)

dn_mean <- list(quantile = as.character(probs),
                method   = as.character(methods),
                mrts_k   = as.character(mrts_ks))
dimnames(bs_mean) <- dimnames(qs_mean) <- dn_mean
dimnames(bs_med)  <- dimnames(qs_med)  <- dn_mean

# ---- score_long: per-(score, method, mrts_k, dataset, quantile) ------
flatten_long <- function(arr, score) {
  grid <- expand.grid(
    quantile = probs,
    dataset  = datasets,
    method   = methods,
    mrts_k   = mrts_ks,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$score <- score
  grid$value <- as.vector(arr)
  grid[, c("score", "method", "mrts_k", "dataset", "quantile", "value")]
}
score_long_table <- rbind(
  flatten_long(brier.score, "brier"),
  flatten_long(quant.score, "quant")
)
write.csv(
  score_long_table,
  file.path(tables_dir, sprintf("score_long%d%s.csv", setting_id, data_suffix)),
  row.names = FALSE
)

# ---- score_mean: aggregated per (score, method, mrts_k, quantile) ----
score_long <- function(arr_mean, arr_med, score) {
  out <- expand.grid(quantile = probs, method = methods, mrts_k = mrts_ks,
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  out$score  <- score
  out$mean   <- as.vector(arr_mean)
  out$median <- as.vector(arr_med)
  out
}
score_table <- rbind(
  score_long(bs_mean, bs_med, "brier"),
  score_long(qs_mean, qs_med, "quant")
)
write.csv(
  score_table[, c("score", "method", "mrts_k", "quantile", "mean", "median")],
  file.path(tables_dir, sprintf("score_mean%d%s.csv", setting_id, data_suffix)),
  row.names = FALSE
)

# ---- score_rel_gauss: relative scores vs method 1 (Gaussian) ---------
rel_to_gauss <- function(arr) {
  if (!"1" %in% dimnames(arr)$method) {
    return(array(NA_real_, dim = dim(arr), dimnames = dimnames(arr)))
  }
  g <- arr[, "1", , drop = FALSE]
  sweep(arr, c(1, 3),
        array(g, dim = c(dim(arr)[1], 1, dim(arr)[3])),
        "/")
}
bs_rel_mean <- rel_to_gauss(bs_mean)
qs_rel_mean <- rel_to_gauss(qs_mean)
bs_rel_med  <- rel_to_gauss(bs_med)
qs_rel_med  <- rel_to_gauss(qs_med)

rel_long <- function(arr_mean, arr_med, score) {
  out <- expand.grid(quantile = probs, method = methods, mrts_k = mrts_ks,
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  out$score      <- score
  out$rel_mean   <- as.vector(arr_mean)
  out$rel_median <- as.vector(arr_med)
  out
}
rel_table <- rbind(
  rel_long(bs_rel_mean, bs_rel_med, "brier"),
  rel_long(qs_rel_mean, qs_rel_med, "quant")
)
write.csv(
  rel_table[, c("score", "method", "mrts_k", "quantile", "rel_mean", "rel_median")],
  file.path(tables_dir, sprintf("score_rel_gauss%d%s.csv", setting_id, data_suffix)),
  row.names = FALSE
)

# ---- multivar_score: energy + variogram, aggregated over datasets ----
# energy.score / vario.score are [dataset, method, mrts_k]. Average and
# median over datasets per (method, mrts_k), and form the relative ratio
# vs method 1 (Gaussian), matching the score_rel_gauss convention.
has_multivar  <- exists("energy.score") && exists("vario.score")
has_pred_rmse <- exists("pred.rmse")
has_recovery  <- exists("recovery.rmse")
mv_summary <- function(arr, score) {
  a_mean <- apply(arr, c(2, 3), mean,   na.rm = TRUE)   # method x mrts_k
  a_med  <- apply(arr, c(2, 3), median, na.rm = TRUE)
  rel_of <- function(m) {
    if (!"1" %in% rownames(m)) return(m * NA_real_)
    sweep(m, 2, m["1", ], "/")
  }
  out <- expand.grid(method = methods, mrts_k = mrts_ks,
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  out$score      <- score
  out$mean       <- as.vector(a_mean)
  out$median     <- as.vector(a_med)
  out$rel_mean   <- as.vector(rel_of(a_mean))
  out$rel_median <- as.vector(rel_of(a_med))
  out[, c("score", "method", "mrts_k", "mean", "median",
          "rel_mean", "rel_median")]
}
mv_parts <- list()
if (has_multivar) {
  mv_parts <- c(mv_parts, list(mv_summary(energy.score, "energy"),
                               mv_summary(vario.score,  "variogram")))
}
if (has_pred_rmse) {                                     # point-prediction RMSE
  mv_parts <- c(mv_parts, list(mv_summary(pred.rmse, "pred_rmse")))
}
if (has_recovery) {                                      # mean-surface recovery
  mv_parts <- c(mv_parts, list(mv_summary(recovery.rmse, "recovery_rmse")))
}
if (length(mv_parts) > 0) {
  multivar_table <- do.call(rbind, mv_parts)
} else {
  multivar_table <- data.frame(
    score = character(0), method = integer(0), mrts_k = integer(0),
    mean = numeric(0), median = numeric(0),
    rel_mean = numeric(0), rel_median = numeric(0)
  )
  cat("  (no energy/variogram/pred-rmse scores in cache; multivar table empty)\n")
}
write.csv(
  multivar_table,
  file.path(tables_dir, sprintf("multivar_score%d%s.csv", setting_id, data_suffix)),
  row.names = FALSE
)

# ---- best_method_per_K: lowest mean per quantile ---------------------
best_combo <- function(arr_mean, score) {
  out <- data.frame(score = score, quantile = probs,
                    best_method = NA_integer_,
                    best_mrts_k = NA_integer_,
                    best_value  = NA_real_)
  for (i in seq_len(n_probs)) {
    # Force 2-D shape so the n_ks == 1 (or n_methods == 1) case still works.
    mat <- matrix(arr_mean[i, , , drop = FALSE],
                  nrow = n_methods, ncol = n_ks)
    if (all(is.na(mat))) next
    idx <- which(mat == min(mat, na.rm = TRUE), arr.ind = TRUE)[1, ]
    out$best_method[i] <- methods[idx[1]]
    out$best_mrts_k[i] <- mrts_ks[idx[2]]
    out$best_value[i]  <- mat[idx[1], idx[2]]
  }
  out
}
best_table <- rbind(best_combo(bs_mean, "brier"),
                    best_combo(qs_mean, "quant"))
write.csv(
  best_table,
  file.path(tables_dir, sprintf("best_method_per_K%d%s.csv", setting_id, data_suffix)),
  row.names = FALSE
)

# ---- lambda_coverage: 95% interval coverage for skew methods ---------
ci_lo_idx <- which(abs(intervals - 0.025) < 1e-12)
ci_hi_idx <- which(abs(intervals - 0.975) < 1e-12)
lambda_true  <- 3
skew_methods <- intersect(c(2L, 4L), methods)
cov_rows <- list()
if (length(ci_lo_idx) == 1L && length(ci_hi_idx) == 1L) {
  for (m in skew_methods) {
    mi <- match(m, methods)
    for (ki in seq_along(mrts_ks)) {
      lo <- lambda[ci_lo_idx, , mi, ki]
      hi <- lambda[ci_hi_idx, , mi, ki]
      cover <- mean(lo <= lambda_true & hi >= lambda_true, na.rm = TRUE)
      width <- mean(hi - lo, na.rm = TRUE)
      cov_rows[[length(cov_rows) + 1L]] <- data.frame(
        method        = m,
        mrts_k        = mrts_ks[ki],
        coverage_95   = cover,
        mean_width_95 = width,
        n             = sum(!is.na(lo) & !is.na(hi))
      )
    }
  }
}
cov_table <- if (length(cov_rows) > 0L) {
  do.call(rbind, cov_rows)
} else {
  data.frame(
    method        = integer(0),
    mrts_k        = integer(0),
    coverage_95   = numeric(0),
    mean_width_95 = numeric(0),
    n             = integer(0)
  )
}
write.csv(
  cov_table,
  file.path(tables_dir, sprintf("lambda_coverage%d%s.csv", setting_id, data_suffix)),
  row.names = FALSE
)

# ---- save aggregated .RData artifact (consumed by plots.R) ----------
# lambda is kept raw for the per-dataset CI figure; energy.score /
# vario.score are added when present so plots.R can render the
# multivariate-score figures without re-reading the Stage-1 cache.
simresults_file <- file.path(
  results_dir,
  sprintf("simresults%d%s.RData", setting_id, data_suffix)
)
simresults_objs <- c(
  "bs_mean", "qs_mean", "bs_med", "qs_med",
  "bs_rel_mean", "qs_rel_mean", "bs_rel_med", "qs_rel_med",
  "score_long_table", "score_table", "rel_table", "multivar_table",
  "best_table", "cov_table", "lambda",
  "probs", "mrts_ks", "methods", "datasets", "intervals", "setting",
  "data_suffix"
)
if (has_multivar) {
  simresults_objs <- c(simresults_objs, "energy.score", "vario.score")
}
if (has_pred_rmse) {
  simresults_objs <- c(simresults_objs, "pred.rmse")
}
if (has_recovery) {
  simresults_objs <- c(simresults_objs, "recovery.rmse")
}
save(list = simresults_objs, file = simresults_file)

cat("\nWrote tables to ", tables_dir, "\n", sep = "")
cat("Saved analysis objects to ", simresults_file, "\n", sep = "")
cat("Next: Rscript plots.R --setting=", setting_id,
    if (nzchar(data_suffix)) sprintf(" --data=%s", flags$data) else "",
    " for the figures\n", sep = "")
