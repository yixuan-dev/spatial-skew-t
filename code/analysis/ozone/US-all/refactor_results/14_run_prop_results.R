rm(list = ls())

source("refactor_results/00_bootstrap.R")
source("refactor_results/01_score_engine.R")
source("refactor_results/02_comparison_tables.R")

set_working_dir_to_script()
ctx <- load_us_all_context()
ensure_us_all_output_dirs()

prop_settings_path <- normalizePath(
  file.path("..", "..", "ozone_prop", "US-all", "settings_prop.csv"),
  winslash = "/", mustWork = TRUE
)
prop_results_dir <- normalizePath(
  file.path("..", "..", "ozone_prop", "US-all", "results"),
  winslash = "/", mustWork = TRUE
)
gaussian_result_path <- normalizePath(
  file.path("results", "us-all-1.RData"),
  winslash = "/", mustWork = TRUE
)

prop_settings_raw <- read.csv(prop_settings_path, stringsAsFactors = FALSE)
prop_settings_raw$setting_num_orig <- safe_as_integer(prop_settings_raw$setting)
if (any(is.na(prop_settings_raw$setting_num_orig))) {
  stop("settings_prop.csv has non-numeric setting IDs.")
}

prop_settings_raw <- prop_settings_raw[order(prop_settings_raw$setting_num_orig), , drop = FALSE]

# Map prop setting k -> unified ID k+1, reserve ID 1 for the Gaussian baseline.
prop_offset <- 1L
prop_settings_raw$setting_num_unified <- prop_settings_raw$setting_num_orig + prop_offset

baseline_row <- ctx$settings[!is.na(ctx$settings$setting_num) & ctx$settings$setting_num == 1L, , drop = FALSE]
if (nrow(baseline_row) == 0) {
  stop("No setting 1 (Gaussian) found in main settings.csv")
}

unified_settings <- data.frame(
  setting     = c(1L, prop_settings_raw$setting_num_unified),
  setting_num = c(1L, prop_settings_raw$setting_num_unified),
  method      = c(meta_value(baseline_row, "method", "gaussian"),
                  prop_settings_raw$method),
  knots       = c(meta_value(baseline_row, "knots", 1L),
                  prop_settings_raw$knots),
  thresh      = c(meta_value(baseline_row, "thresh", 0),
                  prop_settings_raw$thresh),
  CMAQ        = c(meta_value(baseline_row, "CMAQ", "yes"),
                  prop_settings_raw$CMAQ),
  TS          = c(meta_value(baseline_row, "TS", NA),
                  rep(NA_character_, nrow(prop_settings_raw))),
  ar2         = c(meta_value(baseline_row, "ar2", NA),
                  rep(NA_character_, nrow(prop_settings_raw))),
  mrts        = c(NA_character_, as.character(prop_settings_raw$mrts)),
  setting_orig = c(1L, prop_settings_raw$setting_num_orig),
  source      = c("ozone_us_all_results",
                  rep("ozone_prop_us_all_results", nrow(prop_settings_raw))),
  stringsAsFactors = FALSE
)

prop_setting_ids <- prop_settings_raw$setting_num_unified
all_requested <- sort(unique(c(1L, prop_setting_ids)))
baseline_ids <- 1L
proposed_ids <- prop_setting_ids

result_file_map <- character(length(all_requested))
names(result_file_map) <- as.character(all_requested)
result_file_map[["1"]] <- gaussian_result_path
for (k in seq_along(prop_settings_raw$setting_num_orig)) {
  unified_id <- prop_settings_raw$setting_num_unified[k]
  orig_id    <- prop_settings_raw$setting_num_orig[k]
  result_file_map[[as.character(unified_id)]] <- file.path(
    prop_results_dir, sprintf("ozone-prop-%d.RData", orig_id)
  )
}

decorate_setting_table <- function(df, setting_col = "setting") {
  if (nrow(df) == 0) {
    return(df)
  }

  idx <- match(df[[setting_col]], unified_settings$setting_num)
  df$mrts          <- unified_settings$mrts[idx]
  df$setting_orig  <- unified_settings$setting_orig[idx]
  df$source        <- unified_settings$source[idx]
  df$is_baseline   <- df[[setting_col]] %in% baseline_ids
  df$is_proposed   <- df[[setting_col]] %in% proposed_ids
  df$model_lane    <- ifelse(df[[setting_col]] == 1L, "gaussian_reference", "proposed")
  df
}

probs <- default_probs(include_999 = FALSE)
threshold_probs <- default_threshold_probs()
thresholds <- quantile(ctx$Y, probs = threshold_probs, na.rm = TRUE)

summary_draws <- suppressWarnings(as.integer(Sys.getenv("US_ALL_SUMMARY_DRAWS", unset = "400")))
if (!is.finite(summary_draws) || summary_draws <= 0) {
  summary_draws <- 400L
}

score_obj <- compute_us_all_scores(
  setting_ids = all_requested,
  result_path_fn = function(i) result_file_map[[as.character(i)]],
  Y = ctx$Y,
  cv_lst = ctx$cv_lst,
  probs = probs,
  thresholds = thresholds,
  threshold_probs = threshold_probs,
  compute_brier_split_diagnostics = TRUE,
  compute_classification_diagnostics = TRUE,
  trans_setting_ids = integer(0),
  enforce_contract = TRUE,
  compute_uncertainty_diagnostics = TRUE,
  summary_draws = summary_draws
)
summary_obj <- summarize_us_all_scores(score_obj, baseline_setting = 1L)
print_score_summary(score_obj, label = "us-all-results-proposed (ozone_prop fits 1-10 + Gaussian baseline)")

comparison_full_table <- build_comparison_full_table(
  summary_obj = summary_obj,
  settings = unified_settings,
  baseline_ids = baseline_ids,
  proposed_ids = proposed_ids
)
if (nrow(comparison_full_table) > 0) {
  comparison_full_table <- decorate_setting_table(comparison_full_table, "setting")
}

comparison_top2 <- build_comparison_top2_all_metrics(
  summary_obj = summary_obj,
  settings = unified_settings,
  candidate_settings = summary_obj$available_settings
)
if (nrow(comparison_top2) > 0) {
  comparison_top2 <- decorate_setting_table(comparison_top2, "setting")
}

comparison_paired_same_basis <- build_paired_same_basis_table(
  summary_obj = summary_obj,
  settings = unified_settings,
  baseline_ids = baseline_ids,
  proposed_ids = proposed_ids
)

comparison_scalar_metrics <- build_comparison_scalar_metrics_table(
  summary_obj = summary_obj,
  settings = unified_settings,
  baseline_ids = baseline_ids,
  proposed_ids = proposed_ids
)
if (nrow(comparison_scalar_metrics) > 0) {
  comparison_scalar_metrics <- decorate_setting_table(comparison_scalar_metrics, "setting")
}

comparison_classification_metrics <- build_comparison_classification_metrics_table(
  summary_obj = summary_obj,
  settings = unified_settings,
  baseline_ids = baseline_ids,
  proposed_ids = proposed_ids
)
if (nrow(comparison_classification_metrics) > 0) {
  comparison_classification_metrics <- decorate_setting_table(comparison_classification_metrics, "setting")
}

comparison_brier_split <- build_comparison_brier_split_table(
  summary_obj = summary_obj,
  settings = unified_settings,
  baseline_ids = baseline_ids,
  proposed_ids = proposed_ids
)
if (nrow(comparison_brier_split) > 0) {
  comparison_brier_split <- decorate_setting_table(comparison_brier_split, "setting")
}

comparison_uncertainty_summary <- build_comparison_uncertainty_summary_table(
  summary_obj = summary_obj,
  settings = unified_settings,
  baseline_ids = baseline_ids,
  proposed_ids = proposed_ids
)
if (nrow(comparison_uncertainty_summary) > 0) {
  comparison_uncertainty_summary <- decorate_setting_table(comparison_uncertainty_summary, "setting")
}

comparison_calibration_bins <- build_comparison_calibration_bins_table(
  summary_obj = summary_obj,
  settings = unified_settings,
  baseline_ids = baseline_ids,
  proposed_ids = proposed_ids
)
if (nrow(comparison_calibration_bins) > 0) {
  comparison_calibration_bins <- decorate_setting_table(comparison_calibration_bins, "setting")
}

write.csv(comparison_full_table,
          us_all_output_path("comparison_full_table_proposed.csv", subdir = "tables"),
          row.names = FALSE)
write.csv(comparison_top2,
          us_all_output_path("comparison_top2_proposed.csv", subdir = "tables"),
          row.names = FALSE)
write_comparison_top2_workbook(
  comparison_top2 = comparison_top2,
  output_path = us_all_output_path("comparison_top2_proposed.xlsx", subdir = "tables")
)
write.csv(comparison_paired_same_basis,
          us_all_output_path("comparison_paired_same_basis_proposed.csv", subdir = "tables"),
          row.names = FALSE)
write.csv(comparison_scalar_metrics,
          us_all_output_path("comparison_scalar_metrics_proposed.csv", subdir = "tables"),
          row.names = FALSE)
write.csv(comparison_classification_metrics,
          us_all_output_path("comparison_classification_metrics_proposed.csv", subdir = "tables"),
          row.names = FALSE)
write.csv(comparison_brier_split,
          us_all_output_path("comparison_brier_split_proposed.csv", subdir = "tables"),
          row.names = FALSE)
write.csv(comparison_uncertainty_summary,
          us_all_output_path("comparison_uncertainty_summary_proposed.csv", subdir = "tables"),
          row.names = FALSE)
write.csv(comparison_calibration_bins,
          us_all_output_path("comparison_calibration_bins_proposed.csv", subdir = "tables"),
          row.names = FALSE)

available_settings <- summary_obj$available_settings
skipped_missing_file <- summary_obj$skipped_missing_file
skipped_bad_contract <- summary_obj$skipped_bad_contract
skipped_scoring_error <- summary_obj$skipped_scoring_error
quant.score <- summary_obj$quant.score
brier.score <- summary_obj$brier.score
mspe.score <- summary_obj$mspe.score
mape.score <- summary_obj$mape.score
quant.score.mean <- summary_obj$quant.score.mean
brier.score.mean <- summary_obj$brier.score.mean
mspe.mean <- summary_obj$mspe.mean
mape.mean <- summary_obj$mape.mean
quant.score.se <- summary_obj$quant.score.se
brier.score.se <- summary_obj$brier.score.se
mspe.se <- summary_obj$mspe.se
mape.se <- summary_obj$mape.se
bs.mean.ref.gau <- summary_obj$bs.mean.ref.gau
qs.mean.ref.gau <- summary_obj$qs.mean.ref.gau
mspe.mean.ref.gau <- summary_obj$mspe.mean.ref.gau
mape.mean.ref.gau <- summary_obj$mape.mean.ref.gau
brier.split.target_probs <- summary_obj$brier.split.target_probs
brier.split.target_thresholds <- summary_obj$brier.split.target_thresholds
brier.split.band_names <- summary_obj$brier.split.band_names
brier.split.score.mean <- summary_obj$brier.split.score.mean
brier.split.score.se <- summary_obj$brier.split.score.se
brier.split.n_obs.total <- summary_obj$brier.split.n_obs.total
brier.split.n_obs.mean <- summary_obj$brier.split.n_obs.mean
brier.split.obs.share <- summary_obj$brier.split.obs.share
brier.split.rel.ref.gau <- summary_obj$brier.split.rel.ref.gau
classification.target_probs <- summary_obj$classification.target_probs
classification.target_thresholds <- summary_obj$classification.target_thresholds
classification.metric.names <- summary_obj$classification.metric_names
classification.count.names <- summary_obj$classification.count_names
classification.metric.mean <- summary_obj$classification.metric.mean
classification.metric.se <- summary_obj$classification.metric.se
classification.count.total <- summary_obj$classification.count.total
classification.count.mean <- summary_obj$classification.count.mean
classification.n_obs.total <- summary_obj$classification.n_obs.total
classification.n_obs.mean <- summary_obj$classification.n_obs.mean
classification.actual_positive_share <- summary_obj$classification.actual_positive_share
classification.predicted_positive_share <- summary_obj$classification.predicted_positive_share
classification.metric.rel.ref.gau <- summary_obj$classification.metric.rel.ref.gau
classification.metric.delta.ref.gau <- summary_obj$classification.metric.delta.ref.gau
classification.probability_cutoff <- summary_obj$classification.probability_cutoff
crps.mean <- summary_obj$crps.mean
crps.se <- summary_obj$crps.se
crps.mean.ref.gau <- summary_obj$crps.mean.ref.gau
coverage.mean <- summary_obj$coverage.mean
coverage.se <- summary_obj$coverage.se
coverage.gap <- summary_obj$coverage.gap
pit.mean <- summary_obj$pit.mean
pit.variance <- summary_obj$pit.variance
pit.ks <- summary_obj$pit.ks
pit.mae <- summary_obj$pit.mae
pit.rmse <- summary_obj$pit.rmse
pit.bin.share.mean <- summary_obj$pit.bin.share.mean
pit.bin.share.se <- summary_obj$pit.bin.share.se
summary.draws.mean <- summary_obj$summary.draws.mean
summary.n_obs.total <- summary_obj$summary.n_obs.total

settings <- unified_settings
result_file_map_proposed <- result_file_map

save(
  list = c(
    "settings", "baseline_ids", "proposed_ids",
    "all_requested", "result_file_map_proposed", "available_settings",
    "skipped_missing_file", "skipped_bad_contract", "skipped_scoring_error",
    "probs", "threshold_probs", "thresholds", "quant.score", "brier.score",
    "mspe.score", "mape.score",
    "quant.score.mean", "brier.score.mean", "mspe.mean", "mape.mean",
    "quant.score.se", "brier.score.se", "mspe.se", "mape.se",
    "bs.mean.ref.gau", "qs.mean.ref.gau", "mspe.mean.ref.gau", "mape.mean.ref.gau",
    "brier.split.target_probs", "brier.split.target_thresholds", "brier.split.band_names",
    "brier.split.score.mean", "brier.split.score.se", "brier.split.n_obs.total",
    "brier.split.n_obs.mean", "brier.split.obs.share", "brier.split.rel.ref.gau",
    "classification.target_probs", "classification.target_thresholds",
    "classification.metric.names", "classification.count.names",
    "classification.metric.mean", "classification.metric.se",
    "classification.count.total", "classification.count.mean",
    "classification.n_obs.total", "classification.n_obs.mean",
    "classification.actual_positive_share", "classification.predicted_positive_share",
    "classification.metric.rel.ref.gau", "classification.metric.delta.ref.gau",
    "classification.probability_cutoff",
    "crps.mean", "crps.se", "crps.mean.ref.gau",
    "coverage.mean", "coverage.se", "coverage.gap",
    "pit.mean", "pit.variance", "pit.ks", "pit.mae", "pit.rmse",
    "pit.bin.share.mean", "pit.bin.share.se",
    "summary.draws.mean", "summary.n_obs.total",
    "comparison_full_table", "comparison_top2", "comparison_paired_same_basis",
    "comparison_scalar_metrics", "comparison_classification_metrics",
    "comparison_brier_split", "comparison_uncertainty_summary", "comparison_calibration_bins",
    "score_obj", "summary_obj"
  ),
  file = us_all_output_path("us-all-results-proposed.RData", subdir = "results")
)

cat("Result files used:\n")
for (id in all_requested) {
  cat("  setting", id, "->", result_file_map[[as.character(id)]], "\n")
}

cat("Outputs written:\n")
cat("- ", us_all_output_path("us-all-results-proposed.RData", subdir = "results"), "\n", sep = "")
cat("- ", us_all_output_path("comparison_full_table_proposed.csv", subdir = "tables"), "\n", sep = "")
cat("- ", us_all_output_path("comparison_top2_proposed.csv", subdir = "tables"), "\n", sep = "")
cat("- ", us_all_output_path("comparison_paired_same_basis_proposed.csv", subdir = "tables"), "\n", sep = "")
cat("- ", us_all_output_path("comparison_scalar_metrics_proposed.csv", subdir = "tables"), "\n", sep = "")
cat("- ", us_all_output_path("comparison_classification_metrics_proposed.csv", subdir = "tables"), "\n", sep = "")
cat("- ", us_all_output_path("comparison_brier_split_proposed.csv", subdir = "tables"), "\n", sep = "")
cat("- ", us_all_output_path("comparison_uncertainty_summary_proposed.csv", subdir = "tables"), "\n", sep = "")
cat("- ", us_all_output_path("comparison_calibration_bins_proposed.csv", subdir = "tables"), "\n", sep = "")
