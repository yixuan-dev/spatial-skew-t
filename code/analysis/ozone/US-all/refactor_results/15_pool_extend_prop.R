rm(list = ls())

source("refactor_results/00_bootstrap.R")
source("refactor_results/01_score_engine.R")
source("refactor_results/02_comparison_tables.R")

set_working_dir_to_script()
ensure_us_all_output_dirs()

extend_path <- us_all_output_path("us-all-results-extend.RData",   subdir = "results")
prop_path   <- us_all_output_path("us-all-results-proposed.RData", subdir = "results")
if (!file.exists(extend_path)) stop("Missing extend RData: ", extend_path)
if (!file.exists(prop_path))   stop("Missing proposed RData: ", prop_path)

e <- new.env(); load(extend_path, envir = e)
p <- new.env(); load(prop_path,   envir = p)

extend_score <- e$score_obj
prop_score   <- p$score_obj

stopifnot(identical(extend_score$probs, prop_score$probs))
stopifnot(isTRUE(all.equal(unname(extend_score$thresholds), unname(prop_score$thresholds))))
stopifnot(identical(extend_score$threshold_probs, prop_score$threshold_probs))
stopifnot(extend_score$nsets == prop_score$nsets)
stopifnot(identical(extend_score$brier.split.target_probs, prop_score$brier.split.target_probs))
stopifnot(identical(extend_score$pit_breaks, prop_score$pit_breaks))
stopifnot(identical(extend_score$uncertainty_levels, prop_score$uncertainty_levels))
# Classification arrays may exist in only one side (extend predates classification diagnostics).
# When mismatched, take whichever side has them and fill the other side's slots with NA.
classification_present <- !is.null(prop_score$classification.target_probs) ||
                          !is.null(extend_score$classification.target_probs)
classification_target_probs <- prop_score$classification.target_probs
if (is.null(classification_target_probs)) classification_target_probs <- extend_score$classification.target_probs
classification_target_idx <- prop_score$classification.target_idx
if (is.null(classification_target_idx)) classification_target_idx <- extend_score$classification.target_idx
classification_target_thresholds <- prop_score$classification.target_thresholds
if (is.null(classification_target_thresholds)) classification_target_thresholds <- extend_score$classification.target_thresholds
classification_metric_names <- prop_score$classification.metric_names
if (is.null(classification_metric_names)) classification_metric_names <- extend_score$classification.metric_names
classification_count_names <- prop_score$classification.count_names
if (is.null(classification_count_names)) classification_count_names <- extend_score$classification.count_names
classification_probability_cutoff <- prop_score$classification.probability_cutoff
if (is.null(classification_probability_cutoff)) classification_probability_cutoff <- extend_score$classification.probability_cutoff

# --- ID layout ---------------------------------------------------------------
# Extend settings keep their original IDs (1..max_extend, 1 = Gaussian baseline).
# Prop settings 2..11 (skew-t fits in proposed RData) shift to max_extend_id+1..max_extend_id+10.
# We use max ID from BOTH the extend score arrays AND main settings.csv (which has rows for
# unscored IDs up to 219 too) so prop slots don't collide with reserved IDs.
extend_settings_csv <- read.csv("settings.csv", stringsAsFactors = FALSE)
max_csv <- suppressWarnings(max(as.integer(extend_settings_csv$setting), na.rm = TRUE))
max_extend <- max(extend_score$max_setting,
                  max(extend_score$available_settings, na.rm = TRUE),
                  max_csv)
prop_src_ids <- 2:11
prop_id_base <- 301L  # ozone_prop fits live at 301+ in the pool
if (prop_id_base <= max_extend) {
  stop("prop_id_base (", prop_id_base, ") must exceed max_extend (", max_extend, ")")
}
prop_dst_ids <- prop_id_base + seq_along(prop_src_ids) - 1L
new_max <- max(prop_dst_ids)

# --- helpers -----------------------------------------------------------------
expand_last_dim <- function(x, new_last) {
  if (is.null(x)) return(NULL)
  d <- dim(x)
  if (is.null(d)) {
    if (length(x) >= new_last) return(x)
    out <- rep(NA_real_, new_last); out[seq_along(x)] <- x
    return(out)
  }
  if (d[length(d)] >= new_last) return(x)
  ndims <- length(d)
  new_d <- d; new_d[ndims] <- new_last
  out <- array(NA_real_, dim = new_d, dimnames = dimnames(x))
  pre <- paste(rep(",", ndims - 1), collapse = "")
  eval(parse(text = sprintf("out[%sseq_len(d[ndims])] <- x", pre)))
  out
}

copy_prop_to_unified <- function(extended_arr, prop_arr, src_ids, dst_ids) {
  if (is.null(extended_arr) || is.null(prop_arr)) return(extended_arr)
  d <- dim(extended_arr)
  if (is.null(d)) {
    extended_arr[dst_ids] <- prop_arr[src_ids]
    return(extended_arr)
  }
  ndims <- length(d)
  pre <- paste(rep(",", ndims - 1), collapse = "")
  eval(parse(text = sprintf(
    "extended_arr[%sdst_ids] <- prop_arr[%ssrc_ids]", pre, pre
  )))
  extended_arr
}

merge_setting_array <- function(extend_arr, prop_arr, new_last, src_ids, dst_ids) {
  if (is.null(extend_arr) && is.null(prop_arr)) return(NULL)
  if (is.null(extend_arr)) {
    # Extend side missing this array (e.g., classification on legacy run).
    d_p <- dim(prop_arr)
    if (is.null(d_p)) {
      out <- rep(NA_real_, new_last)
    } else {
      new_d <- d_p
      new_d[length(new_d)] <- new_last
      out <- array(NA_real_, dim = new_d, dimnames = dimnames(prop_arr))
    }
    return(copy_prop_to_unified(out, prop_arr, src_ids, dst_ids))
  }
  out <- expand_last_dim(extend_arr, new_last)
  if (is.null(prop_arr)) return(out)
  copy_prop_to_unified(out, prop_arr, src_ids, dst_ids)
}

# --- merge each setting-axis array ------------------------------------------
score_arrays <- c(
  "quant.score", "brier.score",
  "mspe.score", "mape.score", "crps.score",
  "coverage.score",
  "pit.mean.score", "pit.variance.score", "pit.ks.score", "pit.mae.score", "pit.rmse.score",
  "pit.bin.share", "pit.bin.count",
  "summary.draws.used", "summary.n_obs",
  "brier.split.score", "brier.split.count", "brier.split.share",
  "classification.metric", "classification.count",
  "classification.obs.total",
  "classification.actual_positive_share", "classification.predicted_positive_share",
  "classification.draws.used"
)

merged <- extend_score
merged$max_setting <- new_max

for (nm in score_arrays) {
  merged[[nm]] <- merge_setting_array(
    extend_score[[nm]], prop_score[[nm]],
    new_last = new_max,
    src_ids = prop_src_ids, dst_ids = prop_dst_ids
  )
}

if (classification_present) {
  merged$classification.target_probs       <- classification_target_probs
  merged$classification.target_idx         <- classification_target_idx
  merged$classification.target_thresholds  <- classification_target_thresholds
  merged$classification.metric_names       <- classification_metric_names
  merged$classification.count_names        <- classification_count_names
  merged$classification.probability_cutoff <- classification_probability_cutoff
}

merged$setting_ids <- sort(unique(c(extend_score$setting_ids, prop_dst_ids)))
merged$available_settings <- sort(unique(c(extend_score$available_settings, prop_dst_ids)))
merged$skipped_missing_file  <- sort(unique(c(extend_score$skipped_missing_file,  prop_dst_ids[ !(prop_src_ids %in% prop_score$available_settings) ])))
merged$skipped_bad_contract  <- sort(unique(c(extend_score$skipped_bad_contract,  integer(0))))
merged$skipped_scoring_error <- sort(unique(c(extend_score$skipped_scoring_error, integer(0))))

# --- summarize ---------------------------------------------------------------
score_obj <- merged
summary_obj <- summarize_us_all_scores(score_obj, baseline_setting = 1L)
print_score_summary(score_obj, label = "us-all-results-pool (extend + ozone_prop)")

# --- unified settings table --------------------------------------------------
ctx <- load_us_all_context()
extend_settings <- ctx$settings
required_cols <- c("setting", "method", "knots", "thresh", "CMAQ", "TS", "ar2", "mrts")
for (col in required_cols) {
  if (!(col %in% names(extend_settings))) extend_settings[[col]] <- NA
}
extend_settings$setting_num  <- safe_as_integer(extend_settings$setting)
extend_settings$setting_orig <- extend_settings$setting_num
if (!("source" %in% names(extend_settings))) {
  extend_settings$source <- "ozone_us_all_results"
}

prop_settings_raw <- read.csv(
  file.path("..", "..", "ozone_prop", "US-all", "settings_prop.csv"),
  stringsAsFactors = FALSE
)
prop_settings_raw$setting_num_orig <- safe_as_integer(prop_settings_raw$setting)
prop_settings_raw <- prop_settings_raw[order(prop_settings_raw$setting_num_orig), , drop = FALSE]

prop_unified <- data.frame(
  setting       = prop_dst_ids,
  setting_num   = prop_dst_ids,
  setting_orig  = prop_settings_raw$setting_num_orig,
  method        = prop_settings_raw$method,
  knots         = prop_settings_raw$knots,
  thresh        = prop_settings_raw$thresh,
  CMAQ          = prop_settings_raw$CMAQ,
  TS            = NA_character_,
  ar2           = NA_character_,
  mrts          = as.character(prop_settings_raw$mrts),
  source        = "ozone_prop_us_all_results",
  stringsAsFactors = FALSE
)

shared_cols <- intersect(names(extend_settings), names(prop_unified))
extend_keep  <- extend_settings[, shared_cols, drop = FALSE]
prop_keep    <- prop_unified[,   shared_cols, drop = FALSE]
unified_settings <- rbind(extend_keep, prop_keep)

# --- proposed/baseline groupings --------------------------------------------
done_morris      <- e$done_morris
done_ar2         <- e$done_ar2
done_mrts        <- e$done_mrts
done_extensions  <- e$done_extensions
done_prop        <- prop_dst_ids

baseline_ids <- done_morris
proposed_ids <- sort(unique(c(done_extensions, done_prop)))

decorate_setting_table <- function(df, setting_col = "setting") {
  if (nrow(df) == 0) return(df)
  idx <- match(df[[setting_col]], unified_settings$setting_num)
  df$mrts         <- unified_settings$mrts[idx]
  df$setting_orig <- unified_settings$setting_orig[idx]
  df$source       <- unified_settings$source[idx]
  df$is_baseline  <- df[[setting_col]] %in% baseline_ids
  df$is_proposed  <- df[[setting_col]] %in% proposed_ids
  df$is_ar2       <- df[[setting_col]] %in% done_ar2
  df$is_mrts      <- df[[setting_col]] %in% done_mrts
  df$is_prop      <- df[[setting_col]] %in% done_prop
  df$model_lane <- ifelse(
    df[[setting_col]] == 1L, "gaussian_reference",
    ifelse(df$is_prop, "ozone_prop",
    ifelse(df$is_ar2,  "ar2",
    ifelse(df$is_mrts, "mrts",
    ifelse(df[[setting_col]] %in% done_morris, "morris_baseline", "other_numeric"))))
  )
  df
}

# --- comparison tables -------------------------------------------------------
comparison_full_table <- build_comparison_full_table(
  summary_obj = summary_obj, settings = unified_settings,
  baseline_ids = baseline_ids, proposed_ids = proposed_ids
)
if (nrow(comparison_full_table) > 0) {
  comparison_full_table <- decorate_setting_table(comparison_full_table, "setting")
}

comparison_top2 <- build_comparison_top2_all_metrics(
  summary_obj = summary_obj, settings = unified_settings,
  candidate_settings = summary_obj$available_settings
)
if (nrow(comparison_top2) > 0) {
  comparison_top2 <- decorate_setting_table(comparison_top2, "setting")
}

comparison_paired_same_basis <- build_paired_same_basis_table(
  summary_obj = summary_obj, settings = unified_settings,
  baseline_ids = baseline_ids, proposed_ids = proposed_ids
)
if (nrow(comparison_paired_same_basis) > 0) {
  comparison_paired_same_basis$proposed_is_ar2  <- comparison_paired_same_basis$proposed_setting %in% done_ar2
  comparison_paired_same_basis$proposed_is_mrts <- comparison_paired_same_basis$proposed_setting %in% done_mrts
  comparison_paired_same_basis$proposed_is_prop <- comparison_paired_same_basis$proposed_setting %in% done_prop
  comparison_paired_same_basis$proposed_lane <- ifelse(
    comparison_paired_same_basis$proposed_is_prop, "ozone_prop",
    ifelse(comparison_paired_same_basis$proposed_is_ar2, "ar2",
    ifelse(comparison_paired_same_basis$proposed_is_mrts, "mrts", "other"))
  )
}

comparison_scalar_metrics <- build_comparison_scalar_metrics_table(
  summary_obj = summary_obj, settings = unified_settings,
  baseline_ids = baseline_ids, proposed_ids = proposed_ids
)
if (nrow(comparison_scalar_metrics) > 0) {
  comparison_scalar_metrics <- decorate_setting_table(comparison_scalar_metrics, "setting")
}

comparison_classification_metrics <- build_comparison_classification_metrics_table(
  summary_obj = summary_obj, settings = unified_settings,
  baseline_ids = baseline_ids, proposed_ids = proposed_ids
)
if (nrow(comparison_classification_metrics) > 0) {
  comparison_classification_metrics <- decorate_setting_table(comparison_classification_metrics, "setting")
}

comparison_brier_split <- build_comparison_brier_split_table(
  summary_obj = summary_obj, settings = unified_settings,
  baseline_ids = baseline_ids, proposed_ids = proposed_ids
)
if (nrow(comparison_brier_split) > 0) {
  comparison_brier_split <- decorate_setting_table(comparison_brier_split, "setting")
}

comparison_uncertainty_summary <- build_comparison_uncertainty_summary_table(
  summary_obj = summary_obj, settings = unified_settings,
  baseline_ids = baseline_ids, proposed_ids = proposed_ids
)
if (nrow(comparison_uncertainty_summary) > 0) {
  comparison_uncertainty_summary <- decorate_setting_table(comparison_uncertainty_summary, "setting")
}

comparison_calibration_bins <- build_comparison_calibration_bins_table(
  summary_obj = summary_obj, settings = unified_settings,
  baseline_ids = baseline_ids, proposed_ids = proposed_ids
)
if (nrow(comparison_calibration_bins) > 0) {
  comparison_calibration_bins <- decorate_setting_table(comparison_calibration_bins, "setting")
}

# --- write CSVs --------------------------------------------------------------
write.csv(comparison_full_table,             us_all_output_path("comparison_full_table_pool.csv",             subdir = "tables"), row.names = FALSE)
write.csv(comparison_top2,                   us_all_output_path("comparison_top2_pool.csv",                   subdir = "tables"), row.names = FALSE)
write_comparison_top2_workbook(
  comparison_top2 = comparison_top2,
  output_path = us_all_output_path("comparison_top2_pool.xlsx", subdir = "tables")
)
write.csv(comparison_paired_same_basis,      us_all_output_path("comparison_paired_same_basis_pool.csv",      subdir = "tables"), row.names = FALSE)
write.csv(comparison_scalar_metrics,         us_all_output_path("comparison_scalar_metrics_pool.csv",         subdir = "tables"), row.names = FALSE)
write.csv(comparison_classification_metrics, us_all_output_path("comparison_classification_metrics_pool.csv", subdir = "tables"), row.names = FALSE)
write.csv(comparison_brier_split,            us_all_output_path("comparison_brier_split_pool.csv",            subdir = "tables"), row.names = FALSE)
write.csv(comparison_uncertainty_summary,    us_all_output_path("comparison_uncertainty_summary_pool.csv",    subdir = "tables"), row.names = FALSE)
write.csv(comparison_calibration_bins,       us_all_output_path("comparison_calibration_bins_pool.csv",       subdir = "tables"), row.names = FALSE)

# --- save pooled RData -------------------------------------------------------
settings              <- unified_settings
available_settings    <- summary_obj$available_settings
skipped_missing_file  <- summary_obj$skipped_missing_file
skipped_bad_contract  <- summary_obj$skipped_bad_contract
skipped_scoring_error <- summary_obj$skipped_scoring_error

probs <- summary_obj$probs
threshold_probs <- summary_obj$threshold_probs
thresholds <- summary_obj$thresholds

quant.score <- summary_obj$quant.score
brier.score <- summary_obj$brier.score
mspe.score  <- summary_obj$mspe.score
mape.score  <- summary_obj$mape.score
quant.score.mean <- summary_obj$quant.score.mean
brier.score.mean <- summary_obj$brier.score.mean
mspe.mean <- summary_obj$mspe.mean
mape.mean <- summary_obj$mape.mean
quant.score.se <- summary_obj$quant.score.se
brier.score.se <- summary_obj$brier.score.se
mspe.se <- summary_obj$mspe.se
mape.se <- summary_obj$mape.se
bs.mean.ref.gau   <- summary_obj$bs.mean.ref.gau
qs.mean.ref.gau   <- summary_obj$qs.mean.ref.gau
mspe.mean.ref.gau <- summary_obj$mspe.mean.ref.gau
mape.mean.ref.gau <- summary_obj$mape.mean.ref.gau
crps.mean         <- summary_obj$crps.mean
crps.se           <- summary_obj$crps.se
crps.mean.ref.gau <- summary_obj$crps.mean.ref.gau
brier.split.target_probs      <- summary_obj$brier.split.target_probs
brier.split.target_thresholds <- summary_obj$brier.split.target_thresholds
brier.split.band_names        <- summary_obj$brier.split.band_names
brier.split.score.mean        <- summary_obj$brier.split.score.mean
brier.split.score.se          <- summary_obj$brier.split.score.se
brier.split.n_obs.total       <- summary_obj$brier.split.n_obs.total
brier.split.n_obs.mean        <- summary_obj$brier.split.n_obs.mean
brier.split.obs.share         <- summary_obj$brier.split.obs.share
brier.split.rel.ref.gau       <- summary_obj$brier.split.rel.ref.gau
classification.target_probs        <- summary_obj$classification.target_probs
classification.target_thresholds   <- summary_obj$classification.target_thresholds
classification.metric.names        <- summary_obj$classification.metric_names
classification.count.names         <- summary_obj$classification.count_names
classification.metric.mean         <- summary_obj$classification.metric.mean
classification.metric.se           <- summary_obj$classification.metric.se
classification.count.total         <- summary_obj$classification.count.total
classification.count.mean          <- summary_obj$classification.count.mean
classification.n_obs.total         <- summary_obj$classification.n_obs.total
classification.n_obs.mean          <- summary_obj$classification.n_obs.mean
classification.actual_positive_share    <- summary_obj$classification.actual_positive_share
classification.predicted_positive_share <- summary_obj$classification.predicted_positive_share
classification.metric.rel.ref.gau  <- summary_obj$classification.metric.rel.ref.gau
classification.metric.delta.ref.gau <- summary_obj$classification.metric.delta.ref.gau
classification.probability_cutoff  <- summary_obj$classification.probability_cutoff
coverage.mean <- summary_obj$coverage.mean
coverage.se   <- summary_obj$coverage.se
coverage.gap  <- summary_obj$coverage.gap
pit.mean      <- summary_obj$pit.mean
pit.variance  <- summary_obj$pit.variance
pit.ks        <- summary_obj$pit.ks
pit.mae       <- summary_obj$pit.mae
pit.rmse      <- summary_obj$pit.rmse
pit.bin.share.mean <- summary_obj$pit.bin.share.mean
pit.bin.share.se   <- summary_obj$pit.bin.share.se
summary.draws.mean   <- summary_obj$summary.draws.mean
summary.n_obs.total  <- summary_obj$summary.n_obs.total

save(
  list = c(
    "settings", "baseline_ids", "proposed_ids",
    "done_morris", "done_ar2", "done_mrts", "done_extensions", "done_prop",
    "available_settings", "skipped_missing_file", "skipped_bad_contract", "skipped_scoring_error",
    "probs", "threshold_probs", "thresholds",
    "quant.score", "brier.score", "mspe.score", "mape.score",
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
  file = us_all_output_path("us-all-results-pool.RData", subdir = "results")
)

cat("\nPool layout:\n")
cat("- Extend settings: 1..", max_extend, " (kept as-is, ", length(e$available_settings), " scored, baseline=1)\n", sep = "")
cat("- Prop settings:   ", min(prop_dst_ids), "..", max(prop_dst_ids), " (mapped from ozone_prop ", paste(prop_src_ids, collapse=","), ")\n", sep = "")
cat("- Pool max_setting:", new_max, "\n")
cat("- Pool available:  ", length(available_settings), "settings\n")

cat("\nOutputs written:\n")
cat("- ", us_all_output_path("us-all-results-pool.RData", subdir = "results"), "\n", sep = "")
for (suff in c("comparison_full_table_pool.csv","comparison_top2_pool.csv","comparison_top2_pool.xlsx",
               "comparison_paired_same_basis_pool.csv","comparison_scalar_metrics_pool.csv",
               "comparison_classification_metrics_pool.csv","comparison_brier_split_pool.csv",
               "comparison_uncertainty_summary_pool.csv","comparison_calibration_bins_pool.csv")) {
  cat("- ", us_all_output_path(suff, subdir = "tables"), "\n", sep = "")
}
