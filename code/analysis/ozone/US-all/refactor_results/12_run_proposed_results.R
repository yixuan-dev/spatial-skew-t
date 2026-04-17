rm(list = ls())

source("refactor_results/00_bootstrap.R")
source("refactor_results/01_score_engine.R")
source("refactor_results/02_comparison_tables.R")

set_working_dir_to_script()
ctx <- load_us_all_context()
ensure_us_all_output_dirs()

normalize_yes_no <- function(x) {
  txt <- tolower(trimws(as.character(x)))
  txt[!nzchar(txt)] <- NA_character_
  txt
}

parse_mrts_int <- function(x) {
  raw <- trimws(as.character(x))
  raw[!nzchar(raw)] <- NA_character_
  num <- suppressWarnings(as.numeric(raw))
  out <- suppressWarnings(as.integer(round(num)))
  out[!is.finite(num)] <- NA_integer_
  out
}

expand_range_tokens <- function(tokens) {
  out <- integer(0)
  for (tok in tokens) {
    tok <- trimws(tok)
    if (!nzchar(tok)) {
      next
    }

    if (grepl("^[0-9]+:[0-9]+$", tok)) {
      parts <- as.integer(strsplit(tok, ":", fixed = TRUE)[[1]])
      if (length(parts) == 2 && !any(is.na(parts))) {
        if (parts[1] <= parts[2]) {
          out <- c(out, seq(parts[1], parts[2]))
        } else {
          out <- c(out, seq(parts[1], parts[2], by = -1))
        }
      }
    } else {
      val <- suppressWarnings(as.integer(tok))
      if (!is.na(val)) {
        out <- c(out, val)
      }
    }
  }
  sort(unique(out))
}

pick_results_dirs <- function() {
  dirs_env <- trimws(Sys.getenv("US_ALL_RESULTS_DIRS", unset = ""))
  dir_primary <- trimws(Sys.getenv("US_ALL_RESULTS_DIR", unset = ""))

  dirs <- character(0)

  if (nzchar(dirs_env)) {
    dirs <- c(dirs, unlist(strsplit(dirs_env, ",", fixed = TRUE)))
  }
  if (nzchar(dir_primary)) {
    dirs <- c(dirs, dir_primary)
  }

  # Default search order: legacy/main then MRTS lane.
  dirs <- c(dirs, "results", "results_new")
  dirs <- trimws(dirs)
  dirs <- unique(dirs[nzchar(dirs)])

  existing <- dirs[dir.exists(dirs)]
  if (length(existing) == 0) {
    stop("No result directory found. Checked: ", paste(dirs, collapse = ", "))
  }

  existing
}

resolve_result_file <- function(setting_id, result_dirs) {
  candidates <- file.path(result_dirs, sprintf("us-all-%d.RData", setting_id))
  hit <- candidates[file.exists(candidates)]
  if (length(hit) > 0) {
    return(hit[1])
  }

  # Return first candidate so scorer records a clean missing-file path.
  candidates[1]
}

settings <- ctx$settings
if (!("setting_num" %in% names(settings)) && "setting" %in% names(settings)) {
  settings$setting_num <- safe_as_integer(settings$setting)
}

if (!("setting_num" %in% names(settings))) {
  stop("settings.csv must provide a numeric setting namespace (column: setting or setting_num).")
}

all_numeric_settings <- sort(unique(settings$setting_num[!is.na(settings$setting_num)]))
if (length(all_numeric_settings) == 0) {
  stop("No numeric settings found in settings.csv")
}

# Keep the Morris baseline lane frozen for historical comparability flags.
done_morris <- c(1:5, 7:9, 11:13, 15:17, 33:36, 38:41, 43:46, 51:74)

done_ar2 <- integer(0)
if ("ar2" %in% names(settings)) {
  ar2_flag <- normalize_yes_no(settings$ar2)
  done_ar2 <- sort(unique(settings$setting_num[!is.na(settings$setting_num) & ar2_flag == "yes"]))
}

done_mrts <- integer(0)
if ("mrts" %in% names(settings)) {
  mrts_k <- parse_mrts_int(settings$mrts)
  done_mrts <- sort(unique(settings$setting_num[!is.na(settings$setting_num) & !is.na(mrts_k) & mrts_k > 0]))
}

done_extensions <- sort(unique(c(done_ar2, done_mrts)))
all_requested <- all_numeric_settings

# Optional quick filter for debugging/smoke checks.
# Example: US_ALL_COMPARE_SETTINGS="1,51,101,201" or "1:20,101:124,201:209"
compare_filter <- trimws(Sys.getenv("US_ALL_COMPARE_SETTINGS", unset = ""))
if (nzchar(compare_filter)) {
  filter_ids <- expand_range_tokens(unlist(strsplit(compare_filter, ",", fixed = TRUE)))
  if (length(filter_ids) > 0) {
    all_requested <- sort(unique(intersect(all_requested, filter_ids)))
  }
}

# Relative metrics require Gaussian reference.
if (!1L %in% all_requested) {
  all_requested <- sort(unique(c(1L, all_requested)))
}

result_dirs <- pick_results_dirs()

result_file_map <- vapply(
  all_requested,
  function(i) resolve_result_file(i, result_dirs = result_dirs),
  character(1)
)
names(result_file_map) <- as.character(all_requested)

summary_draws <- suppressWarnings(as.integer(Sys.getenv("US_ALL_SUMMARY_DRAWS", unset = "400")))
if (!is.finite(summary_draws) || summary_draws <= 0) {
  summary_draws <- 400L
}

decorate_setting_table <- function(df, setting_col = "setting") {
  if (nrow(df) == 0) {
    return(df)
  }

  idx <- match(df[[setting_col]], settings$setting_num)
  baseline_flag <- if ("is_baseline" %in% names(df)) df$is_baseline else df[[setting_col]] %in% done_morris
  df$mrts <- if ("mrts" %in% names(settings)) settings$mrts[idx] else NA
  df$is_ar2 <- df[[setting_col]] %in% done_ar2
  df$is_mrts <- df[[setting_col]] %in% done_mrts
  df$model_lane <- ifelse(
    df[[setting_col]] == 1,
    "gaussian_reference",
    ifelse(
      df$is_ar2,
      "ar2",
      ifelse(
        df$is_mrts,
        "mrts",
        ifelse(baseline_flag, "morris_baseline", "other_numeric")
      )
    )
  )

  df
}

probs <- default_probs(include_999 = FALSE)
threshold_probs <- default_threshold_probs()
thresholds <- quantile(ctx$Y, probs = threshold_probs, na.rm = TRUE)

score_obj <- compute_us_all_scores(
  setting_ids = all_requested,
  result_path_fn = function(i) result_file_map[[as.character(i)]],
  Y = ctx$Y,
  cv_lst = ctx$cv_lst,
  probs = probs,
  thresholds = thresholds,
  threshold_probs = threshold_probs,
  compute_brier_split_diagnostics = TRUE,
  trans_setting_ids = 2L,
  enforce_contract = TRUE,
  compute_uncertainty_diagnostics = TRUE,
  summary_draws = summary_draws
)
summary_obj <- summarize_us_all_scores(score_obj, baseline_setting = 1L)
print_score_summary(score_obj, label = "us-all-results-proposed (all settings incl. AR2 + MRTS)")

comparison_full_table <- build_comparison_full_table(
  summary_obj = summary_obj,
  settings = settings,
  baseline_ids = done_morris,
  proposed_ids = done_extensions
)

if (nrow(comparison_full_table) > 0) {
  comparison_full_table <- decorate_setting_table(comparison_full_table, setting_col = "setting")
}

comparison_top2 <- build_comparison_top2(
  summary_obj = summary_obj,
  settings = settings,
  target_levels = c(0.95, 0.98, 0.99, 0.995),
  candidate_settings = summary_obj$available_settings,
  metric = "brier"
)

if (nrow(comparison_top2) > 0) {
  comparison_top2 <- decorate_setting_table(comparison_top2, setting_col = "setting")
}

comparison_paired_same_basis <- build_paired_same_basis_table(
  summary_obj = summary_obj,
  settings = settings,
  baseline_ids = done_morris,
  proposed_ids = done_extensions
)

comparison_scalar_metrics <- build_comparison_scalar_metrics_table(
  summary_obj = summary_obj,
  settings = settings,
  baseline_ids = done_morris,
  proposed_ids = done_extensions
)

comparison_brier_split <- build_comparison_brier_split_table(
  summary_obj = summary_obj,
  settings = settings,
  baseline_ids = done_morris,
  proposed_ids = done_extensions
)

comparison_uncertainty_summary <- build_comparison_uncertainty_summary_table(
  summary_obj = summary_obj,
  settings = settings,
  baseline_ids = done_morris,
  proposed_ids = done_extensions
)

comparison_calibration_bins <- build_comparison_calibration_bins_table(
  summary_obj = summary_obj,
  settings = settings,
  baseline_ids = done_morris,
  proposed_ids = done_extensions
)

if (nrow(comparison_paired_same_basis) > 0) {
  comparison_paired_same_basis$proposed_is_ar2 <- comparison_paired_same_basis$proposed_setting %in% done_ar2
  comparison_paired_same_basis$proposed_is_mrts <- comparison_paired_same_basis$proposed_setting %in% done_mrts
  comparison_paired_same_basis$proposed_lane <- ifelse(
    comparison_paired_same_basis$proposed_is_ar2,
    "ar2",
    ifelse(comparison_paired_same_basis$proposed_is_mrts, "mrts", "other")
  )
}

if (nrow(comparison_scalar_metrics) > 0) {
  comparison_scalar_metrics <- decorate_setting_table(comparison_scalar_metrics, setting_col = "setting")
}

if (nrow(comparison_brier_split) > 0) {
  comparison_brier_split <- decorate_setting_table(comparison_brier_split, setting_col = "setting")
}

if (nrow(comparison_uncertainty_summary) > 0) {
  comparison_uncertainty_summary <- decorate_setting_table(comparison_uncertainty_summary, setting_col = "setting")
}

if (nrow(comparison_calibration_bins) > 0) {
  comparison_calibration_bins <- decorate_setting_table(comparison_calibration_bins, setting_col = "setting")
}

write.csv(comparison_full_table, us_all_output_path("comparison_full_table.csv", subdir = "tables"), row.names = FALSE)
write.csv(comparison_top2, us_all_output_path("comparison_top2.csv", subdir = "tables"), row.names = FALSE)
write.csv(comparison_paired_same_basis, us_all_output_path("comparison_paired_same_basis.csv", subdir = "tables"), row.names = FALSE)
write.csv(comparison_scalar_metrics, us_all_output_path("comparison_scalar_metrics.csv", subdir = "tables"), row.names = FALSE)
write.csv(comparison_brier_split, us_all_output_path("comparison_brier_split.csv", subdir = "tables"), row.names = FALSE)
write.csv(comparison_uncertainty_summary, us_all_output_path("comparison_uncertainty_summary.csv", subdir = "tables"), row.names = FALSE)
write.csv(comparison_calibration_bins, us_all_output_path("comparison_calibration_bins.csv", subdir = "tables"), row.names = FALSE)

available_settings <- summary_obj$available_settings
skipped_missing_file <- summary_obj$skipped_missing_file
skipped_bad_contract <- summary_obj$skipped_bad_contract
skipped_scoring_error <- summary_obj$skipped_scoring_error
quant.score <- summary_obj$quant.score
brier.score <- summary_obj$brier.score
quant.score.mean <- summary_obj$quant.score.mean
brier.score.mean <- summary_obj$brier.score.mean
quant.score.se <- summary_obj$quant.score.se
brier.score.se <- summary_obj$brier.score.se
bs.mean.ref.gau <- summary_obj$bs.mean.ref.gau
qs.mean.ref.gau <- summary_obj$qs.mean.ref.gau
brier.split.target_probs <- summary_obj$brier.split.target_probs
brier.split.target_thresholds <- summary_obj$brier.split.target_thresholds
brier.split.band_names <- summary_obj$brier.split.band_names
brier.split.score.mean <- summary_obj$brier.split.score.mean
brier.split.score.se <- summary_obj$brier.split.score.se
brier.split.n_obs.total <- summary_obj$brier.split.n_obs.total
brier.split.n_obs.mean <- summary_obj$brier.split.n_obs.mean
brier.split.obs.share <- summary_obj$brier.split.obs.share
brier.split.rel.ref.gau <- summary_obj$brier.split.rel.ref.gau
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

save(
  list = c(
    "done_morris", "done_ar2", "done_mrts", "done_extensions",
    "all_numeric_settings", "all_requested", "result_dirs", "result_file_map", "available_settings",
    "skipped_missing_file", "skipped_bad_contract", "skipped_scoring_error",
    "probs", "threshold_probs", "thresholds", "quant.score", "brier.score",
    "quant.score.mean", "brier.score.mean", "quant.score.se", "brier.score.se",
    "bs.mean.ref.gau", "qs.mean.ref.gau",
    "brier.split.target_probs", "brier.split.target_thresholds", "brier.split.band_names",
    "brier.split.score.mean", "brier.split.score.se", "brier.split.n_obs.total",
    "brier.split.n_obs.mean", "brier.split.obs.share", "brier.split.rel.ref.gau",
    "crps.mean", "crps.se", "crps.mean.ref.gau",
    "coverage.mean", "coverage.se", "coverage.gap",
    "pit.mean", "pit.variance", "pit.ks", "pit.mae", "pit.rmse",
    "pit.bin.share.mean", "pit.bin.share.se",
    "summary.draws.mean", "summary.n_obs.total",
    "comparison_full_table", "comparison_top2", "comparison_paired_same_basis",
    "comparison_scalar_metrics", "comparison_brier_split", "comparison_uncertainty_summary", "comparison_calibration_bins",
    "score_obj", "summary_obj"
  ),
  file = us_all_output_path("us-all-results-proposed.RData", subdir = "results")
)

cat("Result directories searched (priority order):\n")
for (d in result_dirs) {
  cat("-", d, "\n")
}

cat("Outputs written:\n")
cat("- ", us_all_output_path("us-all-results-proposed.RData", subdir = "results"), "\n", sep = "")
cat("- ", us_all_output_path("comparison_full_table.csv", subdir = "tables"), "\n", sep = "")
cat("- ", us_all_output_path("comparison_top2.csv", subdir = "tables"), "\n", sep = "")
cat("- ", us_all_output_path("comparison_paired_same_basis.csv", subdir = "tables"), "\n", sep = "")
cat("- ", us_all_output_path("comparison_scalar_metrics.csv", subdir = "tables"), "\n", sep = "")
cat("- ", us_all_output_path("comparison_brier_split.csv", subdir = "tables"), "\n", sep = "")
cat("- ", us_all_output_path("comparison_uncertainty_summary.csv", subdir = "tables"), "\n", sep = "")
cat("- ", us_all_output_path("comparison_calibration_bins.csv", subdir = "tables"), "\n", sep = "")
