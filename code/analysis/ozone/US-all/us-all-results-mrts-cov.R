rm(list = ls())

source("refactor_results/00_bootstrap.R")
source("refactor_results/01_score_engine.R")
source("refactor_results/02_comparison_tables.R")

set_working_dir_to_script()
source("mrts_cov_helpers.R")
ctx <- load_us_all_context()

pick_results_dir <- function() {
  from_env <- trimws(Sys.getenv("US_ALL_RESULTS_DIR", unset = ""))
  if (nzchar(from_env)) {
    return(from_env)
  }

  if (dir.exists("results_new")) {
    return("results_new")
  }

  "results"
}

results_dir <- pick_results_dir()
if (!dir.exists(results_dir)) {
  stop("Results directory not found: ", results_dir)
}

target_pairs <- build_mrts_target_pairs(ctx$settings)
baseline_settings <- sort(unique(target_pairs$baseline_setting))
mrts_settings <- sort(unique(target_pairs$mrts_setting))
requested_settings <- sort(unique(c(1L, baseline_settings, mrts_settings)))

probs <- default_probs(include_999 = FALSE)
thresholds <- quantile(ctx$Y, probs = probs, na.rm = TRUE)

score_obj <- compute_us_all_scores(
  setting_ids = requested_settings,
  result_path_fn = function(i) file.path(results_dir, sprintf("us-all-%d.RData", i)),
  Y = ctx$Y,
  cv_lst = ctx$cv_lst,
  probs = probs,
  thresholds = thresholds,
  trans_setting_ids = 2L,
  enforce_contract = TRUE,
  partial_ok_setting_ids = 210L
)
summary_obj <- summarize_us_all_scores(score_obj, baseline_setting = 1L)
print_score_summary(score_obj, label = "us-all-results-mrts-cov")

comparison_mrts_cov_full_table <- build_comparison_full_table(
  summary_obj = summary_obj,
  settings = ctx$settings,
  baseline_ids = baseline_settings,
  proposed_ids = mrts_settings
)

settings_lookup <- ctx$settings
if (!("setting_num" %in% names(settings_lookup)) && "setting" %in% names(settings_lookup)) {
  settings_lookup$setting_num <- safe_as_integer(settings_lookup$setting)
}

add_meta_columns <- function(df, settings_lookup) {
  if (nrow(df) == 0) {
    return(df)
  }

  idx <- match(df$setting, settings_lookup$setting_num)
  if ("mrts" %in% names(settings_lookup)) {
    df$mrts <- settings_lookup$mrts[idx]
  }
  if ("ar2" %in% names(settings_lookup)) {
    df$ar2 <- settings_lookup$ar2[idx]
  }
  df
}

comparison_mrts_cov_full_table <- add_meta_columns(comparison_mrts_cov_full_table, settings_lookup)

candidate_settings <- sort(unique(c(baseline_settings, mrts_settings)))
comparison_mrts_cov_top2 <- build_comparison_top2(
  summary_obj = summary_obj,
  settings = ctx$settings,
  target_levels = c(0.95, 0.98, 0.99, 0.995),
  candidate_settings = candidate_settings,
  metric = "brier"
)
if (nrow(comparison_mrts_cov_top2) > 0 && "setting_num" %in% names(settings_lookup) && "mrts" %in% names(settings_lookup)) {
  idx <- match(comparison_mrts_cov_top2$setting, settings_lookup$setting_num)
  comparison_mrts_cov_top2$mrts <- settings_lookup$mrts[idx]
}

paired_rows <- list()
row_id <- 1

for (k in seq_len(nrow(target_pairs))) {
  base_setting <- target_pairs$baseline_setting[k]
  mrts_setting <- target_pairs$mrts_setting[k]

  if (!(base_setting %in% summary_obj$available_settings && mrts_setting %in% summary_obj$available_settings)) {
    next
  }

  for (q in seq_along(summary_obj$probs)) {
    base_rel_brier <- summary_obj$bs.mean.ref.gau[base_setting, q]
    mrts_rel_brier <- summary_obj$bs.mean.ref.gau[mrts_setting, q]
    base_abs_brier <- summary_obj$brier.score.mean[base_setting, q]
    mrts_abs_brier <- summary_obj$brier.score.mean[mrts_setting, q]

    base_rel_quant <- summary_obj$qs.mean.ref.gau[base_setting, q]
    mrts_rel_quant <- summary_obj$qs.mean.ref.gau[mrts_setting, q]
    base_abs_quant <- summary_obj$quant.score.mean[base_setting, q]
    mrts_abs_quant <- summary_obj$quant.score.mean[mrts_setting, q]

    improve_rel_brier_pct <- if (is.finite(base_rel_brier) && base_rel_brier != 0) {
      100 * (base_rel_brier - mrts_rel_brier) / base_rel_brier
    } else {
      NA_real_
    }

    improve_abs_brier_pct <- if (is.finite(base_abs_brier) && base_abs_brier != 0) {
      100 * (base_abs_brier - mrts_abs_brier) / base_abs_brier
    } else {
      NA_real_
    }

    improve_rel_quant_pct <- if (is.finite(base_rel_quant) && base_rel_quant != 0) {
      100 * (base_rel_quant - mrts_rel_quant) / base_rel_quant
    } else {
      NA_real_
    }

    improve_abs_quant_pct <- if (is.finite(base_abs_quant) && base_abs_quant != 0) {
      100 * (base_abs_quant - mrts_abs_quant) / base_abs_quant
    } else {
      NA_real_
    }

    paired_rows[[row_id]] <- data.frame(
      pair_id = target_pairs$pair_id[k],
      quantile = summary_obj$probs[q],
      method = target_pairs$method[k],
      knots = target_pairs$knots[k],
      thresh = target_pairs$thresh[k],
      baseline_setting = base_setting,
      mrts_setting = mrts_setting,
      mrts_k = target_pairs$mrts_k[k],
      baseline_rel_brier = base_rel_brier,
      mrts_rel_brier = mrts_rel_brier,
      delta_rel_brier = mrts_rel_brier - base_rel_brier,
      improve_rel_brier_pct = improve_rel_brier_pct,
      baseline_abs_brier = base_abs_brier,
      mrts_abs_brier = mrts_abs_brier,
      delta_abs_brier = mrts_abs_brier - base_abs_brier,
      improve_abs_brier_pct = improve_abs_brier_pct,
      baseline_rel_quant = base_rel_quant,
      mrts_rel_quant = mrts_rel_quant,
      delta_rel_quant = mrts_rel_quant - base_rel_quant,
      improve_rel_quant_pct = improve_rel_quant_pct,
      baseline_abs_quant = base_abs_quant,
      mrts_abs_quant = mrts_abs_quant,
      delta_abs_quant = mrts_abs_quant - base_abs_quant,
      improve_abs_quant_pct = improve_abs_quant_pct,
      stringsAsFactors = FALSE
    )
    row_id <- row_id + 1
  }
}

comparison_mrts_cov_paired <- if (length(paired_rows) > 0) do.call(rbind, paired_rows) else data.frame()

if (nrow(comparison_mrts_cov_paired) > 0) {
  split_groups <- split(comparison_mrts_cov_paired, list(comparison_mrts_cov_paired$mrts_k, comparison_mrts_cov_paired$quantile), drop = TRUE)

  summary_rows <- list()
  sid <- 1
  for (nm in names(split_groups)) {
    grp <- split_groups[[nm]]

    summary_rows[[sid]] <- data.frame(
      mrts_k = grp$mrts_k[1],
      quantile = grp$quantile[1],
      n_pairs = nrow(grp),
      n_improved_rel_brier = sum(grp$delta_rel_brier < 0, na.rm = TRUE),
      mean_delta_rel_brier = mean(grp$delta_rel_brier, na.rm = TRUE),
      mean_improve_rel_brier_pct = mean(grp$improve_rel_brier_pct, na.rm = TRUE),
      mean_delta_abs_brier = mean(grp$delta_abs_brier, na.rm = TRUE),
      mean_improve_abs_brier_pct = mean(grp$improve_abs_brier_pct, na.rm = TRUE),
      mean_delta_rel_quant = mean(grp$delta_rel_quant, na.rm = TRUE),
      mean_improve_rel_quant_pct = mean(grp$improve_rel_quant_pct, na.rm = TRUE),
      mean_delta_abs_quant = mean(grp$delta_abs_quant, na.rm = TRUE),
      mean_improve_abs_quant_pct = mean(grp$improve_abs_quant_pct, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    sid <- sid + 1
  }

  comparison_mrts_cov_summary_by_k <- do.call(rbind, summary_rows)
  comparison_mrts_cov_summary_by_k <- comparison_mrts_cov_summary_by_k[
    order(comparison_mrts_cov_summary_by_k$mrts_k, comparison_mrts_cov_summary_by_k$quantile),
    ,
    drop = FALSE
  ]
  rownames(comparison_mrts_cov_summary_by_k) <- NULL
} else {
  comparison_mrts_cov_summary_by_k <- data.frame()
}

tables_dir <- "D:/Github/spatial-skew-t/code/analysis/ozone/US-all/output/us-all/tables"
if (!dir.exists(tables_dir)) {
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
}

write.csv(target_pairs, file.path(tables_dir, "mrts_cov_target_settings.csv"), row.names = FALSE)
write.csv(comparison_mrts_cov_full_table, file.path(tables_dir, "comparison_mrts_cov_full_table.csv"), row.names = FALSE)
write.csv(comparison_mrts_cov_top2, file.path(tables_dir, "comparison_mrts_cov_top2.csv"), row.names = FALSE)
write.csv(comparison_mrts_cov_paired, file.path(tables_dir, "comparison_mrts_cov_paired.csv"), row.names = FALSE)
write.csv(comparison_mrts_cov_summary_by_k, file.path(tables_dir, "comparison_mrts_cov_summary_by_k.csv"), row.names = FALSE)

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

save(
  list = c(
    "results_dir", "target_pairs", "baseline_settings", "mrts_settings", "requested_settings",
    "available_settings", "skipped_missing_file", "skipped_bad_contract", "skipped_scoring_error",
    "probs", "thresholds",
    "quant.score", "brier.score",
    "quant.score.mean", "brier.score.mean", "quant.score.se", "brier.score.se",
    "bs.mean.ref.gau", "qs.mean.ref.gau",
    "comparison_mrts_cov_full_table", "comparison_mrts_cov_top2",
    "comparison_mrts_cov_paired", "comparison_mrts_cov_summary_by_k",
    "score_obj", "summary_obj"
  ),
  file = "us-all-results-mrts-cov.RData"
)

cat("\n=== us-all-results-mrts-cov summary ===\n")
cat("Results dir:", results_dir, "\n")
cat("Tables dir:", tables_dir, "\n")
cat("Baseline settings:", paste(baseline_settings, collapse = ", "), "\n")
cat("MRTS settings:", paste(mrts_settings, collapse = ", "), "\n")
cat("Requested settings:", paste(requested_settings, collapse = ", "), "\n")
cat("Outputs written:\n")
cat("- [tables]/mrts_cov_target_settings.csv\n")
cat("- [tables]/comparison_mrts_cov_full_table.csv\n")
cat("- [tables]/comparison_mrts_cov_top2.csv\n")
cat("- [tables]/comparison_mrts_cov_paired.csv\n")
cat("- [tables]/comparison_mrts_cov_summary_by_k.csv\n")
cat("- us-all-results-mrts-cov.RData\n")
