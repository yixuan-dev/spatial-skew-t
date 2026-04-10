rm(list = ls())

source("refactor_results/00_bootstrap.R")
source("refactor_results/01_score_engine.R")
source("refactor_results/02_comparison_tables.R")

set_working_dir_to_script()
ctx <- load_us_all_context()

# Keep the Morris baseline frozen for reproducibility.
done_morris <- c(1:5, 7:9, 11:13, 15:17, 33:36, 38:41, 43:46, 51:74)

# Proposed model namespace (AR2 lane). The settings file marks these rows.
done_proposed <- integer(0)
if ("ar2" %in% names(ctx$settings) && "setting_num" %in% names(ctx$settings)) {
  done_proposed <- ctx$settings$setting_num[tolower(ctx$settings$ar2) == "yes"]
  done_proposed <- sort(unique(done_proposed[!is.na(done_proposed)]))
}

all_requested <- sort(unique(c(done_morris, done_proposed)))

probs <- default_probs(include_999 = FALSE)
thresholds <- quantile(ctx$Y, probs = probs, na.rm = TRUE)

score_obj <- compute_us_all_scores(
  setting_ids = all_requested,
  result_path_fn = function(i) sprintf("results/us-all-%d.RData", i),
  Y = ctx$Y,
  cv_lst = ctx$cv_lst,
  probs = probs,
  thresholds = thresholds,
  trans_setting_ids = 2L,
  enforce_contract = TRUE
)
summary_obj <- summarize_us_all_scores(score_obj, baseline_setting = 1L)
print_score_summary(score_obj, label = "us-all-results-proposed")

comparison_full_table <- build_comparison_full_table(
  summary_obj = summary_obj,
  settings = ctx$settings,
  baseline_ids = done_morris,
  proposed_ids = done_proposed
)

comparison_top2 <- build_comparison_top2(
  summary_obj = summary_obj,
  settings = ctx$settings,
  target_quantiles = c(0.95, 0.98, 0.99, 0.995),
  candidate_settings = summary_obj$available_settings,
  metric = "brier"
)

comparison_paired_same_basis <- build_paired_same_basis_table(
  summary_obj = summary_obj,
  settings = ctx$settings,
  baseline_ids = done_morris,
  proposed_ids = done_proposed
)

write.csv(comparison_full_table, "comparison_full_table.csv", row.names = FALSE)
write.csv(comparison_top2, "comparison_top2.csv", row.names = FALSE)
write.csv(comparison_paired_same_basis, "comparison_paired_same_basis.csv", row.names = FALSE)

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
    "done_morris", "done_proposed", "all_requested", "available_settings",
    "skipped_missing_file", "skipped_bad_contract", "skipped_scoring_error",
    "probs", "thresholds", "quant.score", "brier.score",
    "quant.score.mean", "brier.score.mean", "quant.score.se", "brier.score.se",
    "bs.mean.ref.gau", "qs.mean.ref.gau",
    "comparison_full_table", "comparison_top2", "comparison_paired_same_basis",
    "score_obj", "summary_obj"
  ),
  file = "us-all-results-proposed.RData"
)

cat("Outputs written:\n")
cat("- us-all-results-proposed.RData\n")
cat("- comparison_full_table.csv\n")
cat("- comparison_top2.csv\n")
cat("- comparison_paired_same_basis.csv\n")
