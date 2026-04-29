rm(list = ls())

source("refactor_results/00_bootstrap.R")
source("refactor_results/01_score_engine.R")

set_working_dir_to_script()
ctx <- load_us_all_context()
ensure_us_all_output_dirs()

probs <- default_probs(include_999 = FALSE)
threshold_probs <- default_threshold_probs()
thresholds <- quantile(ctx$Y, probs = threshold_probs, na.rm = TRUE)

# Original split from us-all-results-a.R
# - CMAQ files: us-all-<id>.RData
# - no-CMAQ files: us-all-<id>-a.RData
done.cmaq <- c(1, 2, 6, 10)
done.nocmaq <- c(4, 5, 8, 9, 12, 13)
done <- sort(unique(c(done.cmaq, done.nocmaq)))

resolve_result_file <- function(setting_id) {
  pick_existing <- function(candidates) {
    hit <- candidates[file.exists(candidates)]
    if (length(hit) > 0) {
      return(hit[1])
    }
    candidates[1]
  }

  if (setting_id %in% done.cmaq) {
    candidates <- c(
      sprintf("us-all-%d.RData", setting_id),
      file.path("results", sprintf("us-all-%d.RData", setting_id)),
      file.path("results_new", sprintf("us-all-%d.RData", setting_id))
    )
    return(pick_existing(candidates))
  }

  if (setting_id %in% done.nocmaq) {
    candidates <- c(
      sprintf("us-all-%d-a.RData", setting_id),
      file.path("results", sprintf("us-all-%d-a.RData", setting_id)),
      file.path("results_new", sprintf("us-all-%d-a.RData", setting_id))
    )
    return(pick_existing(candidates))
  }

  candidates <- c(
    sprintf("us-all-%d.RData", setting_id),
    file.path("results", sprintf("us-all-%d.RData", setting_id)),
    file.path("results_new", sprintf("us-all-%d.RData", setting_id))
  )
  pick_existing(candidates)
}

score_obj <- compute_us_all_scores(
  setting_ids = done,
  result_path_fn = resolve_result_file,
  Y = ctx$Y,
  cv_lst = ctx$cv_lst,
  probs = probs,
  thresholds = thresholds,
  threshold_probs = threshold_probs,
  trans_setting_ids = integer(0),
  enforce_contract = TRUE
)
summary_obj <- summarize_us_all_scores(score_obj, baseline_setting = 1L)
print_score_summary(score_obj, label = "us-all-results-a (CMAQ vs no-CMAQ)")

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
available_settings <- summary_obj$available_settings

# Legacy compatibility object from us-all-results-a.R
beta.0 <- array(NA_real_, dim = c(5000, length(ctx$cv_lst), summary_obj$max_setting))
savelist <- list(quant.score, brier.score, beta.0, probs, thresholds)
save(savelist, file = us_all_output_path("us-all-results-combined.RData", subdir = "results"))

save(
  list = c(
    "done", "done.cmaq", "done.nocmaq", "available_settings",
    "probs", "threshold_probs", "thresholds", "quant.score", "brier.score",
    "mspe.score", "mape.score",
    "quant.score.mean", "brier.score.mean", "mspe.mean", "mape.mean",
    "quant.score.se", "brier.score.se", "mspe.se", "mape.se",
    "bs.mean.ref.gau", "qs.mean.ref.gau", "mspe.mean.ref.gau", "mape.mean.ref.gau",
    "score_obj", "summary_obj"
  ),
  file = us_all_output_path("us-all-results-a.RData", subdir = "results")
)

q995_idx_brier <- which(abs(threshold_probs - 0.995) < 1e-12)
q995_idx_quant <- which(abs(probs - 0.995) < 1e-12)
summary_table <- data.frame(
  setting = available_settings,
  rel_brier_q995 = if (length(q995_idx_brier) > 0) bs.mean.ref.gau[available_settings, q995_idx_brier] else NA_real_,
  rel_quant_q995 = if (length(q995_idx_quant) > 0) qs.mean.ref.gau[available_settings, q995_idx_quant] else NA_real_,
  mspe_mean = mspe.mean[available_settings],
  mspe_rel_to_gaussian = mspe.mean.ref.gau[available_settings],
  mape_mean = mape.mean[available_settings],
  mape_rel_to_gaussian = mape.mean.ref.gau[available_settings],
  stringsAsFactors = FALSE
)
write.csv(summary_table, us_all_output_path("us-all-results-a-summary.csv", subdir = "tables"), row.names = FALSE)

cat("Outputs written:\n")
cat("- ", us_all_output_path("us-all-results-combined.RData", subdir = "results"), "\n", sep = "")
cat("- ", us_all_output_path("us-all-results-a.RData", subdir = "results"), "\n", sep = "")
cat("- ", us_all_output_path("us-all-results-a-summary.csv", subdir = "tables"), "\n", sep = "")
