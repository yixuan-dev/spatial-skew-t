rm(list = ls())

source("refactor_results/00_bootstrap.R")
source("refactor_results/01_score_engine.R")

set_working_dir_to_script()
ctx <- load_us_all_context()
ensure_us_all_output_dirs()

probs <- default_probs(include_999 = FALSE)
threshold_probs <- default_threshold_probs()
thresholds <- quantile(ctx$Y, probs = threshold_probs, na.rm = TRUE)

result.files <- list.files("results", pattern = "^us-all-[0-9]+\\.RData$", full.names = FALSE)
done <- sort(as.integer(sub("^us-all-([0-9]+)\\.RData$", "\\1", result.files)))
done <- done[!is.na(done) & done >= 1 & done <= 74]

score_obj <- compute_us_all_scores(
  setting_ids = done,
  result_path_fn = function(i) sprintf("results/us-all-%d.RData", i),
  Y = ctx$Y,
  cv_lst = ctx$cv_lst,
  probs = probs,
  thresholds = thresholds,
  threshold_probs = threshold_probs,
  trans_setting_ids = 2L,
  enforce_contract = TRUE
)
summary_obj <- summarize_us_all_scores(score_obj, baseline_setting = 1L)
print_score_summary(score_obj, label = "us-all-results (Morris baseline)")

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

# Legacy compatibility object (beta arrays are not reconstructed here).
savelist <- list(
  quant.score,
  brier.score,
  NA,
  NA,
  probs,
  thresholds
)

save(savelist, file = us_all_output_path("us-all-results-0401.RData", subdir = "results"))
save(
  list = c(
    "done", "available_settings", "probs", "threshold_probs", "thresholds",
    "quant.score", "brier.score", "mspe.score", "mape.score",
    "quant.score.mean", "brier.score.mean", "mspe.mean", "mape.mean",
    "quant.score.se", "brier.score.se", "mspe.se", "mape.se",
    "bs.mean.ref.gau", "qs.mean.ref.gau", "mspe.mean.ref.gau", "mape.mean.ref.gau",
    "score_obj", "summary_obj"
  ),
  file = us_all_output_path("us-all-results.RData", subdir = "results")
)

q99_idx_brier <- which(abs(threshold_probs - 0.99) < 1e-12)
q99_idx_quant <- which(abs(probs - 0.99) < 1e-12)
summary_table <- data.frame(
  setting = available_settings,
  rel_brier_q99 = if (length(q99_idx_brier) > 0) bs.mean.ref.gau[available_settings, q99_idx_brier] else NA_real_,
  rel_quant_q99 = if (length(q99_idx_quant) > 0) qs.mean.ref.gau[available_settings, q99_idx_quant] else NA_real_,
  mspe_mean = mspe.mean[available_settings],
  mspe_rel_to_gaussian = mspe.mean.ref.gau[available_settings],
  mape_mean = mape.mean[available_settings],
  mape_rel_to_gaussian = mape.mean.ref.gau[available_settings],
  stringsAsFactors = FALSE
)
write.csv(summary_table, us_all_output_path("us-all-results-summary.csv", subdir = "tables"), row.names = FALSE)

cat("Outputs written:\n")
cat("- ", us_all_output_path("us-all-results-0401.RData", subdir = "results"), "\n", sep = "")
cat("- ", us_all_output_path("us-all-results.RData", subdir = "results"), "\n", sep = "")
cat("- ", us_all_output_path("us-all-results-summary.csv", subdir = "tables"), "\n", sep = "")
