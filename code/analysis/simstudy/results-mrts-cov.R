rm(list = ls())

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = FALSE)
  script_dir <- dirname(script_path)
  if (dir.exists(script_dir)) {
    setwd(script_dir)
  }
}

source("./mrts_cov_helpers.R")

load("simdata.RData")
source("../../R/auxfunctions.R")

method_ids <- resolve_mrts_method_ids()

baseline_plan <- get_mrts_baseline_plan(method_ids = method_ids)
baseline_plan$baseline_analysis <- baseline_plan$analysis_id
baseline_plan$baseline_method_id <- baseline_plan$method_id
baseline_plan$baseline_method_key <- baseline_plan$method_key
baseline_plan$baseline_label <- baseline_plan$label

mrts_plan <- get_mrts_analysis_plan(method_ids = method_ids)

analysis_cols <- c(
  "analysis_id", "method_id", "method_key", "family", "label", "mrts_k",
  "output_tag", "baseline_analysis", "baseline_method_id",
  "baseline_method_key", "baseline_label"
)
analysis_plan <- rbind(
  baseline_plan[, analysis_cols, drop = FALSE],
  mrts_plan[, analysis_cols, drop = FALSE]
)
analysis_plan <- analysis_plan[order(analysis_plan$analysis_id), , drop = FALSE]
analysis_plan$analysis_slot <- seq_len(nrow(analysis_plan))

results_dir <- trimws(Sys.getenv("SIMSTUDY_MRTS_RESULTS_DIR", unset = "results"))
setting_ids <- parse_env_indices("SIMSTUDY_MRTS_SETTINGS", seq_len(dim(y)[4]), dim(y)[4])
dataset_ids <- parse_env_indices("SIMSTUDY_MRTS_DATASETS", seq_len(dim(y)[3]), dim(y)[3])
probs <- c(0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995)
n_datasets <- length(dataset_ids)

obs <- c(rep(TRUE, nrow(y) - ntest), rep(FALSE, ntest))

quant.score <- array(NA_real_, dim = c(length(probs), dim(y)[3], nrow(analysis_plan), dim(y)[4]))
brier.score <- array(NA_real_, dim = c(length(probs), dim(y)[3], nrow(analysis_plan), dim(y)[4]))

available_rows <- list()
missing_rows <- list()
row_id <- 1L
missing_id <- 1L

for (setting_id in setting_ids) {
  cat("start setting", setting_id, "\n")
  for (dataset_id in dataset_ids) {
    thresholds <- quantile(y[, , dataset_id, setting_id], probs = probs, na.rm = TRUE)
    validate <- y[!obs, , dataset_id, setting_id]

    for (ii in seq_len(nrow(analysis_plan))) {
      spec <- analysis_plan[ii, , drop = FALSE]
      analysis_slot <- spec$analysis_slot[1]
      result_file <- build_simstudy_result_file(
        results_dir = results_dir,
        setting_id = setting_id,
        method_id = spec$method_id[1],
        dataset_id = dataset_id,
        mrts_k = spec$mrts_k[1]
      )

      if (!file.exists(result_file)) {
        missing_rows[[missing_id]] <- data.frame(
          setting = setting_id,
          dataset = dataset_id,
          analysis_id = spec$analysis_id[1],
          analysis_slot = analysis_slot,
          method_id = spec$method_id[1],
          method_key = spec$method_key[1],
          mrts_k = spec$mrts_k[1],
          result_file = result_file,
          stringsAsFactors = FALSE
        )
        missing_id <- missing_id + 1L
        next
      }

      rm(list = intersect(c("fit.1", "fit", "analysis_spec", "mrts_meta"), ls()))
      load(result_file)

      fit_obj <- NULL
      if (exists("fit.1")) {
        fit_obj <- fit.1
      } else if (exists("fit")) {
        fit_obj <- fit
      }

      if (is.null(fit_obj) || is.null(fit_obj$yp)) {
        next
      }

      pred <- fit_obj$yp
      quant.score[, dataset_id, analysis_slot, setting_id] <- QuantScore(pred, probs, validate)
      brier.score[, dataset_id, analysis_slot, setting_id] <- BrierScore(pred, thresholds, validate)

      available_rows[[row_id]] <- data.frame(
        setting = setting_id,
        dataset = dataset_id,
        analysis_id = spec$analysis_id[1],
        analysis_slot = analysis_slot,
        method_id = spec$method_id[1],
        method_key = spec$method_key[1],
        mrts_k = spec$mrts_k[1],
        result_file = result_file,
        stringsAsFactors = FALSE
      )
      row_id <- row_id + 1L

      rm(list = intersect(c("fit.1", "fit", "analysis_spec", "mrts_meta"), ls()))
    }
  }
}

available_results <- if (length(available_rows) > 0) do.call(rbind, available_rows) else data.frame()
missing_results <- if (length(missing_rows) > 0) do.call(rbind, missing_rows) else data.frame()

quant.score.mean <- apply(quant.score, c(1, 3, 4), mean, na.rm = TRUE)
brier.score.mean <- apply(brier.score, c(1, 3, 4), mean, na.rm = TRUE)
quant.score.se <- apply(quant.score, c(1, 3, 4), sd, na.rm = TRUE) / sqrt(n_datasets)
brier.score.se <- apply(brier.score, c(1, 3, 4), sd, na.rm = TRUE) / sqrt(n_datasets)

paired_rows <- list()
paired_id <- 1L

for (ii in seq_len(nrow(mrts_plan))) {
  spec <- mrts_plan[ii, , drop = FALSE]
  baseline_slot <- analysis_plan$analysis_slot[match(spec$baseline_analysis[1], analysis_plan$analysis_id)]
  proposal_slot <- analysis_plan$analysis_slot[match(spec$analysis_id[1], analysis_plan$analysis_id)]

  for (setting_id in setting_ids) {
    for (qq in seq_along(probs)) {
      base_brier <- brier.score.mean[qq, baseline_slot, setting_id]
      prop_brier <- brier.score.mean[qq, proposal_slot, setting_id]
      base_quant <- quant.score.mean[qq, baseline_slot, setting_id]
      prop_quant <- quant.score.mean[qq, proposal_slot, setting_id]

      paired_rows[[paired_id]] <- data.frame(
        setting = setting_id,
        quantile = probs[qq],
        family = spec$family[1],
        baseline_analysis = spec$baseline_analysis[1],
        baseline_method_id = spec$baseline_method_id[1],
        baseline_method_key = spec$baseline_method_key[1],
        baseline_label = spec$baseline_label[1],
        mrts_analysis = spec$analysis_id[1],
        mrts_method_id = spec$method_id[1],
        mrts_method_key = spec$method_key[1],
        mrts_label = spec$label[1],
        mrts_k = spec$mrts_k[1],
        baseline_brier = base_brier,
        mrts_brier = prop_brier,
        delta_brier = prop_brier - base_brier,
        improve_brier_pct = if (is.finite(base_brier) && base_brier != 0) 100 * (base_brier - prop_brier) / base_brier else NA_real_,
        baseline_quant = base_quant,
        mrts_quant = prop_quant,
        delta_quant = prop_quant - base_quant,
        improve_quant_pct = if (is.finite(base_quant) && base_quant != 0) 100 * (base_quant - prop_quant) / base_quant else NA_real_,
        stringsAsFactors = FALSE
      )
      paired_id <- paired_id + 1L
    }
  }
}

comparison_mrts_cov_paired <- if (length(paired_rows) > 0) do.call(rbind, paired_rows) else data.frame()

summary_rows <- list()
summary_id <- 1L
if (nrow(comparison_mrts_cov_paired) > 0) {
  groups <- split(
    comparison_mrts_cov_paired,
    list(
      comparison_mrts_cov_paired$family,
      comparison_mrts_cov_paired$mrts_k,
      comparison_mrts_cov_paired$quantile
    ),
    drop = TRUE
  )

  for (group_name in names(groups)) {
    group_df <- groups[[group_name]]
    summary_rows[[summary_id]] <- data.frame(
      family = group_df$family[1],
      mrts_k = group_df$mrts_k[1],
      quantile = group_df$quantile[1],
      n_settings = nrow(group_df),
      n_improved_brier = sum(group_df$delta_brier < 0, na.rm = TRUE),
      mean_delta_brier = mean(group_df$delta_brier, na.rm = TRUE),
      mean_improve_brier_pct = mean(group_df$improve_brier_pct, na.rm = TRUE),
      n_improved_quant = sum(group_df$delta_quant < 0, na.rm = TRUE),
      mean_delta_quant = mean(group_df$delta_quant, na.rm = TRUE),
      mean_improve_quant_pct = mean(group_df$improve_quant_pct, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    summary_id <- summary_id + 1L
  }
}

comparison_mrts_cov_summary <- if (length(summary_rows) > 0) do.call(rbind, summary_rows) else data.frame()

write.csv(analysis_plan, "mrts_cov_analysis_plan.csv", row.names = FALSE)
write.csv(available_results, "mrts_cov_available_results.csv", row.names = FALSE)
write.csv(missing_results, "mrts_cov_missing_results.csv", row.names = FALSE)
write.csv(comparison_mrts_cov_paired, "comparison_mrts_cov_paired.csv", row.names = FALSE)
write.csv(comparison_mrts_cov_summary, "comparison_mrts_cov_summary.csv", row.names = FALSE)

save(
  analysis_plan,
  baseline_plan,
  mrts_plan,
  available_results,
  missing_results,
  results_dir,
  setting_ids,
  dataset_ids,
  method_ids,
  probs,
  quant.score,
  brier.score,
  quant.score.mean,
  brier.score.mean,
  quant.score.se,
  brier.score.se,
  comparison_mrts_cov_paired,
  comparison_mrts_cov_summary,
  file = "results-mrts-cov.RData"
)

cat("Outputs written:\n")
cat("- mrts_cov_analysis_plan.csv\n")
cat("- mrts_cov_available_results.csv\n")
cat("- mrts_cov_missing_results.csv\n")
cat("- comparison_mrts_cov_paired.csv\n")
cat("- comparison_mrts_cov_summary.csv\n")
cat("- results-mrts-cov.RData\n")
