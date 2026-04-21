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
output_dir <- trimws(Sys.getenv("SIMSTUDY_MRTS_OUTPUT_DIR", unset = "comparison_mrts"))
if (!nzchar(output_dir)) {
  output_dir <- "comparison_mrts"
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(output_dir)) {
  stop(sprintf("Unable to create output directory '%s'", output_dir), call. = FALSE)
}

output_path <- function(filename) {
  file.path(output_dir, filename)
}

setting_ids <- parse_env_indices("SIMSTUDY_MRTS_SETTINGS", seq_len(dim(y)[4]), dim(y)[4])
dataset_ids <- parse_env_indices("SIMSTUDY_MRTS_DATASETS", seq_len(dim(y)[3]), dim(y)[3])
probs <- c(0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995)
n_datasets <- length(dataset_ids)

parse_numeric_env <- function(env_name, default_values) {
  raw <- trimws(Sys.getenv(env_name, unset = ""))
  if (!nzchar(raw)) {
    return(as.numeric(default_values))
  }

  expr <- raw
  if (grepl("^\\(.*\\)$", expr)) {
    expr <- paste0("c", expr)
  }

  values <- tryCatch(
    eval(parse(text = expr), envir = baseenv()),
    error = function(e) {
      stop(sprintf("Invalid %s expression '%s'", env_name, raw), call. = FALSE)
    }
  )

  if (!is.numeric(values) || length(values) == 0 || any(!is.finite(values))) {
    stop(sprintf("%s must be a non-empty numeric vector", env_name), call. = FALSE)
  }

  as.numeric(values)
}

is_selected_quantile <- function(quantiles, selected_quantiles) {
  vapply(
    quantiles,
    function(z) any(abs(z - selected_quantiles) < 1e-12),
    logical(1)
  )
}

classify_quantile_band <- function(quantiles, bulk_range, tail_min) {
  bands <- rep("other", length(quantiles))
  bands[quantiles >= bulk_range[1] & quantiles <= bulk_range[2]] <- "bulk"
  bands[quantiles >= tail_min] <- "tail"
  bands
}

bulk_range <- parse_numeric_env("SIMSTUDY_MRTS_BULK_RANGE", c(0.90, 0.95))
if (length(bulk_range) != 2 || any(bulk_range < 0 | bulk_range > 1)) {
  stop("SIMSTUDY_MRTS_BULK_RANGE must contain two probabilities in [0, 1].", call. = FALSE)
}
bulk_range <- sort(bulk_range)

tail_min <- parse_numeric_env("SIMSTUDY_MRTS_TAIL_MIN", 0.98)
if (length(tail_min) != 1 || tail_min < 0 || tail_min > 1) {
  stop("SIMSTUDY_MRTS_TAIL_MIN must be a single probability in [0, 1].", call. = FALSE)
}

focus_quantiles <- sort(unique(parse_numeric_env(
  "SIMSTUDY_MRTS_FOCUS_QUANTILES",
  c(0.95, 0.98, 0.99)
)))
if (any(focus_quantiles < 0 | focus_quantiles > 1)) {
  stop("SIMSTUDY_MRTS_FOCUS_QUANTILES must contain probabilities in [0, 1].", call. = FALSE)
}

objective <- tolower(trimws(Sys.getenv("SIMSTUDY_MRTS_OBJECTIVE", unset = "balanced")))
if (!objective %in% c("balanced", "extreme-first", "bulk-first")) {
  stop("SIMSTUDY_MRTS_OBJECTIVE must be one of: balanced, extreme-first, bulk-first", call. = FALSE)
}

band_weights <- switch(objective,
  "balanced" = c(bulk = 0.5, tail = 0.5),
  "extreme-first" = c(bulk = 0.25, tail = 0.75),
  "bulk-first" = c(bulk = 0.75, tail = 0.25)
)

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

comparison_mrts_cov_summary_by_band <- data.frame()
comparison_mrts_cov_focus_quantiles <- data.frame()
comparison_mrts_cov_best_k <- data.frame()
comparison_mrts_cov_report_config <- data.frame(
  objective = objective,
  bulk_min = bulk_range[1],
  bulk_max = bulk_range[2],
  tail_min = tail_min,
  focus_quantiles = paste(focus_quantiles, collapse = ","),
  weight_bulk = unname(band_weights["bulk"]),
  weight_tail = unname(band_weights["tail"]),
  stringsAsFactors = FALSE
)

if (nrow(comparison_mrts_cov_paired) > 0) {
  comparison_mrts_cov_paired$quantile_band <- classify_quantile_band(
    comparison_mrts_cov_paired$quantile,
    bulk_range = bulk_range,
    tail_min = tail_min
  )

  band_rows <- list()
  band_id <- 1L
  band_groups <- split(
    comparison_mrts_cov_paired,
    list(
      comparison_mrts_cov_paired$family,
      comparison_mrts_cov_paired$mrts_k,
      comparison_mrts_cov_paired$quantile_band
    ),
    drop = TRUE
  )

  for (group_name in names(band_groups)) {
    group_df <- band_groups[[group_name]]
    n_valid_brier <- sum(is.finite(group_df$delta_brier))
    n_valid_quant <- sum(is.finite(group_df$delta_quant))

    band_rows[[band_id]] <- data.frame(
      family = group_df$family[1],
      mrts_k = group_df$mrts_k[1],
      quantile_band = group_df$quantile_band[1],
      n_pairs = nrow(group_df),
      n_settings = length(unique(group_df$setting)),
      n_improved_brier = sum(group_df$delta_brier < 0, na.rm = TRUE),
      win_rate_brier = if (n_valid_brier > 0) sum(group_df$delta_brier < 0, na.rm = TRUE) / n_valid_brier else NA_real_,
      mean_delta_brier = mean(group_df$delta_brier, na.rm = TRUE),
      mean_improve_brier_pct = mean(group_df$improve_brier_pct, na.rm = TRUE),
      n_improved_quant = sum(group_df$delta_quant < 0, na.rm = TRUE),
      win_rate_quant = if (n_valid_quant > 0) sum(group_df$delta_quant < 0, na.rm = TRUE) / n_valid_quant else NA_real_,
      mean_delta_quant = mean(group_df$delta_quant, na.rm = TRUE),
      mean_improve_quant_pct = mean(group_df$improve_quant_pct, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    band_id <- band_id + 1L
  }

  if (length(band_rows) > 0) {
    comparison_mrts_cov_summary_by_band <- do.call(rbind, band_rows)
    rownames(comparison_mrts_cov_summary_by_band) <- NULL
  }

  focus_rows <- list()
  focus_id <- 1L
  focus_pairs <- comparison_mrts_cov_paired[
    is_selected_quantile(comparison_mrts_cov_paired$quantile, focus_quantiles), ,
    drop = FALSE
  ]

  if (nrow(focus_pairs) > 0) {
    focus_groups <- split(
      focus_pairs,
      list(
        focus_pairs$family,
        focus_pairs$mrts_k,
        focus_pairs$quantile
      ),
      drop = TRUE
    )

    for (group_name in names(focus_groups)) {
      group_df <- focus_groups[[group_name]]
      n_valid_brier <- sum(is.finite(group_df$delta_brier))
      n_valid_quant <- sum(is.finite(group_df$delta_quant))

      focus_rows[[focus_id]] <- data.frame(
        family = group_df$family[1],
        mrts_k = group_df$mrts_k[1],
        quantile = group_df$quantile[1],
        quantile_band = group_df$quantile_band[1],
        n_settings = length(unique(group_df$setting)),
        n_improved_brier = sum(group_df$delta_brier < 0, na.rm = TRUE),
        win_rate_brier = if (n_valid_brier > 0) sum(group_df$delta_brier < 0, na.rm = TRUE) / n_valid_brier else NA_real_,
        mean_improve_brier_pct = mean(group_df$improve_brier_pct, na.rm = TRUE),
        n_improved_quant = sum(group_df$delta_quant < 0, na.rm = TRUE),
        win_rate_quant = if (n_valid_quant > 0) sum(group_df$delta_quant < 0, na.rm = TRUE) / n_valid_quant else NA_real_,
        mean_improve_quant_pct = mean(group_df$improve_quant_pct, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      focus_id <- focus_id + 1L
    }

    comparison_mrts_cov_focus_quantiles <- do.call(rbind, focus_rows)
    comparison_mrts_cov_focus_quantiles <- comparison_mrts_cov_focus_quantiles[
      order(
        comparison_mrts_cov_focus_quantiles$family,
        comparison_mrts_cov_focus_quantiles$mrts_k,
        comparison_mrts_cov_focus_quantiles$quantile
      ), ,
      drop = FALSE
    ]
    rownames(comparison_mrts_cov_focus_quantiles) <- NULL
  }

  if (nrow(comparison_mrts_cov_summary_by_band) > 0) {
    candidate_band_summary <- comparison_mrts_cov_summary_by_band[
      comparison_mrts_cov_summary_by_band$quantile_band %in% c("bulk", "tail"), ,
      drop = FALSE
    ]

    if (nrow(candidate_band_summary) > 0) {
      extract_band_metric <- function(df, metric_name) {
        vals <- c(bulk = NA_real_, tail = NA_real_)
        for (band_name in c("bulk", "tail")) {
          band_df <- df[df$quantile_band == band_name, , drop = FALSE]
          if (nrow(band_df) > 0) {
            vals[band_name] <- band_df[[metric_name]][1]
          }
        }
        vals
      }

      weighted_metric <- function(metric_vec, weights) {
        valid <- is.finite(metric_vec)
        if (!any(valid)) {
          return(NA_real_)
        }

        valid_names <- names(metric_vec)[valid]
        w <- weights[valid_names]
        sum(metric_vec[valid] * w) / sum(w)
      }

      best_rows <- list()
      best_id <- 1L
      keys <- unique(candidate_band_summary[, c("family", "mrts_k"), drop = FALSE])

      for (ii in seq_len(nrow(keys))) {
        key_df <- candidate_band_summary[
          candidate_band_summary$family == keys$family[ii] &
            candidate_band_summary$mrts_k == keys$mrts_k[ii], ,
          drop = FALSE
        ]

        win_rate_brier <- extract_band_metric(key_df, "win_rate_brier")
        improve_brier <- extract_band_metric(key_df, "mean_improve_brier_pct")
        improve_quant <- extract_band_metric(key_df, "mean_improve_quant_pct")

        best_rows[[best_id]] <- data.frame(
          family = keys$family[ii],
          mrts_k = keys$mrts_k[ii],
          objective = objective,
          bulk_win_rate_brier = win_rate_brier["bulk"],
          tail_win_rate_brier = win_rate_brier["tail"],
          bulk_mean_improve_brier_pct = improve_brier["bulk"],
          tail_mean_improve_brier_pct = improve_brier["tail"],
          bulk_mean_improve_quant_pct = improve_quant["bulk"],
          tail_mean_improve_quant_pct = improve_quant["tail"],
          objective_score = weighted_metric(win_rate_brier, band_weights),
          tie_break_brier = weighted_metric(improve_brier, band_weights),
          tie_break_quant = weighted_metric(improve_quant, band_weights),
          stringsAsFactors = FALSE
        )
        best_id <- best_id + 1L
      }

      ranking_base <- do.call(rbind, best_rows)
      ranked_rows <- list()
      ranked_id <- 1L

      for (family_id in unique(ranking_base$family)) {
        family_df <- ranking_base[ranking_base$family == family_id, , drop = FALSE]
        ord <- order(
          -family_df$objective_score,
          -family_df$tie_break_brier,
          -family_df$tie_break_quant,
          family_df$mrts_k,
          na.last = TRUE
        )
        family_df <- family_df[ord, , drop = FALSE]
        family_df$rank_within_family <- seq_len(nrow(family_df))
        family_df$is_selected_best <- family_df$rank_within_family == 1L

        ranked_rows[[ranked_id]] <- family_df
        ranked_id <- ranked_id + 1L
      }

      comparison_mrts_cov_best_k <- do.call(rbind, ranked_rows)
      rownames(comparison_mrts_cov_best_k) <- NULL
    }
  }
}

analysis_plan_csv <- output_path("mrts_cov_analysis_plan.csv")
available_results_csv <- output_path("mrts_cov_available_results.csv")
missing_results_csv <- output_path("mrts_cov_missing_results.csv")
comparison_paired_csv <- output_path("comparison_mrts_cov_paired.csv")
comparison_summary_csv <- output_path("comparison_mrts_cov_summary.csv")
comparison_by_band_csv <- output_path("comparison_mrts_cov_summary_by_band.csv")
comparison_focus_quantiles_csv <- output_path("comparison_mrts_cov_focus_quantiles.csv")
comparison_best_k_csv <- output_path("comparison_mrts_cov_best_k.csv")
comparison_report_config_csv <- output_path("comparison_mrts_cov_report_config.csv")
results_bundle_rdata <- output_path("results-mrts-cov.RData")

write.csv(analysis_plan, analysis_plan_csv, row.names = FALSE)
write.csv(available_results, available_results_csv, row.names = FALSE)
write.csv(missing_results, missing_results_csv, row.names = FALSE)
write.csv(comparison_mrts_cov_paired, comparison_paired_csv, row.names = FALSE)
write.csv(comparison_mrts_cov_summary, comparison_summary_csv, row.names = FALSE)
write.csv(comparison_mrts_cov_summary_by_band, comparison_by_band_csv, row.names = FALSE)
write.csv(comparison_mrts_cov_focus_quantiles, comparison_focus_quantiles_csv, row.names = FALSE)
write.csv(comparison_mrts_cov_best_k, comparison_best_k_csv, row.names = FALSE)
write.csv(comparison_mrts_cov_report_config, comparison_report_config_csv, row.names = FALSE)

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
  output_dir,
  quant.score,
  brier.score,
  quant.score.mean,
  brier.score.mean,
  quant.score.se,
  brier.score.se,
  bulk_range,
  tail_min,
  focus_quantiles,
  objective,
  band_weights,
  comparison_mrts_cov_paired,
  comparison_mrts_cov_summary,
  comparison_mrts_cov_summary_by_band,
  comparison_mrts_cov_focus_quantiles,
  comparison_mrts_cov_best_k,
  comparison_mrts_cov_report_config,
  file = results_bundle_rdata
)

written_outputs <- c(
  analysis_plan_csv,
  available_results_csv,
  missing_results_csv,
  comparison_paired_csv,
  comparison_summary_csv,
  comparison_by_band_csv,
  comparison_focus_quantiles_csv,
  comparison_best_k_csv,
  comparison_report_config_csv,
  results_bundle_rdata
)

cat("Outputs written:\n")
for (path in written_outputs) {
  cat("- ", path, "\n", sep = "")
}
