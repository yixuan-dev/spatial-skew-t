# Comparison-table builders shared by the refactored US-all result scripts.

find_setting_meta <- function(settings, setting_id) {
  settings_local <- settings
  if (!("setting_num" %in% names(settings_local)) && "setting" %in% names(settings_local)) {
    settings_local$setting_num <- safe_as_integer(settings_local$setting)
  }

  if (!("setting_num" %in% names(settings_local))) {
    return(data.frame())
  }

  idx_num <- !is.na(settings_local$setting_num) & settings_local$setting_num == setting_id
  meta <- settings_local[idx_num, , drop = FALSE]
  if (nrow(meta) == 0 && "setting" %in% names(settings_local)) {
    meta <- settings_local[as.character(settings_local$setting) == as.character(setting_id), , drop = FALSE]
  }

  if (nrow(meta) == 0) {
    return(data.frame())
  }

  meta[1, , drop = FALSE]
}

meta_value <- function(meta, column_name, default = NA) {
  if (nrow(meta) == 0 || !(column_name %in% names(meta))) {
    return(default)
  }
  meta[[column_name]][1]
}

base_meta_row <- function(settings, setting_id, baseline_ids = integer(0), proposed_ids = integer(0)) {
  meta <- find_setting_meta(settings, setting_id)

  data.frame(
    setting = setting_id,
    method = meta_value(meta, "method"),
    knots = meta_value(meta, "knots"),
    thresh = meta_value(meta, "thresh"),
    CMAQ = meta_value(meta, "CMAQ"),
    TS = meta_value(meta, "TS"),
    ar2 = meta_value(meta, "ar2"),
    is_baseline = setting_id %in% baseline_ids,
    is_proposed = setting_id %in% proposed_ids,
    stringsAsFactors = FALSE
  )
}

build_comparison_full_table <- function(summary_obj, settings, baseline_ids, proposed_ids) {
  rows <- list()
  row_id <- 1

  for (i in summary_obj$available_settings) {
    meta <- find_setting_meta(settings, i)

    for (q in seq_along(summary_obj$threshold_probs)) {
      rows[[row_id]] <- data.frame(
        setting = i,
        metric = "brier",
        quantile = summary_obj$threshold_probs[q],
        target_type = "threshold_quantile",
        target_level = summary_obj$threshold_probs[q],
        threshold_value = summary_obj$thresholds[q],
        score_mean = summary_obj$brier.score.mean[i, q],
        score_se = summary_obj$brier.score.se[i, q],
        rel_to_gaussian = summary_obj$bs.mean.ref.gau[i, q],
        method = meta_value(meta, "method"),
        knots = meta_value(meta, "knots"),
        thresh = meta_value(meta, "thresh"),
        CMAQ = meta_value(meta, "CMAQ"),
        TS = meta_value(meta, "TS"),
        ar2 = meta_value(meta, "ar2"),
        is_baseline = i %in% baseline_ids,
        is_proposed = i %in% proposed_ids,
        stringsAsFactors = FALSE
      )
      row_id <- row_id + 1
    }

    for (q in seq_along(summary_obj$probs)) {
      rows[[row_id]] <- data.frame(
        setting = i,
        metric = "quantile",
        quantile = summary_obj$probs[q],
        target_type = "score_quantile",
        target_level = summary_obj$probs[q],
        threshold_value = NA_real_,
        score_mean = summary_obj$quant.score.mean[i, q],
        score_se = summary_obj$quant.score.se[i, q],
        rel_to_gaussian = summary_obj$qs.mean.ref.gau[i, q],
        method = meta_value(meta, "method"),
        knots = meta_value(meta, "knots"),
        thresh = meta_value(meta, "thresh"),
        CMAQ = meta_value(meta, "CMAQ"),
        TS = meta_value(meta, "TS"),
        ar2 = meta_value(meta, "ar2"),
        is_baseline = i %in% baseline_ids,
        is_proposed = i %in% proposed_ids,
        stringsAsFactors = FALSE
      )
      row_id <- row_id + 1
    }
  }

  if (length(rows) == 0) {
    return(data.frame())
  }

  do.call(rbind, rows)
}

build_comparison_top2 <- function(
    summary_obj,
    settings,
    target_levels = c(0.95, 0.98, 0.99, 0.995),
    candidate_settings = NULL,
    metric = c("brier", "quantile")
) {
  metric <- match.arg(metric)
  score_mat <- if (metric == "brier") summary_obj$bs.mean.ref.gau else summary_obj$qs.mean.ref.gau
  level_grid <- if (metric == "brier") summary_obj$threshold_probs else summary_obj$probs
  target_type <- if (metric == "brier") "threshold_quantile" else "score_quantile"

  if (is.null(candidate_settings)) {
    candidate_settings <- summary_obj$available_settings
  }

  rows <- list()
  row_id <- 1

  for (q in target_levels) {
    q_idx <- which(abs(level_grid - q) < 1e-12)
    if (length(q_idx) == 0) {
      next
    }

    rel_vec <- score_mat[candidate_settings, q_idx]
    ord <- order(rel_vec, decreasing = FALSE, na.last = NA)
    if (length(ord) == 0) {
      next
    }

    winners <- candidate_settings[ord[seq_len(min(2, length(ord)))]]
    for (rank_idx in seq_along(winners)) {
      w <- winners[rank_idx]
      meta <- find_setting_meta(settings, w)

      rows[[row_id]] <- data.frame(
        metric = metric,
        quantile = q,
        target_type = target_type,
        target_level = q,
        threshold_value = if (metric == "brier") summary_obj$thresholds[q_idx] else NA_real_,
        rank = rank_idx,
        setting = w,
        rel_score_to_gaussian = score_mat[w, q_idx],
        method = meta_value(meta, "method"),
        knots = meta_value(meta, "knots"),
        thresh = meta_value(meta, "thresh"),
        CMAQ = meta_value(meta, "CMAQ"),
        TS = meta_value(meta, "TS"),
        ar2 = meta_value(meta, "ar2"),
        stringsAsFactors = FALSE
      )
      row_id <- row_id + 1
    }
  }

  if (length(rows) == 0) {
    return(data.frame())
  }

  do.call(rbind, rows)
}

make_top2_rows <- function(
    settings,
    candidate_settings,
    ranking_values,
    metric_family,
    metric,
    ranking_group,
    ranking_basis,
    score_direction,
    score_values = NULL,
    score_ses = NULL,
    rel_to_gaussian = NULL,
    gap_to_target = NULL,
    target_type = "scalar",
    target_level = NA_real_,
    threshold_value = NA_real_,
    band_type = NA_character_,
    target_value = NA_real_,
    prediction_cutoff = NA_real_,
    decreasing = FALSE,
    n_top = 2L
) {
  candidate_settings <- sort(unique(candidate_settings[!is.na(candidate_settings)]))
  candidate_settings <- candidate_settings[candidate_settings <= length(ranking_values)]
  candidate_settings <- candidate_settings[is.finite(ranking_values[candidate_settings])]

  if (length(candidate_settings) == 0) {
    return(data.frame())
  }

  ord <- order(ranking_values[candidate_settings], decreasing = decreasing, na.last = NA)
  if (length(ord) == 0) {
    return(data.frame())
  }

  winners <- candidate_settings[ord[seq_len(min(n_top, length(ord)))]]
  rows <- vector("list", length(winners))

  for (rank_idx in seq_along(winners)) {
    w <- winners[rank_idx]
    meta <- find_setting_meta(settings, w)

    rows[[rank_idx]] <- data.frame(
      metric_family = metric_family,
      metric = metric,
      ranking_group = ranking_group,
      quantile = if (grepl("quantile", target_type, fixed = TRUE)) target_level else NA_real_,
      target_type = target_type,
      target_level = target_level,
      threshold_value = threshold_value,
      band_type = band_type,
      target_value = target_value,
      prediction_cutoff = prediction_cutoff,
      ranking_basis = ranking_basis,
      rank = rank_idx,
      setting = w,
      rank_value = ranking_values[w],
      score_value = if (is.null(score_values)) ranking_values[w] else score_values[w],
      score_se = if (is.null(score_ses)) NA_real_ else score_ses[w],
      rel_score_to_gaussian = if (is.null(rel_to_gaussian)) NA_real_ else rel_to_gaussian[w],
      rel_to_gaussian = if (is.null(rel_to_gaussian)) NA_real_ else rel_to_gaussian[w],
      gap_to_target = if (is.null(gap_to_target)) NA_real_ else gap_to_target[w],
      score_direction = score_direction,
      method = meta_value(meta, "method"),
      knots = meta_value(meta, "knots"),
      thresh = meta_value(meta, "thresh"),
      CMAQ = meta_value(meta, "CMAQ"),
      TS = meta_value(meta, "TS"),
      ar2 = meta_value(meta, "ar2"),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

build_comparison_top2_all_metrics <- function(
    summary_obj,
    settings,
    candidate_settings = NULL,
    n_top = 2L
) {
  if (is.null(candidate_settings)) {
    candidate_settings <- summary_obj$available_settings
  }
  candidate_settings <- intersect(candidate_settings, summary_obj$available_settings)

  blocks <- list()
  block_id <- 1L

  for (q_idx in seq_along(summary_obj$threshold_probs)) {
    q <- summary_obj$threshold_probs[q_idx]
    blocks[[block_id]] <- make_top2_rows(
      settings = settings,
      candidate_settings = candidate_settings,
      ranking_values = summary_obj$bs.mean.ref.gau[, q_idx],
      score_values = summary_obj$brier.score.mean[, q_idx],
      score_ses = summary_obj$brier.score.se[, q_idx],
      rel_to_gaussian = summary_obj$bs.mean.ref.gau[, q_idx],
      metric_family = "score",
      metric = "brier",
      ranking_group = sprintf("brier__q%s", format(q, trim = TRUE, scientific = FALSE)),
      ranking_basis = "rel_to_gaussian",
      score_direction = "lower_better",
      target_type = "threshold_quantile",
      target_level = q,
      threshold_value = summary_obj$thresholds[q_idx],
      n_top = n_top
    )
    block_id <- block_id + 1L
  }

  for (q_idx in seq_along(summary_obj$probs)) {
    q <- summary_obj$probs[q_idx]
    blocks[[block_id]] <- make_top2_rows(
      settings = settings,
      candidate_settings = candidate_settings,
      ranking_values = summary_obj$qs.mean.ref.gau[, q_idx],
      score_values = summary_obj$quant.score.mean[, q_idx],
      score_ses = summary_obj$quant.score.se[, q_idx],
      rel_to_gaussian = summary_obj$qs.mean.ref.gau[, q_idx],
      metric_family = "score",
      metric = "quantile",
      ranking_group = sprintf("quantile__q%s", format(q, trim = TRUE, scientific = FALSE)),
      ranking_basis = "rel_to_gaussian",
      score_direction = "lower_better",
      target_type = "score_quantile",
      target_level = q,
      n_top = n_top
    )
    block_id <- block_id + 1L
  }

  if (!is.null(summary_obj$crps.mean)) {
    blocks[[block_id]] <- make_top2_rows(
      settings = settings,
      candidate_settings = candidate_settings,
      ranking_values = summary_obj$crps.mean.ref.gau,
      score_values = summary_obj$crps.mean,
      score_ses = summary_obj$crps.se,
      rel_to_gaussian = summary_obj$crps.mean.ref.gau,
      metric_family = "scalar",
      metric = "crps",
      ranking_group = "crps",
      ranking_basis = "rel_to_gaussian",
      score_direction = "lower_better",
      n_top = n_top
    )
    block_id <- block_id + 1L
  }

  if (!is.null(summary_obj$brier.split.score.mean)) {
    for (target_idx in seq_along(summary_obj$brier.split.target_probs)) {
      q <- summary_obj$brier.split.target_probs[target_idx]
      threshold_value <- summary_obj$brier.split.target_thresholds[target_idx]
      for (band_idx in seq_along(summary_obj$brier.split.band_names)) {
        band <- summary_obj$brier.split.band_names[band_idx]
        blocks[[block_id]] <- make_top2_rows(
          settings = settings,
          candidate_settings = candidate_settings,
          ranking_values = summary_obj$brier.split.rel.ref.gau[, band_idx, target_idx],
          score_values = summary_obj$brier.split.score.mean[, band_idx, target_idx],
          score_ses = summary_obj$brier.split.score.se[, band_idx, target_idx],
          rel_to_gaussian = summary_obj$brier.split.rel.ref.gau[, band_idx, target_idx],
          metric_family = "split_score",
          metric = "brier_split",
          ranking_group = sprintf("brier_split__%s__q%s", band, format(q, trim = TRUE, scientific = FALSE)),
          ranking_basis = "rel_to_gaussian",
          score_direction = "lower_better",
          target_type = "threshold_quantile",
          target_level = q,
          threshold_value = threshold_value,
          band_type = band,
          n_top = n_top
        )
        block_id <- block_id + 1L
      }
    }
  }

  if (!is.null(summary_obj$classification.metric.mean)) {
    for (metric_idx in seq_along(summary_obj$classification.metric_names)) {
      metric_name <- summary_obj$classification.metric_names[metric_idx]
      for (target_idx in seq_along(summary_obj$classification.target_probs)) {
        q <- summary_obj$classification.target_probs[target_idx]
        threshold_value <- summary_obj$classification.target_thresholds[target_idx]
        blocks[[block_id]] <- make_top2_rows(
          settings = settings,
          candidate_settings = candidate_settings,
          ranking_values = summary_obj$classification.metric.mean[, metric_idx, target_idx],
          score_values = summary_obj$classification.metric.mean[, metric_idx, target_idx],
          score_ses = summary_obj$classification.metric.se[, metric_idx, target_idx],
          rel_to_gaussian = summary_obj$classification.metric.rel.ref.gau[, metric_idx, target_idx],
          gap_to_target = summary_obj$classification.metric.delta.ref.gau[, metric_idx, target_idx],
          metric_family = "classification",
          metric = metric_name,
          ranking_group = sprintf(
            "classification__%s__q%s__cutoff%s",
            metric_name,
            format(q, trim = TRUE, scientific = FALSE),
            format(summary_obj$classification.probability_cutoff, trim = TRUE, scientific = FALSE)
          ),
          ranking_basis = "raw_score",
          score_direction = "higher_better",
          target_type = "threshold_quantile",
          target_level = q,
          threshold_value = threshold_value,
          prediction_cutoff = summary_obj$classification.probability_cutoff,
          decreasing = TRUE,
          n_top = n_top
        )
        block_id <- block_id + 1L
      }
    }
  }

  if (!is.null(summary_obj$coverage.mean)) {
    for (k in seq_along(summary_obj$uncertainty_levels)) {
      level <- summary_obj$uncertainty_levels[k]
      blocks[[block_id]] <- make_top2_rows(
        settings = settings,
        candidate_settings = candidate_settings,
        ranking_values = abs(summary_obj$coverage.gap[, k]),
        score_values = summary_obj$coverage.mean[, k],
        score_ses = summary_obj$coverage.se[, k],
        gap_to_target = summary_obj$coverage.gap[, k],
        metric_family = "uncertainty",
        metric = "coverage",
        ranking_group = sprintf("coverage__level%s", format(level, trim = TRUE, scientific = FALSE)),
        ranking_basis = "abs_gap_to_target",
        score_direction = "closer_to_target_better",
        target_type = "coverage_level",
        target_level = level,
        target_value = level,
        n_top = n_top
      )
      block_id <- block_id + 1L
    }
  }

  if (!is.null(summary_obj$pit.mean)) {
    blocks[[block_id]] <- make_top2_rows(
      settings = settings,
      candidate_settings = candidate_settings,
      ranking_values = abs(summary_obj$pit.mean - summary_obj$pit.expected_mean),
      score_values = summary_obj$pit.mean,
      score_ses = summary_obj$pit.mean.se,
      gap_to_target = summary_obj$pit.mean - summary_obj$pit.expected_mean,
      metric_family = "uncertainty",
      metric = "pit_mean",
      ranking_group = "pit_mean",
      ranking_basis = "abs_gap_to_target",
      score_direction = "closer_to_target_better",
      target_type = "pit_summary",
      target_value = summary_obj$pit.expected_mean,
      n_top = n_top
    )
    block_id <- block_id + 1L
  }

  if (!is.null(summary_obj$pit.variance)) {
    blocks[[block_id]] <- make_top2_rows(
      settings = settings,
      candidate_settings = candidate_settings,
      ranking_values = abs(summary_obj$pit.variance - summary_obj$pit.expected_variance),
      score_values = summary_obj$pit.variance,
      score_ses = summary_obj$pit.variance.se,
      gap_to_target = summary_obj$pit.variance - summary_obj$pit.expected_variance,
      metric_family = "uncertainty",
      metric = "pit_variance",
      ranking_group = "pit_variance",
      ranking_basis = "abs_gap_to_target",
      score_direction = "closer_to_target_better",
      target_type = "pit_summary",
      target_value = summary_obj$pit.expected_variance,
      n_top = n_top
    )
    block_id <- block_id + 1L
  }

  if (!is.null(summary_obj$pit.ks)) {
    blocks[[block_id]] <- make_top2_rows(
      settings = settings,
      candidate_settings = candidate_settings,
      ranking_values = summary_obj$pit.ks,
      score_values = summary_obj$pit.ks,
      score_ses = summary_obj$pit.ks.se,
      metric_family = "uncertainty",
      metric = "pit_ks",
      ranking_group = "pit_ks",
      ranking_basis = "raw_score",
      score_direction = "lower_better",
      target_type = "pit_summary",
      n_top = n_top
    )
    block_id <- block_id + 1L
  }

  if (!is.null(summary_obj$pit.mae)) {
    blocks[[block_id]] <- make_top2_rows(
      settings = settings,
      candidate_settings = candidate_settings,
      ranking_values = summary_obj$pit.mae,
      score_values = summary_obj$pit.mae,
      score_ses = summary_obj$pit.mae.se,
      metric_family = "uncertainty",
      metric = "pit_uniformity_mae",
      ranking_group = "pit_uniformity_mae",
      ranking_basis = "raw_score",
      score_direction = "lower_better",
      target_type = "pit_summary",
      n_top = n_top
    )
    block_id <- block_id + 1L
  }

  if (!is.null(summary_obj$pit.rmse)) {
    blocks[[block_id]] <- make_top2_rows(
      settings = settings,
      candidate_settings = candidate_settings,
      ranking_values = summary_obj$pit.rmse,
      score_values = summary_obj$pit.rmse,
      score_ses = summary_obj$pit.rmse.se,
      metric_family = "uncertainty",
      metric = "pit_uniformity_rmse",
      ranking_group = "pit_uniformity_rmse",
      ranking_basis = "raw_score",
      score_direction = "lower_better",
      target_type = "pit_summary",
      n_top = n_top
    )
    block_id <- block_id + 1L
  }

  blocks <- blocks[vapply(blocks, nrow, integer(1)) > 0]
  if (length(blocks) == 0) {
    return(data.frame())
  }

  do.call(rbind, blocks)
}

sanitize_excel_sheet_name <- function(x, fallback = "sheet") {
  out <- gsub("[\\\\/:?*\\[\\]]", "_", as.character(x))
  out <- trimws(out)
  out[!nzchar(out)] <- fallback
  out <- substr(out, 1, 31)

  seen <- character(0)
  for (i in seq_along(out)) {
    candidate <- out[i]
    if (!(candidate %in% seen)) {
      seen <- c(seen, candidate)
      out[i] <- candidate
      next
    }

    suffix_id <- 2L
    repeat {
      suffix <- paste0("_", suffix_id)
      stem <- substr(candidate, 1, max(1L, 31L - nchar(suffix)))
      proposal <- paste0(stem, suffix)
      if (!(proposal %in% seen)) {
        seen <- c(seen, proposal)
        out[i] <- proposal
        break
      }
      suffix_id <- suffix_id + 1L
    }
  }

  out
}

write_comparison_top2_workbook <- function(comparison_top2, output_path) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    return(invisible(FALSE))
  }

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "all_metrics")
  openxlsx::writeData(wb, "all_metrics", comparison_top2)

  metric_names <- unique(as.character(comparison_top2$metric))
  metric_sheet_names <- sanitize_excel_sheet_name(metric_names, fallback = "metric")

  for (i in seq_along(metric_names)) {
    metric_name <- metric_names[i]
    sheet_name <- metric_sheet_names[i]
    metric_df <- comparison_top2[as.character(comparison_top2$metric) == metric_name, , drop = FALSE]
    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeData(wb, sheet_name, metric_df)
  }

  openxlsx::saveWorkbook(wb, output_path, overwrite = TRUE)
  invisible(TRUE)
}

build_paired_same_basis_table <- function(summary_obj, settings, baseline_ids, proposed_ids) {
  join_cols <- c("method", "knots", "thresh", "CMAQ", "TS")
  if (!all(join_cols %in% names(settings))) {
    return(data.frame())
  }

  settings_local <- settings
  if (!("setting_num" %in% names(settings_local)) && "setting" %in% names(settings_local)) {
    settings_local$setting_num <- safe_as_integer(settings_local$setting)
  }
  if (!("setting_num" %in% names(settings_local))) {
    return(data.frame())
  }

  base_meta <- settings_local[settings_local$setting_num %in% baseline_ids, , drop = FALSE]
  prop_meta <- settings_local[settings_local$setting_num %in% proposed_ids, , drop = FALSE]

  rows <- list()
  row_id <- 1

  for (k in seq_len(nrow(prop_meta))) {
    p <- prop_meta[k, , drop = FALSE]

    same_basis <- rep(TRUE, nrow(base_meta))
    for (nm in join_cols) {
      lhs <- as.character(base_meta[[nm]])
      rhs <- as.character(p[[nm]][1])
      both_na <- is.na(lhs) & is.na(rhs)
      same_basis <- same_basis & (lhs == rhs | both_na)
    }

    matches <- base_meta[same_basis, , drop = FALSE]
    if (nrow(matches) == 0) {
      next
    }

    b_setting <- matches$setting_num[1]
    p_setting <- p$setting_num[1]

    if (!(b_setting %in% summary_obj$available_settings && p_setting %in% summary_obj$available_settings)) {
      next
    }

    for (q in seq_along(summary_obj$threshold_probs)) {
      rows[[row_id]] <- data.frame(
        metric = "brier",
        quantile = summary_obj$threshold_probs[q],
        target_type = "threshold_quantile",
        target_level = summary_obj$threshold_probs[q],
        threshold_value = summary_obj$thresholds[q],
        baseline_setting = b_setting,
        proposed_setting = p_setting,
        baseline_rel_brier = summary_obj$bs.mean.ref.gau[b_setting, q],
        proposed_rel_brier = summary_obj$bs.mean.ref.gau[p_setting, q],
        delta_rel_brier = summary_obj$bs.mean.ref.gau[p_setting, q] - summary_obj$bs.mean.ref.gau[b_setting, q],
        baseline_rel_score = summary_obj$bs.mean.ref.gau[b_setting, q],
        proposed_rel_score = summary_obj$bs.mean.ref.gau[p_setting, q],
        delta_rel_score = summary_obj$bs.mean.ref.gau[p_setting, q] - summary_obj$bs.mean.ref.gau[b_setting, q],
        baseline_rel_quant = NA_real_,
        proposed_rel_quant = NA_real_,
        delta_rel_quant = NA_real_,
        method = p$method[1],
        knots = p$knots[1],
        thresh = p$thresh[1],
        CMAQ = p$CMAQ[1],
        TS = p$TS[1],
        stringsAsFactors = FALSE
      )
      row_id <- row_id + 1
    }

    for (q in seq_along(summary_obj$probs)) {
      rows[[row_id]] <- data.frame(
        metric = "quantile",
        quantile = summary_obj$probs[q],
        target_type = "score_quantile",
        target_level = summary_obj$probs[q],
        threshold_value = NA_real_,
        baseline_setting = b_setting,
        proposed_setting = p_setting,
        baseline_rel_brier = NA_real_,
        proposed_rel_brier = NA_real_,
        delta_rel_brier = NA_real_,
        baseline_rel_score = summary_obj$qs.mean.ref.gau[b_setting, q],
        proposed_rel_score = summary_obj$qs.mean.ref.gau[p_setting, q],
        delta_rel_score = summary_obj$qs.mean.ref.gau[p_setting, q] - summary_obj$qs.mean.ref.gau[b_setting, q],
        baseline_rel_quant = summary_obj$qs.mean.ref.gau[b_setting, q],
        proposed_rel_quant = summary_obj$qs.mean.ref.gau[p_setting, q],
        delta_rel_quant = summary_obj$qs.mean.ref.gau[p_setting, q] - summary_obj$qs.mean.ref.gau[b_setting, q],
        method = p$method[1],
        knots = p$knots[1],
        thresh = p$thresh[1],
        CMAQ = p$CMAQ[1],
        TS = p$TS[1],
        stringsAsFactors = FALSE
      )
      row_id <- row_id + 1
    }
  }

  if (length(rows) == 0) {
    return(data.frame())
  }

  do.call(rbind, rows)
}

build_comparison_brier_split_table <- function(summary_obj, settings, baseline_ids, proposed_ids) {
  if (
    is.null(summary_obj$brier.split.score.mean) ||
    length(summary_obj$brier.split.target_probs) == 0 ||
    length(summary_obj$brier.split.band_names) == 0
  ) {
    return(data.frame())
  }

  rows <- list()
  row_id <- 1

  for (i in summary_obj$available_settings) {
    base_row <- base_meta_row(settings, i, baseline_ids = baseline_ids, proposed_ids = proposed_ids)

    for (target_idx in seq_along(summary_obj$brier.split.target_probs)) {
      event_quantile <- summary_obj$brier.split.target_probs[target_idx]
      event_threshold_value <- summary_obj$brier.split.target_thresholds[target_idx]

      for (band_idx in seq_along(summary_obj$brier.split.band_names)) {
        band_type <- summary_obj$brier.split.band_names[band_idx]

        rows[[row_id]] <- cbind(
          data.frame(
            band_type = band_type,
            event_quantile = event_quantile,
            event_threshold_value = event_threshold_value,
            score_mean = summary_obj$brier.split.score.mean[i, band_idx, target_idx],
            score_se = summary_obj$brier.split.score.se[i, band_idx, target_idx],
            rel_to_gaussian = summary_obj$brier.split.rel.ref.gau[i, band_idx, target_idx],
            n_obs_total = summary_obj$brier.split.n_obs.total[i, band_idx, target_idx],
            n_obs_mean_per_fold = summary_obj$brier.split.n_obs.mean[i, band_idx, target_idx],
            obs_share = summary_obj$brier.split.obs.share[i, band_idx, target_idx],
            score_direction = "lower_better",
            stringsAsFactors = FALSE
          ),
          base_row
        )
        row_id <- row_id + 1
      }
    }
  }

  if (length(rows) == 0) {
    return(data.frame())
  }

  do.call(rbind, rows)
}

build_comparison_scalar_metrics_table <- function(summary_obj, settings, baseline_ids, proposed_ids) {
  rows <- list()
  row_id <- 1

  for (i in summary_obj$available_settings) {
    base_row <- base_meta_row(settings, i, baseline_ids = baseline_ids, proposed_ids = proposed_ids)

    rows[[row_id]] <- cbind(
      data.frame(
        metric_family = "scalar",
        crps_mean = summary_obj$crps.mean[i],
        crps_se = summary_obj$crps.se[i],
        crps_rel_to_gaussian = summary_obj$crps.mean.ref.gau[i],
        crps_direction = "lower_better",
        crps_draws_used = summary_obj$summary.draws.mean[i],
        n_obs_total = summary_obj$summary.n_obs.total[i],
        loo_elpd = NA_real_,
        loo_elpd_se = NA_real_,
        delta_loo_elpd_to_gaussian = NA_real_,
        loo_direction = "higher_better",
        waic = NA_real_,
        waic_se = NA_real_,
        delta_waic_to_gaussian = NA_real_,
        waic_direction = "lower_better",
        predictive_density_status = summary_obj$predictive_density_status[i],
        predictive_density_note = summary_obj$predictive_density_note[i],
        stringsAsFactors = FALSE
      ),
      base_row
    )
    row_id <- row_id + 1
  }

  if (length(rows) == 0) {
    return(data.frame())
  }

  do.call(rbind, rows)
}

build_comparison_classification_metrics_table <- function(summary_obj, settings, baseline_ids, proposed_ids) {
  if (
    is.null(summary_obj$classification.metric.mean) ||
    length(summary_obj$classification.target_probs) == 0 ||
    length(summary_obj$classification.metric_names) == 0 ||
    length(summary_obj$classification.count_names) == 0
  ) {
    return(data.frame())
  }

  rows <- list()
  row_id <- 1

  for (i in summary_obj$available_settings) {
    base_row <- base_meta_row(settings, i, baseline_ids = baseline_ids, proposed_ids = proposed_ids)

    for (target_idx in seq_along(summary_obj$classification.target_probs)) {
      event_quantile <- summary_obj$classification.target_probs[target_idx]
      event_threshold_value <- summary_obj$classification.target_thresholds[target_idx]

      rows[[row_id]] <- cbind(
        data.frame(
          event_quantile = event_quantile,
          event_threshold_value = event_threshold_value,
          probability_cutoff = summary_obj$classification.probability_cutoff,
          classification_draws_used = summary_obj$classification.draws.mean[i],
          n_obs_total = summary_obj$classification.n_obs.total[i, target_idx],
          n_obs_mean_per_fold = summary_obj$classification.n_obs.mean[i, target_idx],
          actual_positive_share = summary_obj$classification.actual_positive_share[i, target_idx],
          predicted_positive_share = summary_obj$classification.predicted_positive_share[i, target_idx],
          tp_total = summary_obj$classification.count.total[i, "tp", target_idx],
          tn_total = summary_obj$classification.count.total[i, "tn", target_idx],
          fp_total = summary_obj$classification.count.total[i, "fp", target_idx],
          fn_total = summary_obj$classification.count.total[i, "fn", target_idx],
          tp_mean_per_fold = summary_obj$classification.count.mean[i, "tp", target_idx],
          tn_mean_per_fold = summary_obj$classification.count.mean[i, "tn", target_idx],
          fp_mean_per_fold = summary_obj$classification.count.mean[i, "fp", target_idx],
          fn_mean_per_fold = summary_obj$classification.count.mean[i, "fn", target_idx],
          accuracy_mean = summary_obj$classification.metric.mean[i, "accuracy", target_idx],
          accuracy_se = summary_obj$classification.metric.se[i, "accuracy", target_idx],
          accuracy_rel_to_gaussian = summary_obj$classification.metric.rel.ref.gau[i, "accuracy", target_idx],
          accuracy_delta_to_gaussian = summary_obj$classification.metric.delta.ref.gau[i, "accuracy", target_idx],
          precision_mean = summary_obj$classification.metric.mean[i, "precision", target_idx],
          precision_se = summary_obj$classification.metric.se[i, "precision", target_idx],
          precision_rel_to_gaussian = summary_obj$classification.metric.rel.ref.gau[i, "precision", target_idx],
          precision_delta_to_gaussian = summary_obj$classification.metric.delta.ref.gau[i, "precision", target_idx],
          recall_mean = summary_obj$classification.metric.mean[i, "recall", target_idx],
          recall_se = summary_obj$classification.metric.se[i, "recall", target_idx],
          recall_rel_to_gaussian = summary_obj$classification.metric.rel.ref.gau[i, "recall", target_idx],
          recall_delta_to_gaussian = summary_obj$classification.metric.delta.ref.gau[i, "recall", target_idx],
          specificity_mean = summary_obj$classification.metric.mean[i, "specificity", target_idx],
          specificity_se = summary_obj$classification.metric.se[i, "specificity", target_idx],
          specificity_rel_to_gaussian = summary_obj$classification.metric.rel.ref.gau[i, "specificity", target_idx],
          specificity_delta_to_gaussian = summary_obj$classification.metric.delta.ref.gau[i, "specificity", target_idx],
          f1_mean = summary_obj$classification.metric.mean[i, "f1", target_idx],
          f1_se = summary_obj$classification.metric.se[i, "f1", target_idx],
          f1_rel_to_gaussian = summary_obj$classification.metric.rel.ref.gau[i, "f1", target_idx],
          f1_delta_to_gaussian = summary_obj$classification.metric.delta.ref.gau[i, "f1", target_idx],
          accuracy_direction = "higher_better",
          precision_direction = "higher_better",
          recall_direction = "higher_better",
          specificity_direction = "higher_better",
          f1_direction = "higher_better",
          count_direction = "tp_tn_higher_better; fp_fn_lower_better",
          stringsAsFactors = FALSE
        ),
        base_row
      )
      row_id <- row_id + 1
    }
  }

  if (length(rows) == 0) {
    return(data.frame())
  }

  do.call(rbind, rows)
}

build_comparison_uncertainty_summary_table <- function(summary_obj, settings, baseline_ids, proposed_ids) {
  rows <- list()
  row_id <- 1

  for (i in summary_obj$available_settings) {
    base_row <- base_meta_row(settings, i, baseline_ids = baseline_ids, proposed_ids = proposed_ids)
    row <- list(
      uncertainty_family = "coverage_pit",
      summary_draws_used = summary_obj$summary.draws.mean[i],
      n_obs_total = summary_obj$summary.n_obs.total[i],
      pit_mean = summary_obj$pit.mean[i],
      pit_mean_se = summary_obj$pit.mean.se[i],
      pit_mean_target = summary_obj$pit.expected_mean,
      pit_mean_gap = if (is.finite(summary_obj$pit.mean[i])) summary_obj$pit.mean[i] - summary_obj$pit.expected_mean else NA_real_,
      pit_variance = summary_obj$pit.variance[i],
      pit_variance_se = summary_obj$pit.variance.se[i],
      pit_variance_target = summary_obj$pit.expected_variance,
      pit_variance_gap = if (is.finite(summary_obj$pit.variance[i])) summary_obj$pit.variance[i] - summary_obj$pit.expected_variance else NA_real_,
      pit_ks_stat = summary_obj$pit.ks[i],
      pit_ks_se = summary_obj$pit.ks.se[i],
      pit_uniformity_mae = summary_obj$pit.mae[i],
      pit_uniformity_mae_se = summary_obj$pit.mae.se[i],
      pit_uniformity_rmse = summary_obj$pit.rmse[i],
      pit_uniformity_rmse_se = summary_obj$pit.rmse.se[i],
      coverage_direction = "closer_to_target_better",
      pit_direction = "pit_mean_to_0.5; pit_variance_to_1_over_12; ks_mae_rmse_lower_better"
    )

    for (k in seq_along(summary_obj$uncertainty_levels)) {
      level <- summary_obj$uncertainty_levels[k]
      label <- sprintf("%d", as.integer(round(level * 100)))
      row[[paste0("coverage_", label, "_target")]] <- level
      row[[paste0("coverage_", label, "_mean")]] <- summary_obj$coverage.mean[i, k]
      row[[paste0("coverage_", label, "_se")]] <- summary_obj$coverage.se[i, k]
      row[[paste0("coverage_", label, "_gap")]] <- summary_obj$coverage.gap[i, k]
    }

    rows[[row_id]] <- cbind(as.data.frame(row, stringsAsFactors = FALSE), base_row)
    row_id <- row_id + 1
  }

  if (length(rows) == 0) {
    return(data.frame())
  }

  do.call(rbind, rows)
}

build_comparison_calibration_bins_table <- function(summary_obj, settings, baseline_ids, proposed_ids) {
  rows <- list()
  row_id <- 1
  pit_breaks <- summary_obj$pit_breaks
  expected_share <- summary_obj$pit.expected_share

  for (i in summary_obj$available_settings) {
    base_row <- base_meta_row(settings, i, baseline_ids = baseline_ids, proposed_ids = proposed_ids)

    for (b in seq_along(expected_share)) {
      rows[[row_id]] <- cbind(
        data.frame(
          calibration_type = "pit_histogram",
          bin_id = b,
          bin_left = pit_breaks[b],
          bin_right = pit_breaks[b + 1],
          bin_mid = (pit_breaks[b] + pit_breaks[b + 1]) / 2,
          expected_share = expected_share[b],
          observed_share_mean = summary_obj$pit.bin.share.mean[i, b],
          observed_share_se = summary_obj$pit.bin.share.se[i, b],
          share_gap = if (is.finite(summary_obj$pit.bin.share.mean[i, b])) summary_obj$pit.bin.share.mean[i, b] - expected_share[b] else NA_real_,
          count_total = summary_obj$pit.bin.count.total[i, b],
          count_mean_per_fold = summary_obj$pit.bin.count.mean[i, b],
          summary_draws_used = summary_obj$summary.draws.mean[i],
          n_obs_total = summary_obj$summary.n_obs.total[i],
          bin_direction = "closer_to_expected_share_better",
          stringsAsFactors = FALSE
        ),
        base_row
      )
      row_id <- row_id + 1
    }
  }

  if (length(rows) == 0) {
    return(data.frame())
  }

  do.call(rbind, rows)
}
