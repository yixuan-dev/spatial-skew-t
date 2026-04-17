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
            metric = "brier_split",
            split_scheme = "same_threshold",
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
