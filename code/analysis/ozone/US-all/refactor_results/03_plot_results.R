# Plot builders for MRTS-focused US-all result diagnostics.

ensure_plot_packages <- function() {
  pkgs <- c("ggplot2", "dplyr", "tidyr", "scales")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing plotting packages: ",
      paste(missing, collapse = ", "),
      ". Please install them before running the plot runner."
    )
  }

  invisible(pkgs)
}

default_extreme_plot_probs <- function() {
  c(0.9, 0.95, 0.98, 0.99, 0.995)
}

parse_plot_summary_draws <- function(
    primary_env = "US_ALL_PLOT_SUMMARY_DRAWS",
    fallback_env = "US_ALL_SUMMARY_DRAWS",
    default = 400L
) {
  primary <- suppressWarnings(as.integer(Sys.getenv(primary_env, unset = "")))
  if (is.finite(primary) && primary > 0) {
    return(primary)
  }

  fallback <- suppressWarnings(as.integer(Sys.getenv(fallback_env, unset = "")))
  if (is.finite(fallback) && fallback > 0) {
    return(fallback)
  }

  as.integer(default)
}

resolve_plot_result_dirs <- function() {
  dirs_env <- trimws(Sys.getenv("US_ALL_RESULTS_DIRS", unset = ""))
  dir_primary <- trimws(Sys.getenv("US_ALL_RESULTS_DIR", unset = ""))

  dirs <- character(0)
  if (nzchar(dirs_env)) {
    dirs <- c(dirs, unlist(strsplit(dirs_env, ",", fixed = TRUE)))
  }
  if (nzchar(dir_primary)) {
    dirs <- c(dirs, dir_primary)
  }

  dirs <- c(dirs, "results", "results_new")
  dirs <- trimws(dirs)
  dirs <- unique(dirs[nzchar(dirs)])

  existing <- dirs[dir.exists(dirs)]
  if (length(existing) == 0) {
    stop("No result directory found for plot diagnostics. Checked: ", paste(dirs, collapse = ", "))
  }

  existing
}

resolve_plot_result_file <- function(setting_id, result_dirs) {
  candidates <- file.path(result_dirs, sprintf("us-all-%d.RData", setting_id))
  hit <- candidates[file.exists(candidates)]
  if (length(hit) > 0) {
    return(hit[1])
  }

  candidates[1]
}

format_basis_label <- function(method, knots, thresh) {
  sprintf("%s | knots=%s | thresh=%s", method, as.character(knots), as.character(thresh))
}

safe_rel_gap <- function(value, baseline_value) {
  if (!is.finite(value) || !is.finite(baseline_value) || baseline_value == 0) {
    return(NA_real_)
  }
  value / baseline_value - 1
}

build_mrts_focus_settings <- function(settings, comparison_paired_same_basis, available_settings = NULL) {
  settings_local <- settings
  if (!("setting_num" %in% names(settings_local)) && "setting" %in% names(settings_local)) {
    settings_local$setting_num <- safe_as_integer(settings_local$setting)
  }
  if (!("setting_num" %in% names(settings_local))) {
    stop("settings.csv must provide setting or setting_num for plot labeling.")
  }

  paired_cols <- c(
    "baseline_setting", "proposed_setting", "method", "knots",
    "thresh", "CMAQ", "TS", "proposed_lane"
  )
  paired_local <- unique(comparison_paired_same_basis[, paired_cols, drop = FALSE])

  if (!is.null(available_settings)) {
    paired_local <- paired_local[
      paired_local$baseline_setting %in% available_settings &
        paired_local$proposed_setting %in% available_settings,
      ,
      drop = FALSE
    ]
  }

  mrts_pairs <- paired_local[paired_local$proposed_lane == "mrts", , drop = FALSE]
  if (nrow(mrts_pairs) == 0) {
    return(data.frame())
  }

  baseline_ids <- sort(unique(mrts_pairs$baseline_setting))
  ar2_pairs <- paired_local[
    paired_local$proposed_lane == "ar2" &
      paired_local$baseline_setting %in% baseline_ids,
    ,
    drop = FALSE
  ]

  rows <- list()
  row_id <- 1L

  for (baseline_id in baseline_ids) {
    mrts_basis_rows <- mrts_pairs[mrts_pairs$baseline_setting == baseline_id, , drop = FALSE]
    basis_row <- mrts_basis_rows[1, , drop = FALSE]
    basis_label <- format_basis_label(basis_row$method[1], basis_row$knots[1], basis_row$thresh[1])

    rows[[row_id]] <- data.frame(
      setting = baseline_id,
      baseline_setting = baseline_id,
      basis_label = basis_label,
      method = basis_row$method[1],
      knots = basis_row$knots[1],
      thresh = basis_row$thresh[1],
      CMAQ = basis_row$CMAQ[1],
      TS = basis_row$TS[1],
      model_lane_focus = "baseline",
      short_label = "baseline",
      model_label = sprintf("Baseline (%d)", baseline_id),
      model_order = 0L,
      mrts_k = NA_integer_,
      stringsAsFactors = FALSE
    )
    row_id <- row_id + 1L

    ar2_hit <- ar2_pairs[ar2_pairs$baseline_setting == baseline_id, , drop = FALSE]
    if (nrow(ar2_hit) > 0) {
      ar2_setting <- sort(unique(ar2_hit$proposed_setting))[1]
      rows[[row_id]] <- data.frame(
        setting = ar2_setting,
        baseline_setting = baseline_id,
        basis_label = basis_label,
        method = basis_row$method[1],
        knots = basis_row$knots[1],
        thresh = basis_row$thresh[1],
        CMAQ = basis_row$CMAQ[1],
        TS = basis_row$TS[1],
        model_lane_focus = "ar2",
        short_label = "AR2",
        model_label = sprintf("AR2 (%d)", ar2_setting),
        model_order = 1L,
        mrts_k = NA_integer_,
        stringsAsFactors = FALSE
      )
      row_id <- row_id + 1L
    }

    mrts_setting_ids <- sort(unique(mrts_basis_rows$proposed_setting))
    for (mrts_setting in mrts_setting_ids) {
      meta <- find_setting_meta(settings_local, mrts_setting)
      mrts_k <- safe_as_integer(meta_value(meta, "mrts"))
      short_label <- if (is.na(mrts_k)) {
        sprintf("MRTS-%d", mrts_setting)
      } else {
        sprintf("k=%d", mrts_k)
      }
      order_value <- if (is.na(mrts_k)) 10L + mrts_setting else 10L + mrts_k

      rows[[row_id]] <- data.frame(
        setting = mrts_setting,
        baseline_setting = baseline_id,
        basis_label = basis_label,
        method = basis_row$method[1],
        knots = basis_row$knots[1],
        thresh = basis_row$thresh[1],
        CMAQ = basis_row$CMAQ[1],
        TS = basis_row$TS[1],
        model_lane_focus = "mrts",
        short_label = short_label,
        model_label = sprintf("MRTS k=%s (%d)", if (is.na(mrts_k)) "NA" else mrts_k, mrts_setting),
        model_order = order_value,
        mrts_k = mrts_k,
        stringsAsFactors = FALSE
      )
      row_id <- row_id + 1L
    }
  }

  focus <- do.call(rbind, rows)
  focus <- focus[!duplicated(focus$setting), , drop = FALSE]
  focus <- focus[order(focus$baseline_setting, focus$model_order, focus$setting), , drop = FALSE]
  rownames(focus) <- NULL
  focus
}

compute_setting_mean_prediction_metrics <- function(
    setting_id,
    result_path,
    Y,
    cv_lst,
    trans_setting_ids = 2L,
    summary_draws = 400L,
    enforce_contract = TRUE
) {
  if (!file.exists(result_path)) {
    return(NULL)
  }

  load_env <- new.env(parent = emptyenv())
  load(result_path, envir = load_env)

  if (!exists("fit", envir = load_env, inherits = FALSE)) {
    return(NULL)
  }

  fit <- get("fit", envir = load_env, inherits = FALSE)
  nsets <- length(cv_lst)
  if (!is.list(fit) || length(fit) < nsets) {
    return(NULL)
  }

  bias_fold <- rep(NA_real_, nsets)
  mae_fold <- rep(NA_real_, nsets)
  rmse_fold <- rep(NA_real_, nsets)
  pred_mean_fold <- rep(NA_real_, nsets)
  obs_mean_fold <- rep(NA_real_, nsets)
  n_obs_fold <- rep(NA_real_, nsets)
  draws_used_fold <- rep(NA_real_, nsets)

  trans <- setting_id %in% trans_setting_ids

  for (d in seq_len(nsets)) {
    fit_d <- fit[[d]]
    if (is.null(fit_d) || is.null(fit_d$yp)) {
      return(NULL)
    }

    validate <- Y[cv_lst[[d]], , drop = FALSE]
    pred_d <- fit_d$yp[, , ]

    needs_contract <- isTRUE(enforce_contract) && !trans
    if (needs_contract && !validate_pred_contract(pred_d, validate)) {
      return(NULL)
    }

    draw_idx <- select_summary_draw_indices(dim(pred_d)[1], summary_draws)
    pred_subset <- pred_d[draw_idx, , , drop = FALSE]
    pred_mean <- apply(pred_subset, c(2, 3), mean, na.rm = TRUE)

    if (is.null(dim(pred_mean))) {
      pred_mean <- matrix(pred_mean, nrow = dim(pred_subset)[2], ncol = dim(pred_subset)[3])
    }
    if (trans) {
      pred_mean <- t(pred_mean)
    }

    pred_vec <- as.vector(pred_mean)
    y_vec <- as.vector(validate)
    keep <- is.finite(y_vec) & is.finite(pred_vec)
    if (!any(keep)) {
      next
    }

    diff_vec <- pred_vec[keep] - y_vec[keep]
    bias_fold[d] <- mean(diff_vec, na.rm = TRUE)
    mae_fold[d] <- mean(abs(diff_vec), na.rm = TRUE)
    rmse_fold[d] <- sqrt(mean(diff_vec^2, na.rm = TRUE))
    pred_mean_fold[d] <- mean(pred_vec[keep], na.rm = TRUE)
    obs_mean_fold[d] <- mean(y_vec[keep], na.rm = TRUE)
    n_obs_fold[d] <- length(diff_vec)
    draws_used_fold[d] <- length(draw_idx)
  }

  data.frame(
    setting = setting_id,
    mean_bias = mean_if_any(bias_fold),
    mean_bias_se = se_if_any(bias_fold),
    mean_abs_bias = mean_if_any(abs(bias_fold)),
    mean_mae = mean_if_any(mae_fold),
    mean_mae_se = se_if_any(mae_fold),
    mean_rmse = mean_if_any(rmse_fold),
    mean_rmse_se = se_if_any(rmse_fold),
    predicted_mean = mean_if_any(pred_mean_fold),
    observed_mean = mean_if_any(obs_mean_fold),
    n_obs_total = sum_if_any(n_obs_fold),
    n_obs_mean_per_fold = mean_if_any(n_obs_fold),
    draws_used_mean = mean_if_any(draws_used_fold),
    stringsAsFactors = FALSE
  )
}

compute_focus_mean_prediction_metrics <- function(
    focus_settings,
    result_dirs,
    Y,
    cv_lst,
    trans_setting_ids = 2L,
    summary_draws = 400L,
    enforce_contract = TRUE
) {
  rows <- list()
  row_id <- 1L

  for (setting_id in focus_settings$setting) {
    result_path <- resolve_plot_result_file(setting_id, result_dirs = result_dirs)
    metrics <- compute_setting_mean_prediction_metrics(
      setting_id = setting_id,
      result_path = result_path,
      Y = Y,
      cv_lst = cv_lst,
      trans_setting_ids = trans_setting_ids,
      summary_draws = summary_draws,
      enforce_contract = enforce_contract
    )
    if (!is.null(metrics)) {
      rows[[row_id]] <- metrics
      row_id <- row_id + 1L
    }
  }

  if (length(rows) == 0) {
    return(data.frame())
  }

  metrics_df <- do.call(rbind, rows)
  out <- merge(focus_settings, metrics_df, by = "setting", all.x = FALSE, sort = FALSE)

  baseline_cols <- c(
    "baseline_setting", "mean_abs_bias", "mean_mae", "mean_rmse",
    "mean_bias", "predicted_mean", "observed_mean"
  )
  baseline_df <- out[out$model_lane_focus == "baseline", baseline_cols, drop = FALSE]
  names(baseline_df)[names(baseline_df) == "mean_abs_bias"] <- "baseline_abs_bias"
  names(baseline_df)[names(baseline_df) == "mean_mae"] <- "baseline_mae"
  names(baseline_df)[names(baseline_df) == "mean_rmse"] <- "baseline_rmse"
  names(baseline_df)[names(baseline_df) == "mean_bias"] <- "baseline_bias"
  names(baseline_df)[names(baseline_df) == "predicted_mean"] <- "baseline_predicted_mean"
  names(baseline_df)[names(baseline_df) == "observed_mean"] <- "baseline_observed_mean"

  out <- merge(out, baseline_df, by = "baseline_setting", all.x = TRUE, sort = FALSE)
  out$mean_abs_bias_rel_gap <- mapply(safe_rel_gap, out$mean_abs_bias, out$baseline_abs_bias)
  out$mean_mae_rel_gap <- mapply(safe_rel_gap, out$mean_mae, out$baseline_mae)
  out$mean_rmse_rel_gap <- mapply(safe_rel_gap, out$mean_rmse, out$baseline_rmse)
  out$mean_bias_delta <- out$mean_bias - out$baseline_bias
  out$predicted_mean_delta <- out$predicted_mean - out$baseline_predicted_mean
  out$observed_mean_delta <- out$observed_mean - out$baseline_observed_mean

  out <- out[order(out$baseline_setting, out$model_order, out$setting), , drop = FALSE]
  rownames(out) <- NULL
  out
}

plot_theme_us_all <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 8))
    )
}

build_focus_palette <- function(labels) {
  if (is.data.frame(labels)) {
    if ("short_label" %in% names(labels)) {
      labels <- labels$short_label
    } else if ("model_label" %in% names(labels)) {
      labels <- labels$model_label
    }
  }

  labels <- unique(as.character(labels))
  colors <- setNames(rep("#666666", length(labels)), labels)

  for (lbl in labels) {
    if (tolower(lbl) == "baseline" || grepl("^Baseline", lbl)) {
      colors[lbl] <- "#333333"
    } else if (lbl == "AR2" || grepl("^AR2", lbl)) {
      colors[lbl] <- "#1F78B4"
    } else if (grepl("k=5", lbl, fixed = TRUE)) {
      colors[lbl] <- "#D95F02"
    } else if (grepl("k=10", lbl, fixed = TRUE)) {
      colors[lbl] <- "#E6AB02"
    } else if (grepl("k=15", lbl, fixed = TRUE)) {
      colors[lbl] <- "#7570B3"
    } else if (grepl("k=800", lbl, fixed = TRUE)) {
      colors[lbl] <- "#1B9E77"
    }
  }

  colors
}

write_plot_file <- function(plot_obj, output_path, width, height) {
  ggplot2::ggsave(
    filename = output_path,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

plot_mrts_extreme_delta_profiles <- function(
    comparison_paired_same_basis,
    focus_settings,
    output_path,
    target_probs = default_extreme_plot_probs()
) {
  proposed_lookup <- focus_settings[
    focus_settings$model_lane_focus != "baseline",
    c("setting", "baseline_setting", "basis_label", "short_label", "model_order"),
    drop = FALSE
  ]
  if (nrow(proposed_lookup) == 0) {
    return(invisible(FALSE))
  }

  plot_df <- merge(
    comparison_paired_same_basis,
    proposed_lookup[, c("setting", "baseline_setting", "basis_label", "short_label", "model_order"), drop = FALSE],
    by.x = c("proposed_setting", "baseline_setting"),
    by.y = c("setting", "baseline_setting"),
    all = FALSE,
    sort = FALSE
  )
  plot_df <- plot_df[
    plot_df$metric %in% c("brier", "quantile") &
      plot_df$target_level %in% target_probs,
    ,
    drop = FALSE
  ]
  if (nrow(plot_df) == 0) {
    return(invisible(FALSE))
  }

  plot_df$metric_label <- ifelse(plot_df$metric == "brier", "Threshold Brier", "Quantile score")
  plot_df$short_label <- factor(
    plot_df$short_label,
    levels = unique(focus_settings$short_label[order(focus_settings$model_order, focus_settings$setting)])
  )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = target_level,
      y = delta_rel_score,
      color = short_label,
      group = proposed_setting
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "#555555", linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_grid(metric_label ~ basis_label, scales = "free_y") +
    ggplot2::scale_x_continuous(breaks = target_probs, labels = scales::label_number(accuracy = 0.001)) +
    ggplot2::scale_y_continuous(labels = scales::label_number(accuracy = 0.001)) +
    ggplot2::scale_color_manual(values = build_focus_palette(focus_settings$short_label)) +
    ggplot2::labs(
      title = "Extreme-score change versus paired baseline",
      subtitle = "Negative values indicate improvement over the matched TS baseline.",
      x = "Target quantile",
      y = "Delta relative score"
    ) +
    plot_theme_us_all()

  write_plot_file(p, output_path = output_path, width = 13, height = 7)
  invisible(TRUE)
}

plot_mrts_extreme_split_brier <- function(
    comparison_brier_split,
    focus_settings,
    output_path,
    target_probs = c(0.95, 0.98, 0.99, 0.995)
) {
  plot_df <- merge(
    comparison_brier_split,
    focus_settings[, c("setting", "baseline_setting", "basis_label", "short_label", "model_order"), drop = FALSE],
    by = "setting",
    all = FALSE,
    sort = FALSE
  )
  plot_df <- plot_df[
    plot_df$event_quantile %in% target_probs &
      plot_df$band_type %in% c("below_threshold", "above_threshold"),
    ,
    drop = FALSE
  ]
  if (nrow(plot_df) == 0) {
    return(invisible(FALSE))
  }

  baseline_df <- plot_df[
    plot_df$short_label %in% focus_settings$short_label[focus_settings$model_lane_focus == "baseline"],
    c("baseline_setting", "band_type", "event_quantile", "score_mean"),
    drop = FALSE
  ]
  names(baseline_df)[names(baseline_df) == "score_mean"] <- "baseline_score_mean"

  plot_df <- merge(
    plot_df,
    baseline_df,
    by = c("baseline_setting", "band_type", "event_quantile"),
    all.x = TRUE,
    sort = FALSE
  )
  plot_df$score_rel_gap_to_baseline <- mapply(
    safe_rel_gap,
    plot_df$score_mean,
    plot_df$baseline_score_mean
  )
  plot_df$band_label <- ifelse(
    plot_df$band_type == "above_threshold",
    "Above threshold (miss-focused)",
    "Below threshold (false-alarm-focused)"
  )
  plot_df$short_label <- factor(
    plot_df$short_label,
    levels = unique(focus_settings$short_label[order(focus_settings$model_order, focus_settings$setting)])
  )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = event_quantile,
      y = score_rel_gap_to_baseline,
      color = short_label,
      group = setting
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "#555555", linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_grid(band_label ~ basis_label, scales = "free_y") +
    ggplot2::scale_x_continuous(breaks = target_probs, labels = scales::label_number(accuracy = 0.001)) +
    ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
    ggplot2::scale_color_manual(values = build_focus_palette(focus_settings$short_label)) +
    ggplot2::labs(
      title = "Split-Brier change versus paired baseline",
      subtitle = "Negative percentages indicate lower split-Brier score than the matched TS baseline.",
      x = "Event quantile",
      y = "Relative gap to paired baseline"
    ) +
    plot_theme_us_all()

  write_plot_file(p, output_path = output_path, width = 13, height = 8)
  invisible(TRUE)
}

plot_mrts_mean_error_gap <- function(mean_metrics, output_path) {
  if (nrow(mean_metrics) == 0) {
    return(invisible(FALSE))
  }

  error_gap_df <- rbind(
    data.frame(
      basis_label = mean_metrics$basis_label,
      short_label = mean_metrics$short_label,
      model_order = mean_metrics$model_order,
      metric_label = "RMSE gap",
      value = mean_metrics$mean_rmse_rel_gap,
      stringsAsFactors = FALSE
    ),
    data.frame(
      basis_label = mean_metrics$basis_label,
      short_label = mean_metrics$short_label,
      model_order = mean_metrics$model_order,
      metric_label = "MAE gap",
      value = mean_metrics$mean_mae_rel_gap,
      stringsAsFactors = FALSE
    ),
    data.frame(
      basis_label = mean_metrics$basis_label,
      short_label = mean_metrics$short_label,
      model_order = mean_metrics$model_order,
      metric_label = "|Bias| gap",
      value = mean_metrics$mean_abs_bias_rel_gap,
      stringsAsFactors = FALSE
    )
  )

  error_gap_df$short_label <- factor(
    error_gap_df$short_label,
    levels = unique(mean_metrics$short_label[order(mean_metrics$model_order, mean_metrics$setting)])
  )

  p <- ggplot2::ggplot(
    error_gap_df,
    ggplot2::aes(x = short_label, y = value, color = short_label)
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "#555555", linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_segment(
      ggplot2::aes(xend = short_label, y = 0, yend = value),
      linewidth = 0.7,
      alpha = 0.8
    ) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::facet_grid(metric_label ~ basis_label, scales = "free_y") +
    ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
    ggplot2::scale_color_manual(values = build_focus_palette(mean_metrics$short_label)) +
    ggplot2::labs(
      title = "Mean-level error change versus paired baseline",
      subtitle = "Negative percentages indicate lower average error than the matched TS baseline.",
      x = NULL,
      y = "Relative gap to paired baseline"
    ) +
    plot_theme_us_all() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

  write_plot_file(p, output_path = output_path, width = 12, height = 8)
  invisible(TRUE)
}

plot_mrts_mean_bias <- function(mean_metrics, output_path) {
  if (nrow(mean_metrics) == 0) {
    return(invisible(FALSE))
  }

  mean_metrics$short_label <- factor(
    mean_metrics$short_label,
    levels = unique(mean_metrics$short_label[order(mean_metrics$model_order, mean_metrics$setting)])
  )

  p <- ggplot2::ggplot(
    mean_metrics,
    ggplot2::aes(x = short_label, y = mean_bias, color = short_label)
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "#555555", linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = mean_bias - mean_bias_se, ymax = mean_bias + mean_bias_se),
      linewidth = 0.7,
      alpha = 0.9
    ) +
    ggplot2::facet_wrap(~ basis_label, scales = "free_y") +
    ggplot2::scale_color_manual(values = build_focus_palette(mean_metrics$short_label)) +
    ggplot2::labs(
      title = "Mean-prediction bias by paired basis",
      subtitle = "Bias is computed as predictive mean minus observed value.",
      x = NULL,
      y = "Mean bias"
    ) +
    plot_theme_us_all() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

  write_plot_file(p, output_path = output_path, width = 12, height = 4.6)
  invisible(TRUE)
}

plot_mrts_tail_mean_tradeoff <- function(
    comparison_paired_same_basis,
    mean_metrics,
    output_path,
    target_prob = 0.99,
    metric = "brier"
) {
  tail_df <- comparison_paired_same_basis[
    comparison_paired_same_basis$metric == metric &
      abs(comparison_paired_same_basis$target_level - target_prob) < 1e-12,
    ,
    drop = FALSE
  ]
  if (nrow(tail_df) == 0 || nrow(mean_metrics) == 0) {
    return(invisible(FALSE))
  }

  tail_df <- merge(
    tail_df,
    mean_metrics[
      mean_metrics$model_lane_focus != "baseline",
      c("setting", "basis_label", "short_label", "mean_rmse_rel_gap"),
      drop = FALSE
    ],
    by.x = "proposed_setting",
    by.y = "setting",
    all = FALSE,
    sort = FALSE
  )
  if (nrow(tail_df) == 0) {
    return(invisible(FALSE))
  }

  tail_df$short_label <- factor(
    tail_df$short_label,
    levels = unique(mean_metrics$short_label[order(mean_metrics$model_order, mean_metrics$setting)])
  )

  p <- ggplot2::ggplot(
    tail_df,
    ggplot2::aes(
      x = mean_rmse_rel_gap,
      y = delta_rel_score,
      color = short_label
    )
  ) +
    ggplot2::geom_vline(xintercept = 0, color = "#555555", linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = 0, color = "#555555", linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::geom_text(
      ggplot2::aes(label = short_label),
      nudge_y = 0.002,
      size = 3,
      show.legend = FALSE,
      check_overlap = TRUE
    ) +
    ggplot2::facet_wrap(~ basis_label, scales = "free") +
    ggplot2::scale_x_continuous(labels = scales::label_percent(accuracy = 0.01)) +
    ggplot2::scale_y_continuous(labels = scales::label_number(accuracy = 0.001)) +
    ggplot2::scale_color_manual(values = build_focus_palette(mean_metrics$short_label)) +
    ggplot2::labs(
      title = sprintf("Tail-versus-mean trade-off at q = %.3f", target_prob),
      subtitle = "Left/down is better: lower mean RMSE and lower tail relative score than the paired TS baseline.",
      x = "Mean RMSE relative gap to paired baseline",
      y = "Tail relative-score delta"
    ) +
    plot_theme_us_all() +
    ggplot2::theme(legend.position = "none")

  write_plot_file(p, output_path = output_path, width = 12, height = 4.6)
  invisible(TRUE)
}
