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

build_comparison_full_table <- function(summary_obj, settings, baseline_ids, proposed_ids) {
  rows <- list()
  row_id <- 1

  for (i in summary_obj$available_settings) {
    meta <- find_setting_meta(settings, i)

    for (q in seq_along(summary_obj$probs)) {
      rows[[row_id]] <- data.frame(
        setting = i,
        metric = "brier",
        quantile = summary_obj$probs[q],
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

      rows[[row_id]] <- data.frame(
        setting = i,
        metric = "quantile",
        quantile = summary_obj$probs[q],
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
    target_quantiles = c(0.95, 0.98, 0.99, 0.995),
    candidate_settings = NULL,
    metric = c("brier", "quantile")
) {
  metric <- match.arg(metric)
  score_mat <- if (metric == "brier") summary_obj$bs.mean.ref.gau else summary_obj$qs.mean.ref.gau

  if (is.null(candidate_settings)) {
    candidate_settings <- summary_obj$available_settings
  }

  rows <- list()
  row_id <- 1

  for (q in target_quantiles) {
    q_idx <- which(abs(summary_obj$probs - q) < 1e-12)
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

    for (q in seq_along(summary_obj$probs)) {
      rows[[row_id]] <- data.frame(
        quantile = summary_obj$probs[q],
        baseline_setting = b_setting,
        proposed_setting = p_setting,
        baseline_rel_brier = summary_obj$bs.mean.ref.gau[b_setting, q],
        proposed_rel_brier = summary_obj$bs.mean.ref.gau[p_setting, q],
        delta_rel_brier = summary_obj$bs.mean.ref.gau[p_setting, q] - summary_obj$bs.mean.ref.gau[b_setting, q],
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
