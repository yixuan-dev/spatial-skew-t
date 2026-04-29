# Shared scoring engine for US-all result files.

validate_pred_contract <- function(pred, validate) {
  if (is.null(pred)) {
    return(FALSE)
  }
  if (length(dim(pred)) != 3) {
    return(FALSE)
  }

  pred_dims <- dim(pred)
  val_dims <- dim(validate)
  if (length(val_dims) != 2) {
    return(FALSE)
  }

  pred_dims[2] == val_dims[1] && pred_dims[3] == val_dims[2]
}

default_uncertainty_levels <- function() {
  c(0.8, 0.9, 0.95)
}

default_pit_breaks <- function() {
  seq(0, 1, by = 0.1)
}

mean_if_any <- function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

sum_if_any <- function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }
  sum(x, na.rm = TRUE)
}

se_if_any <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  if (length(x) == 1) {
    return(0)
  }
  stats::sd(x) / sqrt(length(x))
}

select_summary_draw_indices <- function(n_draws, max_draws = NA_integer_) {
  if (!is.finite(max_draws) || max_draws <= 0 || max_draws >= n_draws) {
    return(seq_len(n_draws))
  }

  idx <- unique(as.integer(round(seq(1, n_draws, length.out = max_draws))))
  idx[idx >= 1 & idx <= n_draws]
}

sample_quantile_from_sorted <- function(sorted_draws, prob) {
  if (is.null(dim(sorted_draws))) {
    sorted_draws <- matrix(sorted_draws, ncol = 1)
  }

  n_draws <- nrow(sorted_draws)
  if (n_draws == 0) {
    return(rep(NA_real_, ncol(sorted_draws)))
  }
  if (n_draws == 1) {
    return(sorted_draws[1, ])
  }

  h <- (n_draws - 1) * prob + 1
  lo <- floor(h)
  hi <- ceiling(h)
  w_hi <- h - lo

  (1 - w_hi) * sorted_draws[lo, ] + w_hi * sorted_draws[hi, ]
}

pit_uniform_ks_stat <- function(pit) {
  pit <- pit[is.finite(pit)]
  if (length(pit) == 0) {
    return(NA_real_)
  }

  pit <- sort(pit)
  n <- length(pit)
  upper_gap <- abs(seq_len(n) / n - pit)
  lower_gap <- abs((seq_len(n) - 1) / n - pit)
  max(c(upper_gap, lower_gap), na.rm = TRUE)
}

default_brier_split_target_probs <- function() {
  c(0.9, 0.95, 0.98, 0.99, 0.995)
}

default_classification_target_probs <- function() {
  default_brier_split_target_probs()
}

default_classification_probability_cutoff <- function() {
  0.5
}

match_prob_levels <- function(level_grid, target_levels, tol = 1e-12) {
  idx <- rep(NA_integer_, length(target_levels))
  for (k in seq_along(target_levels)) {
    hit <- which(abs(level_grid - target_levels[k]) < tol)
    if (length(hit) > 0) {
      idx[k] <- hit[1]
    }
  }
  idx
}

compute_split_brier_scores <- function(
    preds,
    validate,
    thresholds,
    target_idx,
    trans = FALSE
) {
  band_names <- c("all", "below_threshold", "above_threshold")
  scores <- matrix(
    NA_real_,
    nrow = length(band_names),
    ncol = length(target_idx),
    dimnames = list(band_names, NULL)
  )
  counts <- scores
  shares <- scores

  y_all <- as.vector(validate)

  for (k in seq_along(target_idx)) {
    threshold_id <- target_idx[k]
    threshold_value <- thresholds[threshold_id]
    pat <- apply((preds > threshold_value), c(2, 3), mean)
    if (trans) {
      pat <- t(pat)
    }

    p_all <- as.vector(pat)
    keep <- is.finite(y_all) & is.finite(p_all)
    if (!any(keep)) {
      next
    }

    y_vec <- y_all[keep]
    p_vec <- p_all[keep]
    is_event <- y_vec > threshold_value
    sq_err <- (as.numeric(is_event) - p_vec)^2
    total_n <- length(sq_err)
    if (total_n == 0) {
      next
    }

    scores["all", k] <- mean(sq_err, na.rm = TRUE)
    counts["all", k] <- total_n
    shares["all", k] <- 1

    below_idx <- !is_event
    below_n <- sum(below_idx)
    counts["below_threshold", k] <- below_n
    shares["below_threshold", k] <- below_n / total_n
    if (below_n > 0) {
      scores["below_threshold", k] <- mean(sq_err[below_idx], na.rm = TRUE)
    }

    above_idx <- is_event
    above_n <- sum(above_idx)
    counts["above_threshold", k] <- above_n
    shares["above_threshold", k] <- above_n / total_n
    if (above_n > 0) {
      scores["above_threshold", k] <- mean(sq_err[above_idx], na.rm = TRUE)
    }
  }

  list(
    band_names = band_names,
    scores = scores,
    counts = counts,
    shares = shares
  )
}

compute_classification_metrics <- function(
    preds,
    validate,
    thresholds,
    target_idx,
    trans = FALSE,
    probability_cutoff = default_classification_probability_cutoff(),
    max_draws = NA_integer_
) {
  count_names <- c("tp", "tn", "fp", "fn")
  metric_names <- c("accuracy", "precision", "recall", "specificity", "f1")
  counts <- matrix(
    NA_real_,
    nrow = length(count_names),
    ncol = length(target_idx),
    dimnames = list(count_names, NULL)
  )
  metrics <- matrix(
    NA_real_,
    nrow = length(metric_names),
    ncol = length(target_idx),
    dimnames = list(metric_names, NULL)
  )
  obs_total <- rep(NA_real_, length(target_idx))
  actual_positive_share <- rep(NA_real_, length(target_idx))
  predicted_positive_share <- rep(NA_real_, length(target_idx))
  draw_idx <- select_summary_draw_indices(dim(preds)[1], max_draws)
  preds_local <- preds[draw_idx, , , drop = FALSE]

  y_all <- as.vector(validate)

  safe_ratio <- function(num, denom) {
    if (!is.finite(denom) || denom <= 0) {
      return(NA_real_)
    }
    num / denom
  }

  for (k in seq_along(target_idx)) {
    threshold_id <- target_idx[k]
    threshold_value <- thresholds[threshold_id]
    pat <- apply((preds_local > threshold_value), c(2, 3), mean)
    if (trans) {
      pat <- t(pat)
    }

    p_all <- as.vector(pat)
    keep <- is.finite(y_all) & is.finite(p_all)
    if (!any(keep)) {
      next
    }

    y_vec <- y_all[keep]
    p_vec <- p_all[keep]
    actual_positive <- y_vec > threshold_value
    predicted_positive <- p_vec >= probability_cutoff

    tp <- sum(predicted_positive & actual_positive)
    tn <- sum(!predicted_positive & !actual_positive)
    fp <- sum(predicted_positive & !actual_positive)
    fn <- sum(!predicted_positive & actual_positive)
    total_n <- tp + tn + fp + fn

    counts[, k] <- c(tp, tn, fp, fn)
    obs_total[k] <- total_n
    actual_positive_share[k] <- safe_ratio(tp + fn, total_n)
    predicted_positive_share[k] <- safe_ratio(tp + fp, total_n)

    precision <- safe_ratio(tp, tp + fp)
    recall <- safe_ratio(tp, tp + fn)

    metrics["accuracy", k] <- safe_ratio(tp + tn, total_n)
    metrics["precision", k] <- precision
    metrics["recall", k] <- recall
    metrics["specificity", k] <- safe_ratio(tn, tn + fp)
    metrics["f1", k] <- if (is.finite(precision) && is.finite(recall) && (precision + recall) > 0) {
      2 * precision * recall / (precision + recall)
    } else {
      NA_real_
    }
  }

  list(
    count_names = count_names,
    metric_names = metric_names,
    counts = counts,
    metrics = metrics,
    obs_total = obs_total,
    actual_positive_share = actual_positive_share,
    predicted_positive_share = predicted_positive_share,
    probability_cutoff = probability_cutoff,
    draws_used = length(draw_idx)
  )
}

compute_point_prediction_metrics <- function(
    preds,
    validate,
    trans = FALSE
) {
  pred_mean <- apply(preds, c(2, 3), mean, na.rm = TRUE)
  if (trans) {
    pred_mean <- t(pred_mean)
  }

  y_vec <- as.vector(validate)
  pred_vec <- as.vector(pred_mean)
  keep <- is.finite(y_vec) & is.finite(pred_vec)
  if (!any(keep)) {
    return(list(mspe = NA_real_, mape = NA_real_))
  }

  y_keep <- y_vec[keep]
  pred_keep <- pred_vec[keep]
  err <- pred_keep - y_keep
  denom_keep <- is.finite(y_keep) & y_keep != 0

  list(
    mspe = mean(err^2, na.rm = TRUE),
    mape = if (any(denom_keep)) {
      mean(abs(err[denom_keep] / y_keep[denom_keep]), na.rm = TRUE) * 100
    } else {
      NA_real_
    }
  )
}

compute_predictive_diagnostics <- function(
    preds,
    validate,
    trans = FALSE,
    summary_draws = 400L,
    uncertainty_levels = default_uncertainty_levels(),
    pit_breaks = default_pit_breaks()
) {
  draw_idx <- select_summary_draw_indices(dim(preds)[1], summary_draws)

  pred_subset <- preds[draw_idx, , , drop = FALSE]
  if (trans) {
    pred_subset <- aperm(pred_subset, c(1, 3, 2))
  }

  pred_mat <- matrix(pred_subset, nrow = length(draw_idx))
  y_vec <- as.vector(validate)

  keep_obs <- is.finite(y_vec)
  if (!all(keep_obs)) {
    pred_mat <- pred_mat[, keep_obs, drop = FALSE]
    y_vec <- y_vec[keep_obs]
  }

  n_obs <- length(y_vec)
  draws_used <- nrow(pred_mat)
  nbins <- length(pit_breaks) - 1L

  if (n_obs == 0 || draws_used == 0) {
    return(list(
      crps_mean = NA_real_,
      coverage = rep(NA_real_, length(uncertainty_levels)),
      pit_mean = NA_real_,
      pit_variance = NA_real_,
      pit_ks = NA_real_,
      pit_bin_share = rep(NA_real_, nbins),
      pit_bin_count = rep(NA_real_, nbins),
      pit_mae = NA_real_,
      pit_rmse = NA_real_,
      draws_used = draws_used,
      n_obs = n_obs
    ))
  }

  pred_sorted <- apply(pred_mat, 2, sort, na.last = NA)
  if (n_obs == 1L) {
    pred_sorted <- matrix(pred_sorted, nrow = draws_used, ncol = 1L)
  }

  y_mat <- matrix(y_vec, nrow = draws_used, ncol = n_obs, byrow = TRUE)
  crps_first_term <- colMeans(abs(pred_mat - y_mat), na.rm = TRUE)
  crps_weights <- (2 * seq_len(draws_used) - draws_used - 1) / (draws_used^2)
  crps_second_term <- colSums(pred_sorted * matrix(crps_weights, nrow = draws_used, ncol = n_obs), na.rm = TRUE)
  crps_obs <- crps_first_term - crps_second_term

  coverage <- rep(NA_real_, length(uncertainty_levels))
  for (k in seq_along(uncertainty_levels)) {
    level <- uncertainty_levels[k]
    lower <- sample_quantile_from_sorted(pred_sorted, (1 - level) / 2)
    upper <- sample_quantile_from_sorted(pred_sorted, (1 + level) / 2)
    coverage[k] <- mean(y_vec >= lower & y_vec <= upper, na.rm = TRUE)
  }

  pit <- numeric(n_obs)
  for (j in seq_len(n_obs)) {
    pit[j] <- findInterval(y_vec[j], pred_sorted[, j], rightmost.closed = TRUE) / draws_used
  }

  pit_bin_id <- cut(
    pit,
    breaks = pit_breaks,
    include.lowest = TRUE,
    right = TRUE,
    labels = FALSE
  )
  pit_bin_count <- tabulate(pit_bin_id, nbins = nbins)
  pit_bin_share <- pit_bin_count / sum(pit_bin_count)
  expected_share <- diff(pit_breaks)

  list(
    crps_mean = mean(crps_obs, na.rm = TRUE),
    coverage = coverage,
    pit_mean = mean(pit, na.rm = TRUE),
    pit_variance = if (length(pit) > 1L) stats::var(pit) else 0,
    pit_ks = pit_uniform_ks_stat(pit),
    pit_bin_share = pit_bin_share,
    pit_bin_count = pit_bin_count,
    pit_mae = mean(abs(pit_bin_share - expected_share), na.rm = TRUE),
    pit_rmse = sqrt(mean((pit_bin_share - expected_share)^2, na.rm = TRUE)),
    draws_used = draws_used,
    n_obs = n_obs
  )
}

compute_us_all_scores <- function(
    setting_ids,
    result_path_fn,
    Y,
    cv_lst,
    probs,
    thresholds,
    threshold_probs = NULL,
    compute_brier_split_diagnostics = FALSE,
    brier_split_target_probs = default_brier_split_target_probs(),
    compute_classification_diagnostics = FALSE,
    classification_target_probs = default_classification_target_probs(),
    classification_probability_cutoff = default_classification_probability_cutoff(),
    trans_setting_ids = 2L,
    enforce_contract = TRUE,
    compute_uncertainty_diagnostics = FALSE,
    summary_draws = 400L,
    uncertainty_levels = default_uncertainty_levels(),
    pit_breaks = default_pit_breaks(),
    partial_ok_setting_ids = integer(0)
) {
  setting_ids <- sort(unique(setting_ids[!is.na(setting_ids)]))
  if (length(setting_ids) == 0) {
    stop("No valid setting IDs were provided.")
  }

  nsets <- length(cv_lst)
  max_setting <- max(setting_ids)

  if (is.null(threshold_probs)) {
    threshold_probs <- probs
  }
  if (length(threshold_probs) != length(thresholds)) {
    stop("threshold_probs and thresholds must have the same length.")
  }

  brier_split_band_names <- c("all", "below_threshold", "above_threshold")
  classification_count_names <- c("tp", "tn", "fp", "fn")
  classification_metric_names <- c("accuracy", "precision", "recall", "specificity", "f1")
  brier_split_target_idx <- integer(0)
  if (isTRUE(compute_brier_split_diagnostics)) {
    brier_split_target_probs <- as.numeric(brier_split_target_probs)
    brier_split_target_idx <- match_prob_levels(
      level_grid = threshold_probs,
      target_levels = brier_split_target_probs
    )
    if (anyNA(brier_split_target_idx)) {
      missing_probs <- brier_split_target_probs[is.na(brier_split_target_idx)]
      stop(
        "Missing split-Brier threshold_probs: ",
        paste(missing_probs, collapse = ", ")
      )
    }
  }
  classification_target_idx <- integer(0)
  if (isTRUE(compute_classification_diagnostics)) {
    classification_target_probs <- as.numeric(classification_target_probs)
    classification_target_idx <- match_prob_levels(
      level_grid = threshold_probs,
      target_levels = classification_target_probs
    )
    if (anyNA(classification_target_idx)) {
      missing_probs <- classification_target_probs[is.na(classification_target_idx)]
      stop(
        "Missing classification threshold_probs: ",
        paste(missing_probs, collapse = ", ")
      )
    }
  }

  quant.score <- array(NA_real_, dim = c(length(probs), nsets, max_setting))
  brier.score <- array(NA_real_, dim = c(length(thresholds), nsets, max_setting))
  brier.split.score <- if (isTRUE(compute_brier_split_diagnostics)) {
    array(
      NA_real_,
      dim = c(length(brier_split_band_names), length(brier_split_target_probs), nsets, max_setting)
    )
  } else {
    NULL
  }
  brier.split.count <- if (isTRUE(compute_brier_split_diagnostics)) {
    array(
      NA_real_,
      dim = c(length(brier_split_band_names), length(brier_split_target_probs), nsets, max_setting)
    )
  } else {
    NULL
  }
  brier.split.share <- if (isTRUE(compute_brier_split_diagnostics)) {
    array(
      NA_real_,
      dim = c(length(brier_split_band_names), length(brier_split_target_probs), nsets, max_setting)
    )
  } else {
    NULL
  }
  classification.count <- if (isTRUE(compute_classification_diagnostics)) {
    array(
      NA_real_,
      dim = c(length(classification_count_names), length(classification_target_probs), nsets, max_setting),
      dimnames = list(classification_count_names, as.character(classification_target_probs), NULL, NULL)
    )
  } else {
    NULL
  }
  classification.metric <- if (isTRUE(compute_classification_diagnostics)) {
    array(
      NA_real_,
      dim = c(length(classification_metric_names), length(classification_target_probs), nsets, max_setting),
      dimnames = list(classification_metric_names, as.character(classification_target_probs), NULL, NULL)
    )
  } else {
    NULL
  }
  classification.obs.total <- if (isTRUE(compute_classification_diagnostics)) {
    array(
      NA_real_,
      dim = c(length(classification_target_probs), nsets, max_setting),
      dimnames = list(as.character(classification_target_probs), NULL, NULL)
    )
  } else {
    NULL
  }
  classification.actual_positive_share <- if (isTRUE(compute_classification_diagnostics)) {
    array(
      NA_real_,
      dim = c(length(classification_target_probs), nsets, max_setting),
      dimnames = list(as.character(classification_target_probs), NULL, NULL)
    )
  } else {
    NULL
  }
  classification.predicted_positive_share <- if (isTRUE(compute_classification_diagnostics)) {
    array(
      NA_real_,
      dim = c(length(classification_target_probs), nsets, max_setting),
      dimnames = list(as.character(classification_target_probs), NULL, NULL)
    )
  } else {
    NULL
  }
  classification.draws.used <- if (isTRUE(compute_classification_diagnostics)) {
    matrix(NA_real_, nrow = nsets, ncol = max_setting)
  } else {
    NULL
  }
  crps.score <- if (isTRUE(compute_uncertainty_diagnostics)) matrix(NA_real_, nrow = nsets, ncol = max_setting) else NULL
  mspe.score <- matrix(NA_real_, nrow = nsets, ncol = max_setting)
  mape.score <- matrix(NA_real_, nrow = nsets, ncol = max_setting)
  coverage.score <- if (isTRUE(compute_uncertainty_diagnostics)) array(NA_real_, dim = c(length(uncertainty_levels), nsets, max_setting)) else NULL
  pit.mean.score <- if (isTRUE(compute_uncertainty_diagnostics)) matrix(NA_real_, nrow = nsets, ncol = max_setting) else NULL
  pit.variance.score <- if (isTRUE(compute_uncertainty_diagnostics)) matrix(NA_real_, nrow = nsets, ncol = max_setting) else NULL
  pit.ks.score <- if (isTRUE(compute_uncertainty_diagnostics)) matrix(NA_real_, nrow = nsets, ncol = max_setting) else NULL
  pit.mae.score <- if (isTRUE(compute_uncertainty_diagnostics)) matrix(NA_real_, nrow = nsets, ncol = max_setting) else NULL
  pit.rmse.score <- if (isTRUE(compute_uncertainty_diagnostics)) matrix(NA_real_, nrow = nsets, ncol = max_setting) else NULL
  pit.bin.share <- if (isTRUE(compute_uncertainty_diagnostics)) array(NA_real_, dim = c(length(pit_breaks) - 1L, nsets, max_setting)) else NULL
  pit.bin.count <- if (isTRUE(compute_uncertainty_diagnostics)) array(NA_real_, dim = c(length(pit_breaks) - 1L, nsets, max_setting)) else NULL
  summary.draws.used <- if (isTRUE(compute_uncertainty_diagnostics)) matrix(NA_real_, nrow = nsets, ncol = max_setting) else NULL
  summary.n_obs <- if (isTRUE(compute_uncertainty_diagnostics)) matrix(NA_real_, nrow = nsets, ncol = max_setting) else NULL

  available_settings <- integer(0)
  skipped_missing_file <- integer(0)
  skipped_bad_contract <- integer(0)
  skipped_scoring_error <- integer(0)

  for (i in setting_ids) {
    file <- result_path_fn(i)
    cat("start file", file, "\n")

    if (!file.exists(file)) {
      skipped_missing_file <- c(skipped_missing_file, i)
      cat("skip file (missing):", file, "\n")
      next
    }

    load_env <- new.env(parent = emptyenv())
    load(file, envir = load_env)

    if (!exists("fit", envir = load_env, inherits = FALSE)) {
      skipped_bad_contract <- c(skipped_bad_contract, i)
      cat("skip setting", i, ": object 'fit' not found\n")
      next
    }

    fit <- get("fit", envir = load_env, inherits = FALSE)
    if (!is.list(fit) || length(fit) < nsets) {
      skipped_bad_contract <- c(skipped_bad_contract, i)
      cat("skip setting", i, ": fit list length mismatch\n")
      next
    }

    has_bad_fold <- FALSE
    has_scoring_error <- FALSE
    score_error_message <- ""
    is_partial_ok <- i %in% partial_ok_setting_ids

    for (d in seq_len(nsets)) {
      fit.d <- fit[[d]]
      if (is.null(fit.d) || is.null(fit.d$yp)) {
        if (is_partial_ok) {
          cat("setting", i, "fold", d, ": NULL fold skipped (partial_ok)\n")
          next
        }
        has_bad_fold <- TRUE
        break
      }

      val.idx <- cv_lst[[d]]
      validate <- Y[val.idx, , drop = FALSE]
      pred.d <- fit.d$yp[, , ]

      needs_contract <- isTRUE(enforce_contract) && !(i %in% trans_setting_ids)
      if (needs_contract && !validate_pred_contract(pred.d, validate)) {
        has_bad_fold <- TRUE
        break
      }

      trans <- i %in% trans_setting_ids
      score_try <- tryCatch(
        {
          brier_split <- if (isTRUE(compute_brier_split_diagnostics)) {
            compute_split_brier_scores(
              preds = pred.d,
              validate = validate,
              thresholds = thresholds,
              target_idx = brier_split_target_idx,
              trans = trans
            )
          } else {
            NULL
          }

          classification <- if (isTRUE(compute_classification_diagnostics)) {
            compute_classification_metrics(
              preds = pred.d,
              validate = validate,
              thresholds = thresholds,
              target_idx = classification_target_idx,
              trans = trans,
              probability_cutoff = classification_probability_cutoff,
              max_draws = summary_draws
            )
          } else {
            NULL
          }

          diagnostics <- if (isTRUE(compute_uncertainty_diagnostics)) {
            compute_predictive_diagnostics(
              preds = pred.d,
              validate = validate,
              trans = trans,
              summary_draws = summary_draws,
              uncertainty_levels = uncertainty_levels,
              pit_breaks = pit_breaks
            )
          } else {
            NULL
          }

          list(
            quant = QuantScore(pred.d, probs, validate, trans = trans),
            brier = BrierScore(pred.d, thresholds, validate, trans = trans),
            point = compute_point_prediction_metrics(pred.d, validate, trans = trans),
            brier_split = brier_split,
            classification = classification,
            diagnostics = diagnostics
          )
        },
        error = function(e) {
          e
        }
      )

      if (inherits(score_try, "error")) {
        has_scoring_error <- TRUE
        score_error_message <- conditionMessage(score_try)
        break
      }

      quant.score[, d, i] <- score_try$quant
      brier.score[, d, i] <- score_try$brier
      mspe.score[d, i] <- score_try$point$mspe
      mape.score[d, i] <- score_try$point$mape
      if (isTRUE(compute_brier_split_diagnostics)) {
        brier.split.score[, , d, i] <- score_try$brier_split$scores
        brier.split.count[, , d, i] <- score_try$brier_split$counts
        brier.split.share[, , d, i] <- score_try$brier_split$shares
      }
      if (isTRUE(compute_classification_diagnostics)) {
        classification.count[, , d, i] <- score_try$classification$counts
        classification.metric[, , d, i] <- score_try$classification$metrics
        classification.obs.total[, d, i] <- score_try$classification$obs_total
        classification.actual_positive_share[, d, i] <- score_try$classification$actual_positive_share
        classification.predicted_positive_share[, d, i] <- score_try$classification$predicted_positive_share
        classification.draws.used[d, i] <- score_try$classification$draws_used
      }
      if (isTRUE(compute_uncertainty_diagnostics)) {
        crps.score[d, i] <- score_try$diagnostics$crps_mean
        coverage.score[, d, i] <- score_try$diagnostics$coverage
        pit.mean.score[d, i] <- score_try$diagnostics$pit_mean
        pit.variance.score[d, i] <- score_try$diagnostics$pit_variance
        pit.ks.score[d, i] <- score_try$diagnostics$pit_ks
        pit.mae.score[d, i] <- score_try$diagnostics$pit_mae
        pit.rmse.score[d, i] <- score_try$diagnostics$pit_rmse
        pit.bin.share[, d, i] <- score_try$diagnostics$pit_bin_share
        pit.bin.count[, d, i] <- score_try$diagnostics$pit_bin_count
        summary.draws.used[d, i] <- score_try$diagnostics$draws_used
        summary.n_obs[d, i] <- score_try$diagnostics$n_obs
      }
    }

    if (has_scoring_error) {
      skipped_scoring_error <- c(skipped_scoring_error, i)
      cat("skip setting", i, ": scoring error ->", score_error_message, "\n")
      next
    }

    if (has_bad_fold) {
      skipped_bad_contract <- c(skipped_bad_contract, i)
      cat("skip setting", i, ": fit[[d]]$yp contract mismatch\n")
      next
    }

    available_settings <- c(available_settings, i)
    cat("finish file", file, "\n")
  }

  list(
    setting_ids = setting_ids,
    nsets = nsets,
    max_setting = max_setting,
    probs = probs,
    thresholds = thresholds,
    threshold_probs = threshold_probs,
    quant.score = quant.score,
    brier.score = brier.score,
    brier.split.target_probs = brier_split_target_probs,
    brier.split.target_idx = brier_split_target_idx,
    brier.split.target_thresholds = if (length(brier_split_target_idx) > 0) thresholds[brier_split_target_idx] else numeric(0),
    brier.split.band_names = brier_split_band_names,
    brier.split.score = brier.split.score,
    brier.split.count = brier.split.count,
    brier.split.share = brier.split.share,
    classification.target_probs = classification_target_probs,
    classification.target_idx = classification_target_idx,
    classification.target_thresholds = if (length(classification_target_idx) > 0) thresholds[classification_target_idx] else numeric(0),
    classification.metric_names = classification_metric_names,
    classification.count_names = classification_count_names,
    classification.metric = classification.metric,
    classification.count = classification.count,
    classification.obs.total = classification.obs.total,
    classification.actual_positive_share = classification.actual_positive_share,
    classification.predicted_positive_share = classification.predicted_positive_share,
    classification.draws.used = classification.draws.used,
    classification.probability_cutoff = classification_probability_cutoff,
    uncertainty_levels = uncertainty_levels,
    pit_breaks = pit_breaks,
    crps.score = crps.score,
    mspe.score = mspe.score,
    mape.score = mape.score,
    coverage.score = coverage.score,
    pit.mean.score = pit.mean.score,
    pit.variance.score = pit.variance.score,
    pit.ks.score = pit.ks.score,
    pit.mae.score = pit.mae.score,
    pit.rmse.score = pit.rmse.score,
    pit.bin.share = pit.bin.share,
    pit.bin.count = pit.bin.count,
    summary.draws.used = summary.draws.used,
    summary.n_obs = summary.n_obs,
    available_settings = sort(unique(available_settings)),
    skipped_missing_file = sort(unique(skipped_missing_file)),
    skipped_bad_contract = sort(unique(skipped_bad_contract)),
    skipped_scoring_error = sort(unique(skipped_scoring_error))
  )
}

summarize_us_all_scores <- function(score_obj, baseline_setting = 1L) {
  max_setting <- score_obj$max_setting
  probs <- score_obj$probs
  thresholds <- score_obj$thresholds
  threshold_probs <- score_obj$threshold_probs
  nsets <- score_obj$nsets
  available_settings <- score_obj$available_settings

  if (!baseline_setting %in% available_settings) {
    stop(
      "Baseline setting ",
      baseline_setting,
      " is missing; cannot compute relative scores."
    )
  }

  quant.score.mean <- matrix(NA_real_, nrow = max_setting, ncol = length(probs))
  brier.score.mean <- matrix(NA_real_, nrow = max_setting, ncol = length(thresholds))
  quant.score.se <- matrix(NA_real_, nrow = max_setting, ncol = length(probs))
  brier.score.se <- matrix(NA_real_, nrow = max_setting, ncol = length(thresholds))
  crps.mean <- rep(NA_real_, max_setting)
  crps.se <- rep(NA_real_, max_setting)
  mspe.mean <- rep(NA_real_, max_setting)
  mspe.se <- rep(NA_real_, max_setting)
  mape.mean <- rep(NA_real_, max_setting)
  mape.se <- rep(NA_real_, max_setting)
  coverage.mean <- matrix(NA_real_, nrow = max_setting, ncol = length(score_obj$uncertainty_levels))
  coverage.se <- matrix(NA_real_, nrow = max_setting, ncol = length(score_obj$uncertainty_levels))
  coverage.gap <- matrix(NA_real_, nrow = max_setting, ncol = length(score_obj$uncertainty_levels))
  pit.mean <- rep(NA_real_, max_setting)
  pit.mean.se <- rep(NA_real_, max_setting)
  pit.variance <- rep(NA_real_, max_setting)
  pit.variance.se <- rep(NA_real_, max_setting)
  pit.ks <- rep(NA_real_, max_setting)
  pit.ks.se <- rep(NA_real_, max_setting)
  pit.mae <- rep(NA_real_, max_setting)
  pit.mae.se <- rep(NA_real_, max_setting)
  pit.rmse <- rep(NA_real_, max_setting)
  pit.rmse.se <- rep(NA_real_, max_setting)
  pit.bin.share.mean <- matrix(NA_real_, nrow = max_setting, ncol = length(score_obj$pit_breaks) - 1L)
  pit.bin.share.se <- matrix(NA_real_, nrow = max_setting, ncol = length(score_obj$pit_breaks) - 1L)
  pit.bin.count.total <- matrix(NA_real_, nrow = max_setting, ncol = length(score_obj$pit_breaks) - 1L)
  pit.bin.count.mean <- matrix(NA_real_, nrow = max_setting, ncol = length(score_obj$pit_breaks) - 1L)
  brier.split.score.mean <- if (!is.null(score_obj$brier.split.score)) {
    array(
      NA_real_,
      dim = c(max_setting, length(score_obj$brier.split.band_names), length(score_obj$brier.split.target_probs)),
      dimnames = list(NULL, score_obj$brier.split.band_names, as.character(score_obj$brier.split.target_probs))
    )
  } else {
    NULL
  }
  brier.split.score.se <- if (!is.null(score_obj$brier.split.score)) {
    array(
      NA_real_,
      dim = c(max_setting, length(score_obj$brier.split.band_names), length(score_obj$brier.split.target_probs)),
      dimnames = list(NULL, score_obj$brier.split.band_names, as.character(score_obj$brier.split.target_probs))
    )
  } else {
    NULL
  }
  brier.split.n_obs.total <- if (!is.null(score_obj$brier.split.score)) {
    array(
      NA_real_,
      dim = c(max_setting, length(score_obj$brier.split.band_names), length(score_obj$brier.split.target_probs)),
      dimnames = list(NULL, score_obj$brier.split.band_names, as.character(score_obj$brier.split.target_probs))
    )
  } else {
    NULL
  }
  brier.split.n_obs.mean <- if (!is.null(score_obj$brier.split.score)) {
    array(
      NA_real_,
      dim = c(max_setting, length(score_obj$brier.split.band_names), length(score_obj$brier.split.target_probs)),
      dimnames = list(NULL, score_obj$brier.split.band_names, as.character(score_obj$brier.split.target_probs))
    )
  } else {
    NULL
  }
  brier.split.obs.share <- if (!is.null(score_obj$brier.split.score)) {
    array(
      NA_real_,
      dim = c(max_setting, length(score_obj$brier.split.band_names), length(score_obj$brier.split.target_probs)),
      dimnames = list(NULL, score_obj$brier.split.band_names, as.character(score_obj$brier.split.target_probs))
    )
  } else {
    NULL
  }
  classification.metric.mean <- if (!is.null(score_obj$classification.metric)) {
    array(
      NA_real_,
      dim = c(max_setting, length(score_obj$classification.metric_names), length(score_obj$classification.target_probs)),
      dimnames = list(NULL, score_obj$classification.metric_names, as.character(score_obj$classification.target_probs))
    )
  } else {
    NULL
  }
  classification.metric.se <- if (!is.null(score_obj$classification.metric)) {
    array(
      NA_real_,
      dim = c(max_setting, length(score_obj$classification.metric_names), length(score_obj$classification.target_probs)),
      dimnames = list(NULL, score_obj$classification.metric_names, as.character(score_obj$classification.target_probs))
    )
  } else {
    NULL
  }
  classification.count.total <- if (!is.null(score_obj$classification.count)) {
    array(
      NA_real_,
      dim = c(max_setting, length(score_obj$classification.count_names), length(score_obj$classification.target_probs)),
      dimnames = list(NULL, score_obj$classification.count_names, as.character(score_obj$classification.target_probs))
    )
  } else {
    NULL
  }
  classification.count.mean <- if (!is.null(score_obj$classification.count)) {
    array(
      NA_real_,
      dim = c(max_setting, length(score_obj$classification.count_names), length(score_obj$classification.target_probs)),
      dimnames = list(NULL, score_obj$classification.count_names, as.character(score_obj$classification.target_probs))
    )
  } else {
    NULL
  }
  classification.n_obs.total <- if (!is.null(score_obj$classification.count)) {
    matrix(
      NA_real_,
      nrow = max_setting,
      ncol = length(score_obj$classification.target_probs),
      dimnames = list(NULL, as.character(score_obj$classification.target_probs))
    )
  } else {
    NULL
  }
  classification.n_obs.mean <- if (!is.null(score_obj$classification.count)) {
    matrix(
      NA_real_,
      nrow = max_setting,
      ncol = length(score_obj$classification.target_probs),
      dimnames = list(NULL, as.character(score_obj$classification.target_probs))
    )
  } else {
    NULL
  }
  classification.actual_positive_share <- if (!is.null(score_obj$classification.count)) {
    matrix(
      NA_real_,
      nrow = max_setting,
      ncol = length(score_obj$classification.target_probs),
      dimnames = list(NULL, as.character(score_obj$classification.target_probs))
    )
  } else {
    NULL
  }
  classification.predicted_positive_share <- if (!is.null(score_obj$classification.count)) {
    matrix(
      NA_real_,
      nrow = max_setting,
      ncol = length(score_obj$classification.target_probs),
      dimnames = list(NULL, as.character(score_obj$classification.target_probs))
    )
  } else {
    NULL
  }
  classification.draws.mean <- if (!is.null(score_obj$classification.count)) {
    rep(NA_real_, max_setting)
  } else {
    NULL
  }
  summary.draws.mean <- rep(NA_real_, max_setting)
  summary.n_obs.total <- rep(NA_real_, max_setting)

  for (i in available_settings) {
    quant.score.mean[i, ] <- apply(score_obj$quant.score[, , i, drop = FALSE], 1, mean, na.rm = TRUE)
    quant.score.se[i, ] <- apply(score_obj$quant.score[, , i, drop = FALSE], 1, sd, na.rm = TRUE) / sqrt(nsets)
    brier.score.mean[i, ] <- apply(score_obj$brier.score[, , i, drop = FALSE], 1, mean, na.rm = TRUE)
    brier.score.se[i, ] <- apply(score_obj$brier.score[, , i, drop = FALSE], 1, sd, na.rm = TRUE) / sqrt(nsets)
    mspe.mean[i] <- mean_if_any(score_obj$mspe.score[, i])
    mspe.se[i] <- se_if_any(score_obj$mspe.score[, i])
    mape.mean[i] <- mean_if_any(score_obj$mape.score[, i])
    mape.se[i] <- se_if_any(score_obj$mape.score[, i])

    if (!is.null(score_obj$brier.split.score)) {
      for (band_idx in seq_along(score_obj$brier.split.band_names)) {
        for (target_idx in seq_along(score_obj$brier.split.target_probs)) {
          score_vals <- score_obj$brier.split.score[band_idx, target_idx, , i]
          count_vals <- score_obj$brier.split.count[band_idx, target_idx, , i]

          brier.split.score.mean[i, band_idx, target_idx] <- mean_if_any(score_vals)
          brier.split.score.se[i, band_idx, target_idx] <- se_if_any(score_vals)
          brier.split.n_obs.total[i, band_idx, target_idx] <- sum_if_any(count_vals)
          brier.split.n_obs.mean[i, band_idx, target_idx] <- mean_if_any(count_vals)
        }
      }

      for (target_idx in seq_along(score_obj$brier.split.target_probs)) {
        all_total <- brier.split.n_obs.total[i, "all", target_idx]
        if (is.finite(all_total) && all_total > 0) {
          brier.split.obs.share[i, , target_idx] <- brier.split.n_obs.total[i, , target_idx] / all_total
        }
      }
    }

    if (!is.null(score_obj$classification.metric)) {
      for (metric_idx in seq_along(score_obj$classification.metric_names)) {
        for (target_idx in seq_along(score_obj$classification.target_probs)) {
          metric_vals <- score_obj$classification.metric[metric_idx, target_idx, , i]
          classification.metric.mean[i, metric_idx, target_idx] <- mean_if_any(metric_vals)
          classification.metric.se[i, metric_idx, target_idx] <- se_if_any(metric_vals)
        }
      }

      for (count_idx in seq_along(score_obj$classification.count_names)) {
        for (target_idx in seq_along(score_obj$classification.target_probs)) {
          count_vals <- score_obj$classification.count[count_idx, target_idx, , i]
          classification.count.total[i, count_idx, target_idx] <- sum_if_any(count_vals)
          classification.count.mean[i, count_idx, target_idx] <- mean_if_any(count_vals)
        }
      }

      for (target_idx in seq_along(score_obj$classification.target_probs)) {
        obs_vals <- score_obj$classification.obs.total[target_idx, , i]
        actual_share_vals <- score_obj$classification.actual_positive_share[target_idx, , i]
        predicted_share_vals <- score_obj$classification.predicted_positive_share[target_idx, , i]

        classification.n_obs.total[i, target_idx] <- sum_if_any(obs_vals)
        classification.n_obs.mean[i, target_idx] <- mean_if_any(obs_vals)
        classification.actual_positive_share[i, target_idx] <- mean_if_any(actual_share_vals)
        classification.predicted_positive_share[i, target_idx] <- mean_if_any(predicted_share_vals)
      }

      classification.draws.mean[i] <- mean_if_any(score_obj$classification.draws.used[, i])
    }

    if (!is.null(score_obj$crps.score)) {
      crps.mean[i] <- mean_if_any(score_obj$crps.score[, i])
      crps.se[i] <- se_if_any(score_obj$crps.score[, i])

      for (k in seq_along(score_obj$uncertainty_levels)) {
        coverage.mean[i, k] <- mean_if_any(score_obj$coverage.score[k, , i])
        coverage.se[i, k] <- se_if_any(score_obj$coverage.score[k, , i])
        coverage.gap[i, k] <- coverage.mean[i, k] - score_obj$uncertainty_levels[k]
      }

      pit.mean[i] <- mean_if_any(score_obj$pit.mean.score[, i])
      pit.mean.se[i] <- se_if_any(score_obj$pit.mean.score[, i])
      pit.variance[i] <- mean_if_any(score_obj$pit.variance.score[, i])
      pit.variance.se[i] <- se_if_any(score_obj$pit.variance.score[, i])
      pit.ks[i] <- mean_if_any(score_obj$pit.ks.score[, i])
      pit.ks.se[i] <- se_if_any(score_obj$pit.ks.score[, i])
      pit.mae[i] <- mean_if_any(score_obj$pit.mae.score[, i])
      pit.mae.se[i] <- se_if_any(score_obj$pit.mae.score[, i])
      pit.rmse[i] <- mean_if_any(score_obj$pit.rmse.score[, i])
      pit.rmse.se[i] <- se_if_any(score_obj$pit.rmse.score[, i])

      for (b in seq_len(ncol(pit.bin.share.mean))) {
        pit.bin.share.mean[i, b] <- mean_if_any(score_obj$pit.bin.share[b, , i])
        pit.bin.share.se[i, b] <- se_if_any(score_obj$pit.bin.share[b, , i])
        pit.bin.count.total[i, b] <- sum_if_any(score_obj$pit.bin.count[b, , i])
        pit.bin.count.mean[i, b] <- mean_if_any(score_obj$pit.bin.count[b, , i])
      }

      summary.draws.mean[i] <- mean_if_any(score_obj$summary.draws.used[, i])
      summary.n_obs.total[i] <- sum_if_any(score_obj$summary.n_obs[, i])
    }
  }

  bs.mean.ref.gau <- matrix(NA_real_, nrow = max_setting, ncol = length(thresholds))
  qs.mean.ref.gau <- matrix(NA_real_, nrow = max_setting, ncol = length(probs))
  crps.mean.ref.gau <- rep(NA_real_, max_setting)
  mspe.mean.ref.gau <- rep(NA_real_, max_setting)
  mape.mean.ref.gau <- rep(NA_real_, max_setting)
  brier.split.rel.ref.gau <- if (!is.null(score_obj$brier.split.score)) {
    array(
      NA_real_,
      dim = c(max_setting, length(score_obj$brier.split.band_names), length(score_obj$brier.split.target_probs)),
      dimnames = list(NULL, score_obj$brier.split.band_names, as.character(score_obj$brier.split.target_probs))
    )
  } else {
    NULL
  }
  classification.metric.rel.ref.gau <- if (!is.null(score_obj$classification.metric)) {
    array(
      NA_real_,
      dim = c(max_setting, length(score_obj$classification.metric_names), length(score_obj$classification.target_probs)),
      dimnames = list(NULL, score_obj$classification.metric_names, as.character(score_obj$classification.target_probs))
    )
  } else {
    NULL
  }
  classification.metric.delta.ref.gau <- if (!is.null(score_obj$classification.metric)) {
    array(
      NA_real_,
      dim = c(max_setting, length(score_obj$classification.metric_names), length(score_obj$classification.target_probs)),
      dimnames = list(NULL, score_obj$classification.metric_names, as.character(score_obj$classification.target_probs))
    )
  } else {
    NULL
  }

  for (i in available_settings) {
    bs.mean.ref.gau[i, ] <- brier.score.mean[i, ] / brier.score.mean[baseline_setting, ]
    qs.mean.ref.gau[i, ] <- quant.score.mean[i, ] / quant.score.mean[baseline_setting, ]
    if (is.finite(crps.mean[i]) && is.finite(crps.mean[baseline_setting]) && crps.mean[baseline_setting] != 0) {
      crps.mean.ref.gau[i] <- crps.mean[i] / crps.mean[baseline_setting]
    }
    if (is.finite(mspe.mean[i]) && is.finite(mspe.mean[baseline_setting]) && mspe.mean[baseline_setting] != 0) {
      mspe.mean.ref.gau[i] <- mspe.mean[i] / mspe.mean[baseline_setting]
    }
    if (is.finite(mape.mean[i]) && is.finite(mape.mean[baseline_setting]) && mape.mean[baseline_setting] != 0) {
      mape.mean.ref.gau[i] <- mape.mean[i] / mape.mean[baseline_setting]
    }

    if (!is.null(brier.split.rel.ref.gau)) {
      for (band_idx in seq_along(score_obj$brier.split.band_names)) {
        for (target_idx in seq_along(score_obj$brier.split.target_probs)) {
          base_score <- brier.split.score.mean[baseline_setting, band_idx, target_idx]
          current_score <- brier.split.score.mean[i, band_idx, target_idx]
          if (is.finite(base_score) && base_score != 0 && is.finite(current_score)) {
            brier.split.rel.ref.gau[i, band_idx, target_idx] <- current_score / base_score
          }
        }
      }
    }

    if (!is.null(classification.metric.rel.ref.gau)) {
      for (metric_idx in seq_along(score_obj$classification.metric_names)) {
        for (target_idx in seq_along(score_obj$classification.target_probs)) {
          base_metric <- classification.metric.mean[baseline_setting, metric_idx, target_idx]
          current_metric <- classification.metric.mean[i, metric_idx, target_idx]
          if (is.finite(base_metric) && base_metric != 0 && is.finite(current_metric)) {
            classification.metric.rel.ref.gau[i, metric_idx, target_idx] <- current_metric / base_metric
          }
          if (is.finite(base_metric) && is.finite(current_metric)) {
            classification.metric.delta.ref.gau[i, metric_idx, target_idx] <- current_metric - base_metric
          }
        }
      }
    }
  }

  list(
    baseline_setting = baseline_setting,
    max_setting = max_setting,
    nsets = nsets,
    probs = probs,
    thresholds = thresholds,
    threshold_probs = threshold_probs,
    available_settings = available_settings,
    quant.score = score_obj$quant.score,
    brier.score = score_obj$brier.score,
    mspe.score = score_obj$mspe.score,
    mape.score = score_obj$mape.score,
    brier.split.target_probs = score_obj$brier.split.target_probs,
    brier.split.target_thresholds = score_obj$brier.split.target_thresholds,
    brier.split.band_names = score_obj$brier.split.band_names,
    brier.split.score.mean = brier.split.score.mean,
    brier.split.score.se = brier.split.score.se,
    brier.split.n_obs.total = brier.split.n_obs.total,
    brier.split.n_obs.mean = brier.split.n_obs.mean,
    brier.split.obs.share = brier.split.obs.share,
    brier.split.rel.ref.gau = brier.split.rel.ref.gau,
    classification.target_probs = score_obj$classification.target_probs,
    classification.target_thresholds = score_obj$classification.target_thresholds,
    classification.metric_names = score_obj$classification.metric_names,
    classification.count_names = score_obj$classification.count_names,
    classification.metric.mean = classification.metric.mean,
    classification.metric.se = classification.metric.se,
    classification.count.total = classification.count.total,
    classification.count.mean = classification.count.mean,
    classification.n_obs.total = classification.n_obs.total,
    classification.n_obs.mean = classification.n_obs.mean,
    classification.actual_positive_share = classification.actual_positive_share,
    classification.predicted_positive_share = classification.predicted_positive_share,
    classification.draws.mean = classification.draws.mean,
    classification.metric.rel.ref.gau = classification.metric.rel.ref.gau,
    classification.metric.delta.ref.gau = classification.metric.delta.ref.gau,
    classification.probability_cutoff = score_obj$classification.probability_cutoff,
    quant.score.mean = quant.score.mean,
    brier.score.mean = brier.score.mean,
    quant.score.se = quant.score.se,
    brier.score.se = brier.score.se,
    crps.mean = crps.mean,
    crps.se = crps.se,
    crps.mean.ref.gau = crps.mean.ref.gau,
    mspe.mean = mspe.mean,
    mspe.se = mspe.se,
    mspe.mean.ref.gau = mspe.mean.ref.gau,
    mape.mean = mape.mean,
    mape.se = mape.se,
    mape.mean.ref.gau = mape.mean.ref.gau,
    uncertainty_levels = score_obj$uncertainty_levels,
    coverage.mean = coverage.mean,
    coverage.se = coverage.se,
    coverage.gap = coverage.gap,
    pit_breaks = score_obj$pit_breaks,
    pit.expected_share = diff(score_obj$pit_breaks),
    pit.expected_mean = 0.5,
    pit.expected_variance = 1 / 12,
    pit.mean = pit.mean,
    pit.mean.se = pit.mean.se,
    pit.variance = pit.variance,
    pit.variance.se = pit.variance.se,
    pit.ks = pit.ks,
    pit.ks.se = pit.ks.se,
    pit.mae = pit.mae,
    pit.mae.se = pit.mae.se,
    pit.rmse = pit.rmse,
    pit.rmse.se = pit.rmse.se,
    pit.bin.share.mean = pit.bin.share.mean,
    pit.bin.share.se = pit.bin.share.se,
    pit.bin.count.total = pit.bin.count.total,
    pit.bin.count.mean = pit.bin.count.mean,
    summary.draws.mean = summary.draws.mean,
    summary.n_obs.total = summary.n_obs.total,
    predictive_density_status = rep("unavailable", max_setting),
    predictive_density_note = rep(
      "Pointwise log-likelihood is not stored in fit[[d]]; LOO-ELPD and WAIC are exported as NA placeholders.",
      max_setting
    ),
    bs.mean.ref.gau = bs.mean.ref.gau,
    qs.mean.ref.gau = qs.mean.ref.gau,
    skipped_missing_file = score_obj$skipped_missing_file,
    skipped_bad_contract = score_obj$skipped_bad_contract,
    skipped_scoring_error = score_obj$skipped_scoring_error
  )
}

print_score_summary <- function(score_obj, label = "US-all scoring") {
  cat("\n===", label, "===\n")
  cat("Requested settings:", length(score_obj$setting_ids), "\n")
  cat("Available settings:", length(score_obj$available_settings), "\n")
  cat("Missing result files:", length(score_obj$skipped_missing_file), "\n")
  if (length(score_obj$skipped_missing_file) > 0) {
    cat("  ->", paste(score_obj$skipped_missing_file, collapse = ", "), "\n")
  }
  cat("Contract mismatch settings:", length(score_obj$skipped_bad_contract), "\n")
  if (length(score_obj$skipped_bad_contract) > 0) {
    cat("  ->", paste(score_obj$skipped_bad_contract, collapse = ", "), "\n")
  }
  cat("Scoring error settings:", length(score_obj$skipped_scoring_error), "\n")
  if (length(score_obj$skipped_scoring_error) > 0) {
    cat("  ->", paste(score_obj$skipped_scoring_error, collapse = ", "), "\n")
  }
}
