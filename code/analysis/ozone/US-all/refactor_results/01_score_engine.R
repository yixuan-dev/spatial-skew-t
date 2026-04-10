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

compute_us_all_scores <- function(
    setting_ids,
    result_path_fn,
    Y,
    cv_lst,
    probs,
    thresholds,
    trans_setting_ids = 2L,
    enforce_contract = TRUE
) {
  setting_ids <- sort(unique(setting_ids[!is.na(setting_ids)]))
  if (length(setting_ids) == 0) {
    stop("No valid setting IDs were provided.")
  }

  nsets <- length(cv_lst)
  max_setting <- max(setting_ids)

  quant.score <- array(NA_real_, dim = c(length(probs), nsets, max_setting))
  brier.score <- array(NA_real_, dim = c(length(thresholds), nsets, max_setting))

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

    for (d in seq_len(nsets)) {
      fit.d <- fit[[d]]
      if (is.null(fit.d) || is.null(fit.d$yp)) {
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
          list(
            quant = QuantScore(pred.d, probs, validate, trans = trans),
            brier = BrierScore(pred.d, thresholds, validate, trans = trans)
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
    quant.score = quant.score,
    brier.score = brier.score,
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

  for (i in available_settings) {
    quant.score.mean[i, ] <- apply(score_obj$quant.score[, , i, drop = FALSE], 1, mean, na.rm = TRUE)
    quant.score.se[i, ] <- apply(score_obj$quant.score[, , i, drop = FALSE], 1, sd, na.rm = TRUE) / sqrt(nsets)
    brier.score.mean[i, ] <- apply(score_obj$brier.score[, , i, drop = FALSE], 1, mean, na.rm = TRUE)
    brier.score.se[i, ] <- apply(score_obj$brier.score[, , i, drop = FALSE], 1, sd, na.rm = TRUE) / sqrt(nsets)
  }

  bs.mean.ref.gau <- matrix(NA_real_, nrow = max_setting, ncol = length(thresholds))
  qs.mean.ref.gau <- matrix(NA_real_, nrow = max_setting, ncol = length(probs))

  for (i in available_settings) {
    bs.mean.ref.gau[i, ] <- brier.score.mean[i, ] / brier.score.mean[baseline_setting, ]
    qs.mean.ref.gau[i, ] <- quant.score.mean[i, ] / quant.score.mean[baseline_setting, ]
  }

  list(
    baseline_setting = baseline_setting,
    max_setting = max_setting,
    nsets = nsets,
    probs = probs,
    thresholds = thresholds,
    available_settings = available_settings,
    quant.score = score_obj$quant.score,
    brier.score = score_obj$brier.score,
    quant.score.mean = quant.score.mean,
    brier.score.mean = brier.score.mean,
    quant.score.se = quant.score.se,
    brier.score.se = brier.score.se,
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
