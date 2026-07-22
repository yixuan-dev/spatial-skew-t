# MRTS covariate construction shared by us-all-run.R (CV folds) and the
# full-data map pipeline (us-all-full-204.R / predict-cincy-204.R).
# Moved from us-all-run.R (unchanged except the added keep_idx return field)
# so both callers use one definition of the basis and of the near-constant
# column-drop rule (sd <= 1e-10). The kept
# columns depend only on (S_train, k); S_pred is just evaluation points, so a
# fit script and a later prediction script that pass the same S_train recover
# the same columns.

build_mrts_covariates <- function(S_train, S_pred, nt, k) {
  if (is.na(k) || k <= 0) {
    stop("MRTS k must be a positive integer.", call. = FALSE)
  }

  safe_mrts_call <- function(fn, label) {
    train_obj <- tryCatch(fn(S_train, k = k), error = function(e) NULL)
    pred_obj <- tryCatch(fn(S_train, k = k, x = S_pred), error = function(e) NULL)

    if (!is.null(train_obj) && !is.null(pred_obj)) {
      train_mat <- as.matrix(train_obj)
      pred_mat <- as.matrix(pred_obj)
      return(list(train = train_mat, pred = pred_mat, source = label))
    }

    all_obj <- tryCatch(fn(rbind(S_train, S_pred), k = k), error = function(e) NULL)
    if (is.null(all_obj)) {
      return(NULL)
    }

    all_mat <- as.matrix(all_obj)
    if (nrow(all_mat) != nrow(S_train) + nrow(S_pred)) {
      return(NULL)
    }

    n_train <- nrow(S_train)
    list(
      train = all_mat[seq_len(n_train), , drop = FALSE],
      pred = all_mat[(n_train + 1):nrow(all_mat), , drop = FALSE],
      source = paste0(label, " [combined fallback]")
    )
  }

  basis_obj <- NULL
  if (requireNamespace("autoFRK", quietly = TRUE) && exists("mrts", where = asNamespace("autoFRK"), inherits = FALSE)) {
    basis_obj <- safe_mrts_call(autoFRK::mrts, "autoFRK::mrts")
  }

  if (is.null(basis_obj) && exists("mrts", mode = "function")) {
    basis_obj <- safe_mrts_call(mrts, "mrts()")
  }

  if (is.null(basis_obj)) {
    stop("Unable to build MRTS covariates for this setting (mrts function unavailable or failed).", call. = FALSE)
  }

  basis_train <- as.matrix(basis_obj$train)
  basis_pred <- as.matrix(basis_obj$pred)

  if (nrow(basis_train) != nrow(S_train) || nrow(basis_pred) != nrow(S_pred)) {
    stop("MRTS basis row count does not match train/pred locations.", call. = FALSE)
  }

  if (ncol(basis_train) != ncol(basis_pred)) {
    stop("MRTS basis column mismatch between train and pred matrices.", call. = FALSE)
  }

  original_cols <- ncol(basis_train)
  col_sd <- apply(basis_train, 2, sd)
  keep <- which(is.finite(col_sd) & col_sd > 1e-10)

  if (length(keep) == 0) {
    stop("All MRTS columns are near-constant after construction.", call. = FALSE)
  }

  basis_train_kept <- basis_train[, keep, drop = FALSE]
  basis_pred_kept <- basis_pred[, keep, drop = FALSE]

  mrts_train <- array(0, dim = c(nrow(S_train), nt, ncol(basis_train_kept)))
  mrts_pred <- array(0, dim = c(nrow(S_pred), nt, ncol(basis_pred_kept)))

  for (t in seq_len(nt)) {
    mrts_train[, t, ] <- basis_train_kept
    mrts_pred[, t, ] <- basis_pred_kept
  }

  list(
    train = mrts_train,
    pred = mrts_pred,
    source = basis_obj$source,
    original_cols = original_cols,
    kept_cols = ncol(basis_train_kept),
    dropped_cols = original_cols - ncol(basis_train_kept),
    keep_idx = keep
  )
}
