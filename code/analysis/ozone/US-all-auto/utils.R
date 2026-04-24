# utils.R — shared helpers for the US-all-auto pipeline
# Source this file; do not run it directly.

# ---- String / env helpers ----

truthy <- function(x) tolower(trimws(as.character(x))) %in% c("1", "true", "yes", "y")

parse_optional_int <- function(x) {
  txt <- trimws(as.character(x[1]))
  if (!nzchar(txt)) return(NA_integer_)
  val <- suppressWarnings(as.numeric(txt))
  if (!is.finite(val)) return(NA_integer_)
  as.integer(round(val))
}

parse_env_int <- function(name, default) {
  val <- suppressWarnings(as.integer(trimws(Sys.getenv(name, unset = ""))))
  if (length(val) == 1L && is.finite(val) && val > 0L) val else as.integer(default)
}

expand_range_tokens <- function(raw_str) {
  tokens <- unlist(strsplit(raw_str, ",", fixed = TRUE))
  out <- integer(0L)
  for (tok in tokens) {
    tok <- trimws(tok)
    if (!nzchar(tok)) next
    if (grepl("^[0-9]+:[0-9]+$", tok)) {
      parts <- as.integer(strsplit(tok, ":", fixed = TRUE)[[1]])
      if (length(parts) == 2L && !anyNA(parts)) {
        step <- if (parts[1L] <= parts[2L]) 1L else -1L
        out  <- c(out, seq.int(parts[1L], parts[2L], by = step))
      }
    } else {
      val <- suppressWarnings(as.integer(tok))
      if (!is.na(val)) out <- c(out, val)
    }
  }
  sort(unique(out))
}

# ---- Model spec helpers ----

resolve_model <- function(method_raw) {
  switch(tolower(trimws(method_raw)),
    "skew-t"     = list(method = "t",          skew = TRUE,  is_maxstable = FALSE),
    "t"          = list(method = "t",          skew = FALSE, is_maxstable = FALSE),
    "gaussian"   = list(method = "gaussian",   skew = FALSE, is_maxstable = FALSE),
    "max-stable" = list(method = "max-stable", skew = FALSE, is_maxstable = TRUE),
    stop(sprintf("Unsupported method: %s", method_raw), call. = FALSE)
  )
}

pick_tau_init <- function(method_name, tau_init) {
  if (tolower(method_name) == "gaussian") tau_init else 0.05
}

# ---- MRTS covariate builder ----

build_mrts_covariates <- function(S_train, S_pred, nt, k) {
  if (is.na(k) || k <= 0) stop("MRTS k must be a positive integer.", call. = FALSE)

  safe_mrts_call <- function(fn, label) {
    train_obj <- tryCatch(fn(S_train, k = k), error = function(e) NULL)
    pred_obj  <- tryCatch(fn(S_train, k = k, x = S_pred), error = function(e) NULL)
    if (!is.null(train_obj) && !is.null(pred_obj))
      return(list(train = as.matrix(train_obj), pred = as.matrix(pred_obj), source = label))
    all_obj <- tryCatch(fn(rbind(S_train, S_pred), k = k), error = function(e) NULL)
    if (is.null(all_obj)) return(NULL)
    all_mat <- as.matrix(all_obj)
    if (nrow(all_mat) != nrow(S_train) + nrow(S_pred)) return(NULL)
    n_tr <- nrow(S_train)
    list(train = all_mat[seq_len(n_tr), , drop = FALSE],
         pred  = all_mat[(n_tr + 1):nrow(all_mat), , drop = FALSE],
         source = paste0(label, " [combined fallback]"))
  }

  basis_obj <- NULL
  if (requireNamespace("autoFRK", quietly = TRUE) &&
      exists("mrts", where = asNamespace("autoFRK"), inherits = FALSE))
    basis_obj <- safe_mrts_call(autoFRK::mrts, "autoFRK::mrts")
  if (is.null(basis_obj) && exists("mrts", mode = "function"))
    basis_obj <- safe_mrts_call(mrts, "mrts()")
  if (is.null(basis_obj))
    stop("Unable to build MRTS covariates (mrts unavailable or failed).", call. = FALSE)

  col_sd <- apply(basis_obj$train, 2, sd)
  keep   <- which(is.finite(col_sd) & col_sd > 1e-10)
  if (length(keep) == 0L)
    stop("All MRTS columns are near-constant after construction.", call. = FALSE)

  basis_train <- basis_obj$train[, keep, drop = FALSE]
  basis_pred  <- basis_obj$pred[, keep, drop = FALSE]

  mrts_train <- array(0, dim = c(nrow(S_train), nt, ncol(basis_train)))
  mrts_pred  <- array(0, dim = c(nrow(S_pred),  nt, ncol(basis_pred)))
  for (t in seq_len(nt)) {
    mrts_train[, t, ] <- basis_train
    mrts_pred[, t, ]  <- basis_pred
  }

  list(train = mrts_train, pred = mrts_pred, source = basis_obj$source,
       original_cols = ncol(basis_obj$train),
       kept_cols     = ncol(basis_train),
       dropped_cols  = ncol(basis_obj$train) - ncol(basis_train))
}
