prop_normalize_coords <- function(coords, mins = NULL, maxs = NULL) {
  coords <- as.matrix(coords)
  if (is.null(mins)) {
    mins <- apply(coords, 2, min)
  }
  if (is.null(maxs)) {
    maxs <- apply(coords, 2, max)
  }
  spans <- pmax(maxs - mins, 1e-8)
  sweep(sweep(coords, 2, mins, "-"), 2, spans, "/")
}

prop_try_mrts_basis <- function(s_obs, s_pred = NULL, k) {
  safe_mrts_call <- function(fn, label) {
    obs_obj <- tryCatch(fn(s_obs, k = k), error = function(e) NULL)
    pred_obj <- NULL
    if (!is.null(s_pred)) {
      pred_obj <- tryCatch(fn(s_obs, k = k, x = s_pred), error = function(e) NULL)
    }

    if (!is.null(obs_obj) && (is.null(s_pred) || !is.null(pred_obj))) {
      return(list(
        obs = as.matrix(obs_obj),
        pred = if (is.null(s_pred)) NULL else as.matrix(pred_obj),
        source = label
      ))
    }

    combined <- if (is.null(s_pred)) s_obs else rbind(s_obs, s_pred)
    all_obj <- tryCatch(fn(combined, k = k), error = function(e) NULL)
    if (is.null(all_obj)) {
      return(NULL)
    }

    all_mat <- as.matrix(all_obj)
    if (nrow(all_mat) != nrow(combined)) {
      return(NULL)
    }

    n_obs <- nrow(s_obs)
    list(
      obs = all_mat[seq_len(n_obs), , drop = FALSE],
      pred = if (is.null(s_pred)) NULL else all_mat[(n_obs + 1L):nrow(all_mat), , drop = FALSE],
      source = paste0(label, " [combined fallback]")
    )
  }

  basis_obj <- NULL
  if (requireNamespace("autoFRK", quietly = TRUE) &&
      exists("mrts", where = asNamespace("autoFRK"), inherits = FALSE)) {
    basis_obj <- safe_mrts_call(autoFRK::mrts, "autoFRK::mrts")
  }

  if (is.null(basis_obj) && exists("mrts", mode = "function")) {
    basis_obj <- safe_mrts_call(mrts, "mrts()")
  }

  basis_obj
}

prop_fallback_basis <- function(s_obs, s_pred = NULL, k) {
  coords_all <- if (is.null(s_pred)) as.matrix(s_obs) else rbind(as.matrix(s_obs), as.matrix(s_pred))
  coords_scaled <- prop_normalize_coords(coords_all)

  centers <- matrix(NA_real_, nrow = 0L, ncol = ncol(coords_scaled))
  scales <- numeric(0L)
  level <- 0L
  while (nrow(centers) < k) {
    side <- 2 ^ level
    axes <- replicate(ncol(coords_scaled), (seq_len(side) - 0.5) / side, simplify = FALSE)
    level_centers <- as.matrix(do.call(expand.grid, axes))
    centers <- rbind(centers, level_centers)
    scales <- c(scales, rep(max(0.12, 0.9 / side), nrow(level_centers)))
    level <- level + 1L
  }

  centers <- centers[seq_len(k), , drop = FALSE]
  scales <- scales[seq_len(k)]
  dmat <- fields::rdist(coords_scaled, centers)
  basis_all <- exp(-0.5 * sweep(dmat, 2, scales, "/") ^ 2)

  n_obs <- nrow(s_obs)
  list(
    obs = basis_all[seq_len(n_obs), , drop = FALSE],
    pred = if (is.null(s_pred)) NULL else basis_all[(n_obs + 1L):nrow(basis_all), , drop = FALSE],
    source = "radial multiresolution fallback"
  )
}

prop_trim_basis <- function(obs_basis, pred_basis = NULL, tol = 1e-10) {
  col_sd <- apply(obs_basis, 2, sd)
  keep <- which(is.finite(col_sd) & col_sd > tol)
  if (length(keep) == 0L) {
    stop("No non-degenerate basis columns remain after trimming.", call. = FALSE)
  }

  list(
    obs = obs_basis[, keep, drop = FALSE],
    pred = if (is.null(pred_basis)) NULL else pred_basis[, keep, drop = FALSE],
    original_cols = ncol(obs_basis),
    kept_cols = length(keep)
  )
}

prop_orthonormalize_basis <- function(obs_basis, pred_basis = NULL, tol = 1e-8) {
  gram <- prop_symmetrize(crossprod(obs_basis))
  eig <- eigen(gram, symmetric = TRUE)
  positive <- eig$values > (tol * max(eig$values, 1))
  if (!any(positive)) {
    stop("Basis Gram matrix is numerically singular.", call. = FALSE)
  }

  vecs <- eig$vectors[, positive, drop = FALSE]
  vals <- eig$values[positive]
  transform <- vecs %*% diag(1 / sqrt(vals), nrow = length(vals))

  list(
    q_obs = obs_basis %*% transform,
    q_pred = if (is.null(pred_basis)) NULL else pred_basis %*% transform,
    rank_kept = length(vals),
    eigenvalues = vals
  )
}

build_prop_basis <- function(s_obs, s_pred = NULL, k, tol = 1e-8) {
  if (is.na(k) || k < 1L) {
    stop("prop_k must be a positive integer.", call. = FALSE)
  }

  basis_obj <- prop_try_mrts_basis(s_obs = s_obs, s_pred = s_pred, k = k)
  if (is.null(basis_obj)) {
    basis_obj <- prop_fallback_basis(s_obs = s_obs, s_pred = s_pred, k = k)
  }

  trimmed <- prop_trim_basis(
    obs_basis = as.matrix(basis_obj$obs),
    pred_basis = if (is.null(basis_obj$pred)) NULL else as.matrix(basis_obj$pred)
  )
  ortho <- prop_orthonormalize_basis(trimmed$obs, trimmed$pred, tol = tol)

  list(
    q_obs = ortho$q_obs,
    q_pred = ortho$q_pred,
    source = basis_obj$source,
    requested_k = as.integer(k),
    original_cols = trimmed$original_cols,
    kept_cols = trimmed$kept_cols,
    rank_kept = ortho$rank_kept
  )
}
