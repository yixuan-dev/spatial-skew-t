build_basis_matrix <- function(s_obs, s_pred = NULL, k, tol = 1e-8) {
  s_obs <- prop_check_numeric_matrix(s_obs, "s_obs")
  if (nrow(s_obs) < 2L) {
    stop("s_obs must contain at least two observation locations.", call. = FALSE)
  }
  if (ncol(s_obs) < 1L) {
    stop("s_obs must contain at least one coordinate column.", call. = FALSE)
  }
  if (!is.null(s_pred)) {
    s_pred <- prop_check_numeric_matrix(
      s_pred,
      "s_pred",
      ncol_expected = ncol(s_obs)
    )
  }
  k <- prop_check_scalar_numeric(k, "k", lower = 1, integerish = TRUE)
  tol <- prop_check_scalar_numeric(tol, "tol", lower = 0, strict_lower = TRUE)

  basis_obj <- build_prop_basis(
    s_obs = s_obs,
    s_pred = s_pred,
    k = k,
    tol = tol
  )
  basis_obj <- prop_require_names(
    basis_obj,
    "basis_obj",
    c("q_obs", "q_pred", "requested_k", "original_cols", "kept_cols", "rank_kept", "source")
  )

  list(
    F_obs = basis_obj$q_obs,
    F_pred = basis_obj$q_pred,
    requested_k = basis_obj$requested_k,
    original_cols = basis_obj$original_cols,
    kept_cols = basis_obj$kept_cols,
    rank_kept = basis_obj$rank_kept,
    source = basis_obj$source
  )
}

update_residuals <- function(y, x_beta, lambda = 0, zg = NULL, taug = NULL) {
  y <- prop_check_numeric_matrix(y, "y")
  x_beta <- prop_check_numeric_matrix(
    x_beta,
    "x_beta",
    nrow_expected = nrow(y),
    ncol_expected = ncol(y)
  )
  lambda <- prop_check_scalar_numeric(lambda, "lambda")

  if (is.null(zg)) {
    zg <- matrix(0, nrow(y), ncol(y))
  } else {
    zg <- prop_check_numeric_matrix(
      zg,
      "zg",
      nrow_expected = nrow(y),
      ncol_expected = ncol(y)
    )
  }

  mu <- x_beta + lambda * zg
  residuals_raw <- y - mu
  residuals_std <- if (is.null(taug)) {
    NULL
  } else {
    taug <- prop_check_numeric_matrix(
      taug,
      "taug",
      nrow_expected = nrow(y),
      ncol_expected = ncol(y),
      positive = TRUE
    )
    sqrt(taug) * residuals_raw
  }

  list(
    mu = mu,
    raw = residuals_raw,
    std = residuals_std
  )
}

update_M_closed_form <- function(F_obs, residuals_std, F_pred = NULL,
                                 sigma_eps_sq = 0, sigma_floor = 1e-6) {
  F_obs <- prop_check_numeric_matrix(F_obs, "F_obs")
  if (ncol(F_obs) < 1L) {
    stop("F_obs must contain at least one basis column.", call. = FALSE)
  }
  if (ncol(F_obs) > nrow(F_obs)) {
    stop("F_obs must not have more columns than rows.", call. = FALSE)
  }
  residuals_std <- prop_check_numeric_matrix(
    residuals_std,
    "residuals_std",
    nrow_expected = nrow(F_obs)
  )
  if (ncol(residuals_std) < 1L) {
    stop("residuals_std must contain at least one replicate.", call. = FALSE)
  }
  if (!is.null(F_pred)) {
    F_pred <- prop_check_numeric_matrix(
      F_pred,
      "F_pred",
      ncol_expected = ncol(F_obs)
    )
  }
  sigma_eps_sq <- prop_check_scalar_numeric(sigma_eps_sq, "sigma_eps_sq", lower = 0)
  sigma_floor <- prop_check_scalar_numeric(sigma_floor, "sigma_floor", lower = 0)

  update_obj <- prop_closed_form_update(
    residuals = residuals_std,
    q_obs = F_obs,
    sigma_eps_sq = sigma_eps_sq,
    sigma_floor = sigma_floor
  )
  update_obj <- prop_require_names(
    update_obj,
    "update_obj",
    c("rotation", "sigma_xi_sq", "lowrank_eigs", "projected_eigs", "effective_rank", "trace_s", "lstar")
  )

  rotation <- if (is.null(update_obj$rotation)) diag(ncol(F_obs)) else update_obj$rotation
  update_obj$M_hat <- prop_symmetrize(
    rotation %*%
      diag(update_obj$lowrank_eigs, nrow = length(update_obj$lowrank_eigs)) %*%
      t(rotation)
  )
  update_obj$F_rot_obs <- as.matrix(F_obs) %*% rotation
  update_obj$F_rot_pred <- if (is.null(F_pred)) NULL else as.matrix(F_pred) %*% rotation
  update_obj
}

build_R_from_G <- function(F_obs, M_update, F_pred = NULL) {
  F_obs <- prop_check_numeric_matrix(F_obs, "F_obs")
  if (!is.null(F_pred)) {
    F_pred <- prop_check_numeric_matrix(
      F_pred,
      "F_pred",
      ncol_expected = ncol(F_obs)
    )
  }
  M_update <- prop_require_names(
    M_update,
    "M_update",
    c("rotation", "lowrank_eigs", "sigma_xi_sq")
  )
  M_update$rotation <- prop_check_numeric_matrix(
    M_update$rotation,
    "M_update$rotation",
    nrow_expected = ncol(F_obs),
    ncol_expected = ncol(F_obs)
  )
  M_update$lowrank_eigs <- prop_check_numeric_vector(
    M_update$lowrank_eigs,
    "M_update$lowrank_eigs",
    length_expected = ncol(F_obs)
  )
  if (any(M_update$lowrank_eigs < 0)) {
    stop("M_update$lowrank_eigs must be nonnegative.", call. = FALSE)
  }
  M_update$sigma_xi_sq <- prop_check_scalar_numeric(
    M_update$sigma_xi_sq,
    "M_update$sigma_xi_sq",
    lower = 0,
    strict_lower = TRUE
  )

  state <- prop_make_cov_state(
    q_obs = F_obs,
    q_pred = F_pred,
    update_obj = M_update
  )
  prop_prepare_prediction_state(state)
}

quadform_logdet <- function(cov_state, residuals_std) {
  cov_state <- prop_require_names(
    cov_state,
    "cov_state",
    c("prec_obs", "logdet_r_obs", "logdet_prec_obs")
  )
  cov_state$prec_obs <- prop_check_numeric_matrix(cov_state$prec_obs, "cov_state$prec_obs")
  if (nrow(cov_state$prec_obs) != ncol(cov_state$prec_obs)) {
    stop("cov_state$prec_obs must be square.", call. = FALSE)
  }
  cov_state$logdet_r_obs <- prop_check_scalar_numeric(cov_state$logdet_r_obs, "cov_state$logdet_r_obs")
  cov_state$logdet_prec_obs <- prop_check_scalar_numeric(cov_state$logdet_prec_obs, "cov_state$logdet_prec_obs")
  residuals_std <- prop_check_numeric_matrix(
    residuals_std,
    "residuals_std",
    nrow_expected = nrow(cov_state$prec_obs)
  )

  quad_sum <- 0
  for (tt in seq_len(ncol(residuals_std))) {
    quad_sum <- quad_sum + emulator::quad.form(cov_state$prec_obs, residuals_std[, tt])
  }

  list(
    quad_sum = quad_sum,
    logdet_R = cov_state$logdet_r_obs,
    logdet_prec = cov_state$logdet_prec_obs,
    loglik = ncol(residuals_std) * cov_state$logdet_prec_obs - 0.5 * quad_sum
  )
}

predict_at_new_sites <- function(cov_state, residuals_raw, taug, x_pred, beta,
                                 lambda = 0, z_pred = NULL, sigma_pred = NULL,
                                 draw = TRUE) {
  cov_state <- prop_require_names(
    cov_state,
    "cov_state",
    c("q_pred", "r_cross_pred_obs", "prec_obs", "pred_cond_chol")
  )
  if (is.null(cov_state$q_pred)) {
    stop("cov_state does not contain prediction basis rows.", call. = FALSE)
  }
  cov_state$q_pred <- prop_check_numeric_matrix(cov_state$q_pred, "cov_state$q_pred")
  cov_state$r_cross_pred_obs <- prop_check_numeric_matrix(
    cov_state$r_cross_pred_obs,
    "cov_state$r_cross_pred_obs",
    nrow_expected = nrow(cov_state$q_pred)
  )
  cov_state$prec_obs <- prop_check_numeric_matrix(cov_state$prec_obs, "cov_state$prec_obs")
  cov_state$pred_cond_chol <- prop_check_numeric_matrix(
    cov_state$pred_cond_chol,
    "cov_state$pred_cond_chol",
    nrow_expected = nrow(cov_state$q_pred),
    ncol_expected = nrow(cov_state$q_pred)
  )

  residuals_raw <- prop_check_numeric_matrix(residuals_raw, "residuals_raw")
  taug <- prop_check_numeric_matrix(
    taug,
    "taug",
    nrow_expected = nrow(residuals_raw),
    ncol_expected = ncol(residuals_raw),
    positive = TRUE
  )
  x_pred <- prop_check_numeric_array3(
    x_pred,
    "x_pred",
    dim1_expected = nrow(cov_state$q_pred),
    dim2_expected = ncol(residuals_raw)
  )
  beta <- prop_check_numeric_vector(beta, "beta", length_expected = dim(x_pred)[3])
  lambda <- prop_check_scalar_numeric(lambda, "lambda")
  draw <- prop_check_flag(draw, "draw")

  nt <- ncol(residuals_raw)
  np <- nrow(x_pred)
  if (ncol(cov_state$r_cross_pred_obs) != nrow(residuals_raw)) {
    stop("cov_state$r_cross_pred_obs must align with the observed residual dimension.", call. = FALSE)
  }

  if (is.null(z_pred)) {
    z_pred <- matrix(0, np, nt)
  } else {
    z_pred <- prop_check_numeric_matrix(
      z_pred,
      "z_pred",
      nrow_expected = np,
      ncol_expected = nt
    )
  }

  if (is.null(sigma_pred)) {
    sigma_pred <- matrix(1, np, nt)
  } else if (length(sigma_pred) == 1L) {
    sigma_scalar <- prop_check_scalar_numeric(
      sigma_pred,
      "sigma_pred",
      lower = 0,
      strict_lower = TRUE
    )
    sigma_pred <- matrix(sigma_scalar, np, nt)
  } else {
    sigma_pred <- prop_check_numeric_matrix(
      sigma_pred,
      "sigma_pred",
      nrow_expected = np,
      ncol_expected = nt,
      positive = TRUE
    )
  }

  pred_draws <- matrix(NA_real_, nrow = np, ncol = nt)
  for (tt in seq_len(nt)) {
    std_res_t <- sqrt(taug[, tt]) * residuals_raw[, tt]
    xp_beta <- prop_extract_xt(x_pred, tt) %*% beta
    cond_mean_v <- as.vector(
      cov_state$r_cross_pred_obs %*% (cov_state$prec_obs %*% std_res_t)
    )
    cond_noise_v <- if (isTRUE(draw)) {
      as.vector(t(cov_state$pred_cond_chol) %*% rnorm(np))
    } else {
      rep(0, np)
    }

    pred_draws[, tt] <- as.vector(xp_beta) +
      lambda * z_pred[, tt] +
      sigma_pred[, tt] * (cond_mean_v + cond_noise_v)
  }

  pred_draws
}
