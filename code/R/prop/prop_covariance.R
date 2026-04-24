prop_select_lstar <- function(d_vals, trace_s, n, sigma_eps_sq) {
  if (length(d_vals) == 0L) {
    return(0L)
  }

  lstar <- 0L
  cum_vals <- cumsum(d_vals)
  for (ll in seq_along(d_vals)) {
    denom <- n - ll
    if (denom <= 0) {
      next
    }
    rhs <- max((trace_s - cum_vals[ll]) / denom, sigma_eps_sq)
    if (d_vals[ll] > rhs) {
      lstar <- ll
    }
  }
  as.integer(lstar)
}

prop_closed_form_update <- function(residuals, q_obs, sigma_eps_sq = 0, sigma_floor = 1e-6) {
  residuals <- as.matrix(residuals)
  n <- nrow(residuals)
  tt <- ncol(residuals)

  projected <- crossprod(q_obs, residuals)
  h_mat <- tcrossprod(projected) / tt
  h_mat <- prop_symmetrize(h_mat)
  eig_h <- eigen(h_mat, symmetric = TRUE)
  d_vals <- pmax(eig_h$values, 0)

  trace_s <- mean(colSums(residuals ^ 2))
  lstar <- prop_select_lstar(d_vals = d_vals, trace_s = trace_s, n = n, sigma_eps_sq = sigma_eps_sq)
  retained_sum <- if (lstar > 0L) sum(d_vals[seq_len(lstar)]) else 0
  sigma_xi_sq <- max((trace_s - retained_sum) / (n - lstar) - sigma_eps_sq, 0)
  sigma_xi_sq <- max(sigma_xi_sq, sigma_floor)
  lowrank_eigs <- pmax(d_vals - sigma_xi_sq - sigma_eps_sq, 0)

  list(
    rotation = eig_h$vectors,
    sigma_xi_sq = sigma_xi_sq,
    lowrank_eigs = lowrank_eigs,
    projected_eigs = d_vals,
    effective_rank = sum(lowrank_eigs > 1e-10),
    trace_s = trace_s,
    lstar = lstar
  )
}

prop_make_cov_state <- function(q_obs, q_pred, update_obj) {
  q_obs <- as.matrix(q_obs)
  rotation <- update_obj$rotation
  if (is.null(rotation)) {
    rotation <- diag(ncol(q_obs))
  }
  q_obs <- q_obs %*% rotation
  q_pred <- if (is.null(q_pred)) NULL else as.matrix(q_pred) %*% rotation
  d_vals <- as.numeric(update_obj$lowrank_eigs)
  sigma_xi_sq <- as.numeric(update_obj$sigma_xi_sq)

  g_obs <- tcrossprod(sweep(q_obs, 2, d_vals, "*"), q_obs) +
    diag(sigma_xi_sq, nrow(q_obs))
  obs_diag <- sigma_xi_sq + rowSums(sweep(q_obs ^ 2, 2, d_vals, "*"))
  obs_diag <- pmax(obs_diag, 1e-10)
  obs_sqrt <- sqrt(obs_diag)
  r_obs <- sweep(sweep(g_obs, 1, obs_sqrt, "/"), 2, obs_sqrt, "/")
  r_obs <- prop_symmetrize(r_obs)
  chol_r_obs <- prop_stable_chol(r_obs)
  prec_obs <- chol2inv(chol_r_obs)
  logdet_r_obs <- 2 * sum(log(diag(chol_r_obs)))
  logdet_prec_obs <- -sum(log(diag(chol_r_obs)))

  state <- list(
    q_obs = q_obs,
    q_pred = q_pred,
    basis_rotation = rotation,
    lowrank_eigs = d_vals,
    sigma_xi_sq = sigma_xi_sq,
    g_obs = g_obs,
    r_obs = r_obs,
    prec_obs = prec_obs,
    a_obs_diag = obs_diag,
    a_obs_sqrt = obs_sqrt,
    logdet_r_obs = logdet_r_obs,
    logdet_prec_obs = logdet_prec_obs,
    effective_rank = update_obj$effective_rank,
    trace_s = update_obj$trace_s,
    lstar = update_obj$lstar
  )

  positive_idx <- which(d_vals > 1e-10)
  noise_var_obs <- pmax(sigma_xi_sq / obs_diag, 1e-10)
  state$noise_var_obs <- noise_var_obs
  state$noise_sd_obs <- sqrt(noise_var_obs)
  state$inv_noise_var_obs <- 1 / noise_var_obs
  state$factor_rank <- length(positive_idx)
  state$factor_indices <- positive_idx

  if (length(positive_idx) == 0L) {
    state$factor_loadings_obs <- matrix(0, nrow(q_obs), 0)
    state$factor_post_prec <- matrix(0, 0, 0)
    state$factor_post_cov <- matrix(0, 0, 0)
    state$factor_post_chol <- matrix(0, 0, 0)
  } else {
    factor_loadings_obs <- sweep(q_obs[, positive_idx, drop = FALSE], 1, obs_sqrt, "/")
    factor_loadings_obs <- sweep(
      factor_loadings_obs,
      2,
      sqrt(d_vals[positive_idx]),
      "*"
    )
    factor_post_prec <- diag(1, length(positive_idx)) +
      crossprod(
        factor_loadings_obs,
        sweep(factor_loadings_obs, 1, noise_var_obs, "/")
      )
    factor_post_cov <- chol2inv(prop_stable_chol(factor_post_prec))

    state$factor_loadings_obs <- factor_loadings_obs
    state$factor_post_prec <- factor_post_prec
    state$factor_post_cov <- prop_symmetrize(factor_post_cov)
    state$factor_post_chol <- prop_stable_chol(state$factor_post_cov)
  }

  if (!is.null(state$q_pred)) {
    pred_diag <- sigma_xi_sq + rowSums(sweep(state$q_pred ^ 2, 2, d_vals, "*"))
    pred_diag <- pmax(pred_diag, 1e-10)
    pred_sqrt <- sqrt(pred_diag)

    g_cross <- sweep(state$q_pred, 2, d_vals, "*") %*% t(state$q_obs)
    g_pred <- tcrossprod(sweep(state$q_pred, 2, d_vals, "*"), state$q_pred) +
      diag(sigma_xi_sq, nrow(state$q_pred))

    r_cross <- sweep(sweep(g_cross, 1, pred_sqrt, "/"), 2, state$a_obs_sqrt, "/")
    r_pred <- sweep(sweep(g_pred, 1, pred_sqrt, "/"), 2, pred_sqrt, "/")

    state$a_pred_diag <- pred_diag
    state$a_pred_sqrt <- pred_sqrt
    state$r_cross_pred_obs <- r_cross
    state$r_pred <- prop_symmetrize(r_pred)
  }

  state
}

prop_apply_ginv_obs <- function(state, x) {
  x_mat <- prop_as_matrix(x)
  sigma_xi_sq <- state$sigma_xi_sq
  out <- x_mat / sigma_xi_sq

  if (length(state$lowrank_eigs) > 0L) {
    coeff <- state$lowrank_eigs / (sigma_xi_sq * (sigma_xi_sq + state$lowrank_eigs))
    out <- out - state$q_obs %*% sweep(crossprod(state$q_obs, x_mat), 1, coeff, "*")
  }

  if (is.null(dim(x))) {
    return(as.vector(out))
  }
  out
}

prop_apply_rinv_obs <- function(state, x) {
  x_mat <- prop_as_matrix(x)
  scaled <- sweep(x_mat, 1, state$a_obs_sqrt, "*")
  out <- prop_apply_ginv_obs(state, scaled)
  out <- sweep(out, 1, state$a_obs_sqrt, "*")

  if (is.null(dim(x))) {
    return(as.vector(out))
  }
  out
}

prop_prepare_prediction_state <- function(state) {
  if (is.null(state$q_pred)) {
    return(state)
  }

  cross_rinv <- prop_apply_rinv_obs(state, t(state$r_cross_pred_obs))
  cond_cov <- prop_symmetrize(state$r_pred - state$r_cross_pred_obs %*% cross_rinv)
  state$pred_cond_cov <- cond_cov
  state$pred_cond_chol <- prop_stable_chol(cond_cov)
  state
}
