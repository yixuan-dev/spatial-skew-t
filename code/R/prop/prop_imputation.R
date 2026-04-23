prop_draw_truncnorm <- function(mn, sd, lower = -Inf, upper = Inf, prob_tol = 1e-10) {
  mn <- as.numeric(mn)
  sd <- pmax(as.numeric(sd), prob_tol)
  n <- length(mn)

  if (length(lower) == 1L) {
    lower <- rep(lower, n)
  }
  if (length(upper) == 1L) {
    upper <- rep(upper, n)
  }

  lower_u <- pnorm(lower, mean = mn, sd = sd)
  upper_u <- pnorm(upper, mean = mn, sd = sd)
  lower_u[is.infinite(lower) & lower < 0] <- 0
  upper_u[is.infinite(upper) & upper > 0] <- 1

  width_u <- pmax(upper_u - lower_u, 0)
  u <- lower_u + runif(n) * width_u
  u <- pmin(pmax(u, prob_tol), 1 - prob_tol)
  draws <- qnorm(u, mean = mn, sd = sd)

  degenerate <- !is.finite(draws) | (width_u <= prob_tol)
  if (any(degenerate)) {
    fallback <- mn
    has_upper <- is.finite(upper)
    has_lower <- is.finite(lower)

    fallback[has_upper] <- pmin(
      fallback[has_upper],
      upper[has_upper] - prob_tol * pmax(1, abs(upper[has_upper]))
    )
    fallback[has_lower] <- pmax(
      fallback[has_lower],
      lower[has_lower] + prob_tol * pmax(1, abs(lower[has_lower]))
    )
    fallback <- pmin(fallback, upper)
    fallback <- pmax(fallback, lower)
    draws[degenerate] <- fallback[degenerate]
  }

  draws
}

prop_sample_factor_scores <- function(std_residual, cov_state) {
  if (cov_state$factor_rank < 1L) {
    return(numeric(0))
  }

  weighted_residual <- std_residual / cov_state$noise_var_obs
  post_mean <- cov_state$factor_post_cov %*%
    crossprod(cov_state$factor_loadings_obs, weighted_residual)

  as.vector(post_mean + t(cov_state$factor_post_chol) %*% rnorm(cov_state$factor_rank))
}

prop_impute_latent_y <- function(y, mu, taug, cov_state,
                                 censor_obs = NULL, thresh.mtx = NULL,
                                 missing_obs = NULL) {
  y <- as.matrix(y)
  mu <- as.matrix(mu)
  taug <- as.matrix(taug)
  ns <- nrow(y)
  nt <- ncol(y)

  if (is.null(censor_obs)) {
    censor_obs <- matrix(FALSE, ns, nt)
  } else {
    censor_obs <- as.matrix(censor_obs)
  }

  if (is.null(missing_obs)) {
    missing_obs <- matrix(FALSE, ns, nt)
  } else {
    missing_obs <- as.matrix(missing_obs)
  }

  if (!any(censor_obs) && !any(missing_obs)) {
    return(y)
  }

  factor_loadings <- cov_state$factor_loadings_obs
  noise_sd <- cov_state$noise_sd_obs
  y_impute <- y

  for (tt in seq_len(nt)) {
    flagged <- censor_obs[, tt] | missing_obs[, tt]
    if (!any(flagged)) {
      next
    }

    sqrt_taug_t <- sqrt(taug[, tt])
    std_residual_t <- sqrt_taug_t * (y_impute[, tt] - mu[, tt])
    eta_t <- prop_sample_factor_scores(std_residual_t, cov_state)

    flagged_idx <- which(flagged)
    cond_mean <- rep(0, length(flagged_idx))
    if (cov_state$factor_rank > 0L) {
      cond_mean <- as.vector(factor_loadings[flagged_idx, , drop = FALSE] %*% eta_t)
    }

    upper_std <- rep(Inf, length(flagged_idx))
    censored_idx <- which(censor_obs[, tt])
    if (length(censored_idx) > 0L) {
      censored_pos <- match(censored_idx, flagged_idx)
      upper_std[censored_pos] <- sqrt_taug_t[censored_idx] *
        (thresh.mtx[censored_idx, tt] - mu[censored_idx, tt])
    }

    std_draw <- prop_draw_truncnorm(
      mn = cond_mean,
      sd = noise_sd[flagged_idx],
      upper = upper_std
    )
    y_impute[flagged_idx, tt] <- mu[flagged_idx, tt] + std_draw / sqrt_taug_t[flagged_idx]
  }

  y_impute
}
