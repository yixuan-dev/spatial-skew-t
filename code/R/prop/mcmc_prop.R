mcmc <- function(y, s, x, s.pred = NULL, x.pred = NULL,
                 min.s, max.s,
                 thresh.all = 0, thresh.quant = TRUE,
                 thresh.site.specific = FALSE, thresh.site = NULL,
                 nknots = 1, fixknots = FALSE,
                 keep.knots = FALSE,
                 skew = TRUE, method = "t",
                 cov.model = "matern",
                 temporaltau = FALSE, temporalw = FALSE, temporalz = FALSE,
                 ar2_tau = FALSE, ar2_w = FALSE, ar2_z = FALSE,
                 beta.init = NULL, tau.init = 1,
                 tau.alpha.init = 0.1, tau.beta.init = 0.1,
                 rho.init = 5, nu.init = 0.5, gamma.init = 0.5,
                 z.init = 0, lambda.init = NULL, knots.init = NULL,
                 phi.tau.init = NULL, phi.w.init = NULL, phi.z.init = NULL,
                 beta.m = 0, beta.s = 20,
                 tau.alpha.min = 0.1, tau.alpha.max = 10, tau.alpha.by = 0.1,
                 tau.beta.a = 1, tau.beta.b = 1,
                 rho.upper = NULL,
                 lognu.m = -1.2, lognu.s = 1, nu.upper = NULL, fix.nu = FALSE,
                 lambda.m = 0, lambda.s = 20,
                 iters = 5000, burn = 1000, update = 100, thin = 1,
                 iterplot = FALSE,
                 beta_update_method = c("joint", "split"),
                 prop_k = 10L,
                 prop_sigma_floor = 1e-6,
                 prop_cov_update_every = 1L,
                 prop_basis_tol = 1e-8,
                 prop_sigma_eps_sq = 0) {
  beta_update_method <- match.arg(beta_update_method)

  if (!method %in% c("gaussian", "t")) {
    stop("prop backend currently supports method = 'gaussian' or 't' only.", call. = FALSE)
  }
  if (isTRUE(thresh.site.specific)) {
    stop("prop backend does not yet support site-specific thresholds.", call. = FALSE)
  }
  if (isTRUE(temporaltau) || isTRUE(temporalw) || isTRUE(temporalz) ||
      isTRUE(ar2_tau) || isTRUE(ar2_w) || isTRUE(ar2_z)) {
    stop("prop backend does not yet support temporal dependence blocks.", call. = FALSE)
  }
  if (method == "gaussian" && (skew || nknots != 1L)) {
    stop("Gaussian prop phase supports only skew = FALSE and nknots = 1.", call. = FALSE)
  }

  y <- as.matrix(y)
  s <- as.matrix(s)
  x <- array(x, dim = dim(x))
  ns <- nrow(y)
  nt <- ncol(y)
  p <- dim(x)[3]

  predictions <- !is.null(s.pred) && !is.null(x.pred)
  np <- if (predictions) nrow(s.pred) else 0L
  if (predictions) {
    s.pred <- as.matrix(s.pred)
    x.pred <- array(x.pred, dim = dim(x.pred))
  }

  if (thresh.quant) {
    if (thresh.all < 0 || thresh.all > 1) {
      stop("When thresh.quant = TRUE, thresh.all must be between 0 and 1.", call. = FALSE)
    }
    thresh.all.q <- quantile(y, thresh.all, na.rm = TRUE)
  } else {
    thresh.all.q <- thresh.all
  }

  thresh.mtx <- matrix(thresh.all.q, ns, nt)
  thresh.obs <- !is.na(y) & (y < thresh.mtx)
  thresholded <- any(thresh.obs)

  missing.obs <- is.na(y)
  missing <- any(missing.obs)
  if (missing) {
    y[missing.obs] <- mean(y, na.rm = TRUE)
  }

  basis_obj <- build_prop_basis(
    s_obs = s,
    s_pred = if (predictions) s.pred else NULL,
    k = as.integer(prop_k),
    tol = prop_basis_tol
  )
  q_obs <- basis_obj$q_obs
  q_pred <- basis_obj$q_pred

  if (!fixknots) {
    if (is.null(knots.init)) {
      lower.1 <- min.s[1] + (max.s[1] - min.s[1]) / 6
      upper.1 <- max.s[1] - (max.s[1] - min.s[1]) / 6
      lower.2 <- min.s[2] + (max.s[2] - min.s[2]) / 6
      upper.2 <- max.s[2] - (max.s[2] - min.s[2]) / 6
      knots <- array(NA_real_, dim = c(nknots, 2, nt))
      knots[, 1, ] <- runif(nknots * nt, lower.1, upper.1)
      knots[, 2, ] <- runif(nknots * nt, lower.2, upper.2)
    } else {
      knots <- array(knots.init, dim = c(nknots, 2, nt))
    }
  } else {
    if (is.null(knots.init)) {
      stop("You must specify knots.init when fixknots = TRUE.", call. = FALSE)
    }
    knots <- array(knots.init, dim = c(nknots, 2, nt))
  }

  if (nknots == 1L) {
    g <- matrix(1L, ns, nt)
  } else {
    g <- matrix(0L, ns, nt)
    for (tt in seq_len(nt)) {
      g[, tt] <- mem(s, knots[, , tt])
    }
  }

  if (is.null(beta.init)) {
    x_stack <- prop_build_design_stack(x)
    beta <- tryCatch(
      qr.solve(x_stack, as.vector(y)),
      error = function(e) c(mean(y), rep(0, p - 1L))
    )
  } else {
    beta <- prop_expand_scalar(beta.init, p, "beta.init")
  }

  x.beta <- matrix(0, ns, nt)
  for (tt in seq_len(nt)) {
    x.beta[, tt] <- prop_extract_xt(x, tt) %*% beta
  }

  if (method == "gaussian") {
    tau <- matrix(as.numeric(tau.init)[1], nknots, nt)
    taug <- matrix(as.numeric(tau.init)[1], ns, nt)
  } else {
    tau <- matrix(tau.init, nknots, nt)
    taug <- matrix(0, ns, nt)
    for (tt in seq_len(nt)) {
      taug[, tt] <- tau[g[, tt], tt]
    }
  }

  if (skew) {
    z <- matrix(z.init, nknots, nt)
    zg <- matrix(0, ns, nt)
    for (tt in seq_len(nt)) {
      zg[, tt] <- z[g[, tt], tt]
    }
    lambda <- if (is.null(lambda.init)) e1071::skewness(y, na.rm = TRUE) else lambda.init
  } else {
    z <- matrix(0, nknots, nt)
    zg <- matrix(0, ns, nt)
    lambda <- 0
  }

  tau.alpha <- tau.alpha.init
  tau.beta <- tau.beta.init
  prop_cov_update_every <- max(1L, as.integer(prop_cov_update_every))

  save_iters <- seq.int(from = burn + 1L, to = iters, by = thin)
  nsave <- length(save_iters)
  beta_keep <- matrix(NA_real_, nrow = nsave, ncol = p)
  tau_keep <- array(NA_real_, dim = c(nsave, nknots, nt))
  tau.alpha_keep <- rep(NA_real_, nsave)
  tau.beta_keep <- rep(NA_real_, nsave)
  rho_keep <- rep(NA_real_, nsave)
  nu_keep <- rep(NA_real_, nsave)
  gamma_keep <- rep(NA_real_, nsave)
  lambda_keep <- if (skew) rep(NA_real_, nsave) else NULL
  z_keep <- if (skew) array(NA_real_, dim = c(nsave, nknots, nt)) else NULL
  yp_keep <- if (predictions) array(NA_real_, dim = c(nsave, np, nt)) else NULL
  knots_keep <- if (keep.knots && nknots > 1L) array(NA_real_, dim = c(nsave, nknots, 2, nt)) else NULL

  prop_sigma_keep <- rep(NA_real_, nsave)
  prop_rank_keep <- rep(NA_integer_, nsave)
  prop_loglik_keep <- rep(NA_real_, nsave)

  acc.w <- att.w <- mh.w <- matrix(0.1, nknots, nt)
  acc.tau <- att.tau <- matrix(1, nknots, nt)
  mh.tau <- matrix(0.05, nknots, nt)
  nparts.tau <- matrix(1, nknots, nt)

  prop_refresh_xbeta <- function(beta_vec) {
    for (tt in seq_len(nt)) {
      x.beta[, tt] <<- prop_extract_xt(x, tt) %*% beta_vec
    }
  }

  prop_current_residuals <- function() {
    y - x.beta - lambda * zg
  }

  prop_standardized_residuals <- function(residuals) {
    sqrt(taug) * residuals
  }

  prop_update_covariance <- function(std_residuals) {
    update_obj <- prop_closed_form_update(
      residuals = std_residuals,
      q_obs = q_obs,
      sigma_eps_sq = prop_sigma_eps_sq,
      sigma_floor = prop_sigma_floor
    )
    state <- prop_make_cov_state(q_obs = q_obs, q_pred = q_pred, update_obj = update_obj)
    prop_prepare_prediction_state(state)
  }

  prop_compute_loglik <- function(cov_state, std_residuals) {
    quad_sum <- 0
    for (tt in seq_len(ncol(std_residuals))) {
      quad_sum <- quad_sum + quad.form(cov_state$prec_obs, std_residuals[, tt])
    }
    nt * cov_state$logdet_prec_obs - 0.5 * quad_sum
  }

  prop_predict_draws <- function(cov_state, residuals, beta_vec) {
    if (!predictions) {
      return(NULL)
    }

    pred_draws <- matrix(NA_real_, nrow = np, ncol = nt)
    for (tt in seq_len(nt)) {
      gp <- if (nknots == 1L) rep(1L, np) else mem(s.pred, knots[, , tt])
      zgp <- if (skew) z[gp, tt] else rep(0, np)
      siggp <- 1 / sqrt(tau[gp, tt])
      std_res_t <- sqrt(taug[, tt]) * residuals[, tt]
      xp.beta <- prop_extract_xt(x.pred, tt) %*% beta_vec
      cond_mean_v <- as.vector(cov_state$r_cross_pred_obs %*% (cov_state$prec_obs %*% std_res_t))
      cond_noise_v <- as.vector(t(cov_state$pred_cond_chol) %*% rnorm(np))
      pred_draws[, tt] <- as.vector(xp.beta) + lambda * zgp + siggp * (cond_mean_v + cond_noise_v)
    }

    pred_draws
  }

  prop_refresh_xbeta(beta)
  residuals <- prop_current_residuals()
  std_residuals <- prop_standardized_residuals(residuals)
  cov_state <- prop_update_covariance(std_residuals)
  save_idx <- 1L

  for (iter in seq_len(iters)) {
    prec <- cov_state$prec_obs

    if (thresholded || missing) {
      mu <- x.beta + lambda * zg
      y <- prop_impute_latent_y(
        y = y,
        mu = mu,
        taug = taug,
        cov_state = cov_state,
        censor_obs = if (thresholded) thresh.obs else NULL,
        thresh.mtx = if (thresholded) thresh.mtx else NULL,
        missing_obs = if (missing) missing.obs else NULL
      )
    }

    if (beta_update_method == "joint") {
      this.update <- updateBeta(
        beta.m = beta.m, beta.s = beta.s,
        x = x, y = y, zg = zg, taug = taug,
        prec = prec, skew = skew,
        lambda = lambda, lambda.m = lambda.m,
        lambda.s = lambda.s
      )
      beta <- as.vector(this.update$beta)
      lambda <- this.update$lambda
    } else {
      this.update <- updateBetaonly(
        beta.m = beta.m, beta.s = beta.s,
        x = x, y = y, zg = zg, taug = taug,
        prec = prec, lambda = lambda
      )
      beta <- as.vector(this.update$beta)
      if (skew) {
        this.update <- updateZonly(
          lambda.m = lambda.m, lambda.s = lambda.s,
          x = x, y = y, beta = beta, zg = zg,
          taug = taug, prec = prec
        )
        lambda <- this.update$lambda
      } else {
        lambda <- 0
      }
    }

    prop_refresh_xbeta(beta)

    if ((nknots > 1L) && (!fixknots)) {
      this.update <- updateKnotsTS(
        phi = 0, knots = knots,
        g = g, ts = FALSE,
        tau = tau, z = z, s = s,
        min.s = min.s, max.s = max.s,
        x.beta = x.beta,
        lambda = lambda, y = y,
        prec = prec, att = att.w,
        acc = acc.w, mh = mh.w,
        att.phi = 0, acc.phi = 0, mh.phi = 0
      )
      knots <- this.update$knots
      g <- this.update$g
      taug <- this.update$taug
      zg <- this.update$zg
      acc.w <- this.update$acc
      att.w <- this.update$att

      if (iter < max(1, burn %/% 2)) {
        this.update <- mhupdate(
          acc = acc.w, att = att.w, mh = mh.w,
          nattempts = nknots * nt
        )
        acc.w <- this.update$acc
        att.w <- this.update$att
        mh.w <- this.update$mh
      }
    }

    residuals <- prop_current_residuals()
    if (method == "gaussian") {
      tau_scalar <- updateTauGaus(
        res = residuals,
        prec = prec, tau.alpha = tau.alpha,
        tau.beta = tau.beta
      )
      tau <- matrix(tau_scalar, nknots, nt)
      taug <- matrix(tau_scalar, ns, nt)
    } else {
      this.update <- updateTau(
        tau = tau, taug = taug, g = g,
        res = residuals, nparts.tau = nparts.tau,
        prec = prec, z = z,
        tau.alpha = tau.alpha,
        tau.beta = tau.beta, skew = skew,
        att = att.tau, acc = acc.tau,
        mh = mh.tau
      )
      tau <- this.update$tau
      taug <- this.update$taug
      acc.tau <- this.update$acc
      att.tau <- this.update$att

      tau.alpha <- updateTauAlpha(
        tau = tau, tau.beta = tau.beta,
        tau.alpha.min = tau.alpha.min,
        tau.alpha.max = tau.alpha.max,
        tau.alpha.by = tau.alpha.by
      )
      tau.beta <- updateTauBeta(
        tau = tau, tau.alpha = tau.alpha,
        tau.beta.a = tau.beta.a,
        tau.beta.b = tau.beta.b
      )
    }

    if (skew) {
      residuals <- prop_current_residuals()
      this.update <- updateZ(
        y = y, x.beta = x.beta, zg = zg,
        prec = prec, tau = tau, mu = x.beta + lambda * zg,
        taug = taug, g = g, lambda = lambda
      )
      z <- this.update$z
      zg <- this.update$zg
    }

    residuals <- prop_current_residuals()
    std_residuals <- prop_standardized_residuals(residuals)
    if (iter %% prop_cov_update_every == 0L) {
      cov_state <- prop_update_covariance(std_residuals)
    }

    if (iter %in% save_iters) {
      beta_keep[save_idx, ] <- beta
      tau_keep[save_idx, , ] <- tau
      tau.alpha_keep[save_idx] <- tau.alpha
      tau.beta_keep[save_idx] <- tau.beta
      rho_keep[save_idx] <- cov_state$sigma_xi_sq
      nu_keep[save_idx] <- cov_state$effective_rank
      gamma_keep[save_idx] <- cov_state$logdet_prec_obs
      prop_sigma_keep[save_idx] <- cov_state$sigma_xi_sq
      prop_rank_keep[save_idx] <- cov_state$effective_rank
      prop_loglik_keep[save_idx] <- prop_compute_loglik(cov_state, std_residuals)

      if (skew) {
        lambda_keep[save_idx] <- lambda
        z_keep[save_idx, , ] <- z
      }
      if (!is.null(knots_keep)) {
        knots_keep[save_idx, , , ] <- knots
      }
      if (predictions) {
        yp_keep[save_idx, , ] <- prop_predict_draws(cov_state, residuals, beta)
      }
      save_idx <- save_idx + 1L
    }

    if (iter %% update == 0L) {
      cat(sprintf(
        "[prop phase3] iter=%d/%d sigma_xi=%.4f rank=%d lambda=%.4f\n",
        iter, iters, cov_state$sigma_xi_sq, cov_state$effective_rank, lambda
      ))
    }
  }

  list(
    tau = tau_keep,
    beta = beta_keep,
    tau.alpha = tau.alpha_keep,
    tau.beta = tau.beta_keep,
    rho = rho_keep,
    nu = nu_keep,
    gamma = gamma_keep,
    y = NULL,
    yp = yp_keep,
    lambda = lambda_keep,
    z = z_keep,
    knots = knots_keep,
    avgparts = NULL,
    phi.z = NULL,
    phi.w = NULL,
    phi.tau = NULL,
    prop = list(
      phase = 3L,
      requested_k = basis_obj$requested_k,
      basis_original_cols = basis_obj$original_cols,
      basis_kept_cols = basis_obj$kept_cols,
      basis_rank_kept = basis_obj$rank_kept,
      basis_source = basis_obj$source,
      cov_update_every = prop_cov_update_every,
      thresholded = thresholded,
      censored_cells = sum(thresh.obs),
      missing = missing,
      sigma_xi = prop_sigma_keep,
      effective_rank = prop_rank_keep,
      loglik = prop_loglik_keep,
      logdet_prec = gamma_keep
    )
  )
}
