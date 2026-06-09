#########################################################################
# time_block_helpers.R - shared CLI / filename / seed / catalog helpers
# for the time-block forecast simulation study.
#
# Mirrors code/analysis/simstudy/helpers.R, with the MRTS-basis
# axis replaced by the time-block axis. The generic CLI / path helpers are
# copied verbatim so the two studies parse arguments identically.
#########################################################################

# ---- generic CLI / path helpers (shared with the Morris study) --------

parse_index_expr <- function(expr_str, arg_name) {
  expr <- trimws(expr_str)
  if (is.na(expr) || !nzchar(expr)) {
    stop(sprintf("%s specification cannot be empty", arg_name), call. = FALSE)
  }
  if (grepl("^\\(.*\\)$", expr)) {
    expr <- paste0("c", expr)
  }
  values_raw <- tryCatch(
    eval(parse(text = expr), envir = baseenv()),
    error = function(e) {
      stop(sprintf(
        "Invalid %s expression '%s'. Use forms like 1:5, c(1,2,4), (1,2,4)",
        arg_name, expr_str
      ), call. = FALSE)
    }
  )
  if (!is.numeric(values_raw) || length(values_raw) == 0) {
    stop(sprintf("%s must evaluate to a non-empty numeric vector", arg_name), call. = FALSE)
  }
  values_int <- as.integer(values_raw)
  if (any(is.na(values_int)) || any(values_raw != values_int)) {
    stop(sprintf("%s must contain integers only", arg_name), call. = FALSE)
  }
  sort(unique(values_int))
}

# Generic leading --key=value (or --key value) parser.
extract_leading_flags <- function(args, flag_names) {
  values <- setNames(vector("list", length(flag_names)), flag_names)
  while (length(args) > 0L && startsWith(args[1], "--")) {
    first <- args[1]
    matched <- FALSE
    for (name in flag_names) {
      flag_eq <- paste0("--", name, "=")
      flag_word <- paste0("--", name)
      if (startsWith(first, flag_eq)) {
        if (!is.null(values[[name]])) {
          stop(sprintf("Specify --%s only once", name), call. = FALSE)
        }
        values[[name]] <- sub(paste0("^", flag_eq), "", first)
        args <- if (length(args) > 1L) args[-1] else character(0)
        matched <- TRUE
        break
      } else if (first == flag_word) {
        if (!is.null(values[[name]])) {
          stop(sprintf("Specify --%s only once", name), call. = FALSE)
        }
        if (length(args) < 2L) {
          stop(sprintf("Missing value after --%s", name), call. = FALSE)
        }
        values[[name]] <- args[2]
        args <- if (length(args) > 2L) args[-c(1, 2)] else character(0)
        matched <- TRUE
        break
      }
    }
    if (!matched) {
      stop(sprintf(
        "Unknown leading flag: %s. Place known flags (%s) before positional args.",
        first, paste0("--", flag_names, collapse = ", ")
      ), call. = FALSE)
    }
  }
  list(args = args, values = values)
}

# Derive results / output suffix from the dataset filename, exactly as the
# Morris study does (simdata.RData -> "", simdata_def.RData -> "_def").
derive_results_dir <- function(data_path, default_dir = "results") {
  base <- tools::file_path_sans_ext(basename(data_path))
  if (identical(base, "simdata")) {
    return(default_dir)
  }
  suffix <- sub("^simdata", "", base)
  if (!nzchar(suffix)) paste0(default_dir, "_", base) else paste0(default_dir, suffix)
}

derive_data_suffix <- function(data_path) {
  base <- tools::file_path_sans_ext(basename(data_path))
  if (identical(base, "simdata")) return("")
  s <- sub("^simdata", "", base)
  if (!nzchar(s)) paste0("_", base) else s
}

resolve_simstudy_data_path <- function(data_arg = NULL) {
  if (!is.null(data_arg) && nzchar(data_arg)) {
    if (file.exists(data_arg)) {
      return(normalizePath(data_arg, winslash = "/", mustWork = TRUE))
    }
    stop(sprintf("Data file not found: %s", data_arg), call. = FALSE)
  }
  if (!file.exists("./simdata.RData")) {
    stop("Unable to locate simdata.RData; run setup.R first.", call. = FALSE)
  }
  normalizePath("./simdata.RData", winslash = "/", mustWork = TRUE)
}

format_runtime_timestamp <- function(timestamp) {
  if (length(timestamp) == 0 || anyNA(timestamp)) return(NA_character_)
  format(as.POSIXct(timestamp, tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

build_tbf_runtime_info <- function(started_at, elapsed_sec, outputfile,
                                   setting_id, dataset_id, runner,
                                   method_id = NA_integer_,
                                   method_key = NA_character_,
                                   control = NULL) {
  finished_at <- started_at + as.numeric(elapsed_sec)
  list(
    schema_version = 1L,
    runner = as.character(runner),
    output_file = outputfile,
    setting_id = as.integer(setting_id),
    dataset_id = as.integer(dataset_id),
    method_id = if (is.na(method_id)) NA_integer_ else as.integer(method_id),
    method_key = if (is.na(method_key)) NA_character_ else as.character(method_key),
    started_at_utc = format_runtime_timestamp(started_at),
    finished_at_utc = format_runtime_timestamp(finished_at),
    elapsed_sec = unname(as.numeric(elapsed_sec)),
    control = control
  )
}

# ---- filenames + seeds ------------------------------------------------
# fit/forecast cache: <results_dir>/<setting>-<method>-<dataset>.RData
build_tbf_result_file <- function(results_dir, setting_id, method_id, dataset_id) {
  file.path(results_dir, sprintf(
    "%d-%d-%d.RData",
    as.integer(setting_id), as.integer(method_id), as.integer(dataset_id)
  ))
}

# Deterministic seed for one (setting, method, dataset) fit.
get_tbf_seed <- function(setting_id, method_id, dataset_id) {
  as.integer(method_id) * 1000L + as.integer(setting_id) * 100L + as.integer(dataset_id)
}

# ---- method catalog (analysis methods, the --methods axis) ------------
# Method 1 is the i.i.d.-in-time baseline; methods 2-3 activate AR(2)
# temporal pooling. Method 3 freezes the Voronoi membership at the seam
# (the fixed-membership ablation of Remark 8) to separate the tau/z
# signal from knot-flip noise; with K = 1 it coincides with method 2.
get_tbf_method_catalog <- function() {
  data.frame(
    method_id = 1:3,
    method_key = as.character(1:3),
    label = c(
      "Skew-t, K=1, i.i.d. in time",
      "Skew-t, K=1, AR(2) temporal (tau,z,w)",
      "Skew-t, K=1, AR(2) temporal, fixed membership"
    ),
    skew = c(TRUE, TRUE, TRUE),
    nknots = c(1L, 1L, 1L),
    temporal = c(FALSE, TRUE, TRUE),
    ar2 = c(FALSE, TRUE, TRUE),
    fixed_membership = c(FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
}

# ---- setting catalog (data-generating AR(2) strength, --setting axis) -
get_tbf_setting_catalog <- function() {
  data.frame(
    setting_id = 1:4,
    label = c("iid", "weak", "moderate", "strong"),
    phi1 = c(0.00, 0.12, 0.60, 0.80),
    phi2 = c(0.00, -0.05, -0.30, -0.35),
    stringsAsFactors = FALSE
  )
}

# ---- time-block geometry ---------------------------------------------
# AR(2) companion-matrix spectral radius (Definition 5 of the note):
# rho(F) where F = [[phi1, phi2], [1, 0]].
ar2_spectral_radius <- function(phi1, phi2) {
  Fm <- matrix(c(phi1, 1, phi2, 0), 2, 2)
  max(Mod(eigen(Fm, only.values = TRUE)$values))
}

# Effective memory horizon h*(epsilon) = ceil(log eps / log rho)
# (eq. (14) of the note). Returns Inf for the i.i.d. control.
ar2_memory_horizon <- function(phi1, phi2, epsilon = 0.05) {
  if (phi1 == 0 && phi2 == 0) return(0L)
  sr <- ar2_spectral_radius(phi1, phi2)
  if (sr <= 0 || sr >= 1) return(NA_integer_)
  as.integer(ceiling(log(epsilon) / log(sr)))
}

# Expand the stored seam vector into per-block window descriptors.
# Each block b fits on the prefix 1:seam and forecasts leads 1..H, i.e.
# the calendar times (seam, seam + H].
tbf_blocks <- function(block_seams, block_H, nt) {
  lapply(seq_along(block_seams), function(b) {
    seam <- block_seams[b]
    leads <- seq_len(block_H)
    list(
      block_id = b,
      seam = as.integer(seam),
      leads = leads,
      test_times = as.integer(seam + leads),
      train_times = seq_len(seam)
    )
  })
}

#########################################################################
# forecast_block() - posterior predictive forecast over one held-out
# time block (Algorithm 1 of the note).
#
# Implements the three requirements of Section 4.2:
#   (i)   condition on the posterior of the seam state, extracted from
#         each draw's imputed latent trajectory -- never the stationary
#         distribution (Corollary 1 warns this erases the signal);
#   (ii)  propagate parameter + state uncertainty by iterating over all
#         MCMC draws and averaging, with no plug-in of posterior means;
#   (iii) recurse on the latent Gaussian scale, applying the copula
#         transforms (eqs. 8-10) only after the AR(2) recursion.
#
# Arguments:
#   fit       - object returned by mcmc() (AR(2) backend). Must carry the
#               imputed latent trajectories tau [iters x K x T_o] and
#               z [iters x K x T_o], the AR(2) coefficient draws
#               phi.tau / phi.z [iters x 2], and the observation-model
#               draws beta, lambda, tau.alpha, tau.beta, rho, nu, gamma.
#   seam      - last observed time T_o (fit covers times 1:T_o).
#   H         - forecast horizon (number of lead times).
#   x_block   - covariates at the forecast window [ns x H x p].
#   s         - site coordinates [ns x 2].
#   ar2       - TRUE for AR(2) methods, FALSE for the i.i.d. baseline
#               (which draws each latent slice from the N(0,1) marginal,
#               reproducing Corollary 1's stationary predictive).
#
# Returns: predictive sample array yhat [iters x ns x H].
#########################################################################
forecast_block <- function(fit, seam, H, x_block, s, ar2 = TRUE) {
  # The AR(2) backend stores tau/z as keepers[iters, nknots, T_o] but
  # drops the singleton knot axis on return when nknots == 1, so fit$tau
  # arrives as a 2-D [iters, T_o] matrix. Restore the knot axis so the
  # [iter, knot, time] indexing below is uniform for any K.
  tau_traj <- fit$tau
  z_traj <- fit$z
  if (length(dim(tau_traj)) == 2L) {
    d2 <- dim(tau_traj)
    tau_traj <- array(tau_traj, dim = c(d2[1], 1L, d2[2]))
    if (!is.null(z_traj)) z_traj <- array(z_traj, dim = c(d2[1], 1L, d2[2]))
  }
  iters <- dim(tau_traj)[1]
  K <- dim(tau_traj)[2]
  ns <- nrow(s)
  p <- dim(x_block)[3]
  yhat <- array(NA_real_, dim = c(iters, ns, H))

  seam_cols <- c(seam - 1L, seam)  # the two seam states, T_o-1 and T_o

  for (m in seq_len(iters)) {
    a.m <- fit$tau.alpha[m] / 2   # rpotspatTS_arp uses the /2 reparam
    b.m <- fit$tau.beta[m] / 2
    beta.m <- fit$beta[m, ]
    lambda.m <- if (is.null(fit$lambda)) 0 else fit$lambda[m]

    # ---- forecast the latent Gaussian processes tau*, z* -------------
    if (ar2) {
      phi.tau <- fit$phi.tau[m, ]
      phi.z <- fit$phi.z[m, ]
      # seam state on the Gaussian scale: invert the copula at T_o-1, T_o.
      # drop = FALSE then reshape to [K, 2] so the K = 1 case does not
      # collapse to a vector (ar2_recurse indexes seam_state[, j]).
      tau.obs <- matrix(tau_traj[m, , seam_cols, drop = FALSE], K, 2L)
      z.obs <- matrix(z_traj[m, , seam_cols, drop = FALSE], K, 2L)
      tau.seam <- gamma.cop(tau.obs, a.m, b.m)
      sd.z <- 1 / sqrt(tau.obs)
      z.seam <- hn.cop(z.obs, sig = sd.z)
      tau.star <- ar2_recurse(tau.seam, phi.tau, H, K)
      z.star <- ar2_recurse(z.seam, phi.z, H, K)
    } else {
      # i.i.d. baseline: each forecast slice is a fresh N(0,1) draw,
      # i.e. the stationary marginal (Corollary 1).
      tau.star <- matrix(rnorm(K * H), K, H)
      z.star <- matrix(rnorm(K * H), K, H)
    }

    # ---- copula transforms back to observable scale (eqs. 8-9) -------
    tau.fc <- matrix(gamma.invcop(as.vector(tau.star), a.m, b.m), K, H)
    sd.z.fc <- 1 / sqrt(tau.fc)
    z.fc <- matrix(hn.invcop(as.vector(z.star), sig = as.vector(sd.z.fc)), K, H)

    # ---- spatial field at each lead time -----------------------------
    # K = 1: membership is constant, knot 1 serves every site.
    d <- as.matrix(dist(s))
    C <- CorFx(d = d, gamma = fit$gamma[m], rho = fit$rho[m], nu = fit$nu[m])
    t.chol.C <- t(chol(C))
    for (h in seq_len(H)) {
      g <- rep(1L, ns)            # constant Voronoi membership for K = 1
      x.beta <- if (p == 1) {
        matrix(x_block[, h, ], ns, 1) * beta.m
      } else {
        x_block[, h, ] %*% beta.m
      }
      mu <- x.beta + lambda.m * z.fc[g, h]
      sdg <- 1 / sqrt(tau.fc[g, h])
      yhat[m, , h] <- mu + t.chol.C %*% matrix(rnorm(ns, 0, sdg), ns, 1)
    }
  }
  yhat
}

# AR(2) latent recursion: given the 2-column seam state (columns T_o-1,
# T_o) advance H steps with innovation variance fixed by Yule-Walker so
# the marginal variance stays 1 (eq. 5 of the note).
ar2_recurse <- function(seam_state, phi, H, K) {
  phi1 <- phi[1]; phi2 <- phi[2]
  innov.sd <- if (check_ar2_stability(phi1, phi2)) {
    get_ar2_standard_moments(phi1, phi2)$innov_sd
  } else {
    # non-stationary posterior tail draw: clamp to the i.i.d. marginal
    # rather than let the recursion explode (Remark 7d).
    1
  }
  out <- matrix(NA_real_, K, H)
  prev2 <- seam_state[, 1]
  prev1 <- seam_state[, 2]
  for (h in seq_len(H)) {
    cur <- phi1 * prev1 + phi2 * prev2 + rnorm(K, 0, innov.sd)
    out[, h] <- cur
    prev2 <- prev1
    prev1 <- cur
  }
  out
}
