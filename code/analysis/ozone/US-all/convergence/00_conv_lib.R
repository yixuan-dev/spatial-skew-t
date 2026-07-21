# convergence/00_conv_lib.R
# ---------------------------------------------------------------------------
# Shared library for the US-all MCMC convergence diagnostics. Safe to source
# from any script: defines functions and metadata only, no side effects.
#
# All paths returned by conv_paths() are relative to the US-all directory,
# so every runner sets its working directory there first (see the preamble in
# 01_extract_chains.R).
# ---------------------------------------------------------------------------

# The nine settings named in the thesis table tab:ozone-top2-brier, plus the
# metadata needed to cross-check what we find inside each fit object.
# `backend` is the expectation; extraction re-detects it from the phi shape
# and warns on mismatch.
conv_settings <- data.frame(
  setting  = c(8L, 51L, 55L, 58L, 68L, 111L, 120L, 124L, 204L),
  method   = "t",
  skew     = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE),
  nknots   = c(5L, 1L, 5L, 6L, 9L, 7L, 10L, 15L, 1L),
  thresh   = c(50, 0, 50, 50, 75, 50, 50, 75, 0),
  temporal = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  backend  = c("legacy", "legacy", "legacy", "legacy", "legacy",
               "ar2", "ar2", "ar2", "legacy"),
  mrts_k   = c(NA, NA, NA, NA, NA, NA, NA, NA, 10L),
  stringsAsFactors = FALSE
)

conv_paths <- function(create = TRUE) {
  p <- list(
    results_raw = "results",
    cache  = file.path("output", "us-all", "results", "convergence"),
    tables = file.path("output", "us-all", "tables"),
    plots  = file.path("output", "us-all", "plots", "convergence"),
    logs   = file.path("output", "us-all", "logs", "convergence")
  )
  if (create) {
    for (d in p[c("cache", "tables", "plots", "logs")]) {
      if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
    }
  }
  p
}

# --------------------------------------------------------------------------
# Scalar-chain extraction (shape-driven parameter registry)
#
# Takes one fit element (fit[[fold]] from a saved run, or the direct return
# of mcmc()/mcmc_ar2()) and returns a [ndraws x P] matrix of the scalar
# chains worth diagnosing. Per-knot tau/z chains are deliberately excluded:
# the knots are exchangeable, so those chains label-switch; only the
# permutation-invariant field summaries mean.log.tau / mean.z are kept.
# --------------------------------------------------------------------------

detect_backend <- function(f) {
  if (is.matrix(f$phi.tau) || is.matrix(f$phi.z) || is.matrix(f$phi.w)) return("ar2")
  if (!is.null(f$phi.tau) || !is.null(f$phi.z) || !is.null(f$phi.w)) return("legacy")
  "legacy-nots"
}

extract_scalar_chains <- function(f) {
  b <- as.matrix(f$beta)
  ndraws <- nrow(b)
  labs <- c("beta[int]", "beta[cmaq]")
  if (ncol(b) > 2) labs <- c(labs, sprintf("beta[mrts%d]", seq_len(ncol(b) - 2)))
  colnames(b) <- labs[seq_len(ncol(b))]
  cols <- list(b)

  add_vec <- function(v, name) {
    if (is.null(v)) return(invisible(NULL))
    m <- matrix(as.numeric(v), ncol = 1, dimnames = list(NULL, name))
    stopifnot(nrow(m) == ndraws)
    cols[[length(cols) + 1L]] <<- m
  }
  add_phi <- function(v, name) {
    if (is.null(v)) return(invisible(NULL))
    if (is.matrix(v) && ncol(v) == 2) {
      colnames(v) <- sprintf("%s[%d]", name, 1:2)
      stopifnot(nrow(v) == ndraws)
      cols[[length(cols) + 1L]] <<- v
    } else {
      add_vec(as.numeric(v), name)
    }
  }
  # rowMeans over all non-draw dims; column-major reshape keeps draws as rows
  field_summary <- function(a, fun = identity) {
    rowMeans(matrix(fun(as.numeric(a)), nrow = ndraws))
  }

  add_vec(f$tau.alpha, "tau.alpha")
  add_vec(f$tau.beta, "tau.beta")
  add_vec(f$rho, "rho")
  add_vec(f$nu, "nu")
  add_vec(f$gamma, "gamma")
  add_vec(f$lambda, "lambda")
  add_phi(f$phi.z, "phi.z")
  add_phi(f$phi.w, "phi.w")
  add_phi(f$phi.tau, "phi.tau")
  if (!is.null(f$tau)) add_vec(field_summary(f$tau, log), "mean.log.tau")
  if (!is.null(f$z)) add_vec(field_summary(f$z), "mean.z")

  do.call(cbind, cols)
}

# Cross-check the extracted chains against conv_settings expectations.
check_expected_params <- function(chains, setting_row, backend_detected) {
  expected <- c("beta[int]", "beta[cmaq]", "tau.alpha", "tau.beta",
                "rho", "nu", "gamma", "mean.log.tau")
  if (setting_row$skew) expected <- c(expected, "lambda", "mean.z")
  if (setting_row$temporal) {
    phis <- c("phi.z", "phi.w", "phi.tau")
    if (setting_row$backend == "ar2") {
      phis <- as.vector(outer(phis, 1:2, function(p, k) sprintf("%s[%d]", p, k)))
    }
    expected <- c(expected, phis)
  }
  if (!is.na(setting_row$mrts_k)) {
    n_mrts <- sum(grepl("^beta\\[mrts", colnames(chains)))
    if (n_mrts < 1) {
      warning(sprintf("setting %s: mrts_k=%d but no beta[mrts*] chains found",
                      setting_row$setting, setting_row$mrts_k))
    }
    expected <- c(expected, sprintf("beta[mrts%d]", seq_len(n_mrts)))
  }
  missing <- setdiff(expected, colnames(chains))
  extra <- setdiff(colnames(chains), expected)
  if (length(missing)) {
    warning(sprintf("setting %s: expected but missing chains: %s",
                    setting_row$setting, paste(missing, collapse = ", ")))
  }
  if (length(extra)) {
    warning(sprintf("setting %s: unexpected chains: %s",
                    setting_row$setting, paste(extra, collapse = ", ")))
  }
  if (backend_detected != setting_row$backend &&
      !(backend_detected == "legacy-nots" && !setting_row$temporal)) {
    warning(sprintf("setting %s: expected backend %s but detected %s",
                    setting_row$setting, setting_row$backend, backend_detected))
  }
  invisible(NULL)
}

# --------------------------------------------------------------------------
# Rank-normalized split R-hat (Vehtari, Gelman, Simpson, Carpenter, Buerkner
# 2021), plus coda-based ESS and Geweke wrappers.
# --------------------------------------------------------------------------

# Reshape one chain into nsplit contiguous pseudo-chains (drops leading draws
# if the length is not divisible).
split_chain <- function(v, nsplit = 4L) {
  n <- length(v)
  keep <- n - n %% nsplit
  matrix(v[(n - keep + 1L):n], ncol = nsplit)
}

# Split each column (chain) of a matrix in half -> 2M columns.
split_chains_half <- function(m) {
  n <- nrow(m) - nrow(m) %% 2L
  m <- m[(nrow(m) - n + 1L):nrow(m), , drop = FALSE]
  h <- n %/% 2L
  do.call(cbind, lapply(seq_len(ncol(m)), function(j) {
    cbind(m[1:h, j], m[(h + 1L):n, j])
  }))
}

rank_normalize <- function(m) {
  S <- length(m)
  z <- qnorm((rank(as.numeric(m), ties.method = "average") - 3 / 8) / (S + 1 / 4))
  matrix(z, nrow = nrow(m))
}

rhat_basic <- function(m) {
  N <- nrow(m)
  W <- mean(apply(m, 2, var))
  B <- N * var(colMeans(m))
  if (!is.finite(W) || W <= 0) return(NA_real_)
  sqrt(((N - 1) / N * W + B / N) / W)
}

# Input: matrix with columns = (already split) chains, raw draws.
# Returns c(rhat_bulk, rhat_tail, rhat) with NAs for degenerate chains.
rhat_rank_norm <- function(m) {
  if (sd(as.numeric(m)) == 0 || any(apply(m, 2, sd) == 0)) {
    return(c(rhat_bulk = NA_real_, rhat_tail = NA_real_, rhat = NA_real_))
  }
  bulk <- rhat_basic(rank_normalize(m))
  folded <- abs(m - median(m))
  tail <- rhat_basic(rank_normalize(folded))
  c(rhat_bulk = bulk, rhat_tail = tail, rhat = max(bulk, tail))
}

ess_coda <- function(v) {
  tryCatch(unname(coda::effectiveSize(coda::mcmc(v))), error = function(e) NA_real_)
}

ess_coda_multi <- function(chain_list) {
  tryCatch(
    unname(coda::effectiveSize(coda::mcmc.list(lapply(chain_list, coda::mcmc)))),
    error = function(e) NA_real_
  )
}

geweke_z <- function(v) {
  z <- tryCatch(unname(coda::geweke.diag(coda::mcmc(v))$z), error = function(e) NA_real_)
  if (!is.finite(z)) NA_real_ else z
}

# Diagnose every column of a chains matrix as a single chain (pass 1).
diagnose_single_chain <- function(chains, nsplit = 4L) {
  out <- lapply(colnames(chains), function(pn) {
    v <- chains[, pn]
    degenerate <- sd(v) == 0
    r <- if (degenerate) c(rhat_bulk = NA_real_, rhat_tail = NA_real_, rhat = NA_real_)
         else rhat_rank_norm(split_chain(v, nsplit))
    data.frame(
      param = pn, mean = mean(v), sd = sd(v),
      rhat_bulk = r[["rhat_bulk"]], rhat_tail = r[["rhat_tail"]], rhat = r[["rhat"]],
      ess = if (degenerate) NA_real_ else ess_coda(v),
      geweke_z = if (degenerate) NA_real_ else geweke_z(v),
      degenerate = degenerate,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

# pass/marginal/fail per Vehtari et al. thresholds (see README).
classify_flag <- function(rhat, ess, degenerate) {
  ifelse(degenerate, "degenerate",
    ifelse(!is.finite(rhat) | !is.finite(ess), "NA",
      ifelse(rhat > 1.05 | ess < 100, "fail",
        ifelse(rhat > 1.01 | ess < 400, "marginal", "pass"))))
}

# --------------------------------------------------------------------------
# Cache handling
# --------------------------------------------------------------------------

cache_file <- function(setting, paths = conv_paths()) {
  file.path(paths$cache, sprintf("chains-us-all-%d.rds", setting))
}

extract_setting_from_rdata <- function(setting, paths = conv_paths()) {
  src <- file.path(paths$results_raw, sprintf("us-all-%d.RData", setting))
  if (!file.exists(src)) stop(sprintf("missing fit file: %s", src), call. = FALSE)
  row <- conv_settings[conv_settings$setting == setting, , drop = FALSE]
  stopifnot(nrow(row) == 1L)

  e <- new.env()
  cat(sprintf("loading %s (%.0f MB) ...\n", src, file.size(src) / 2^20))
  t0 <- proc.time()[3]
  load(src, envir = e)
  fit <- get("fit", envir = e)
  runtime_info <- if (exists("runtime_info", envir = e, inherits = FALSE)) {
    get("runtime_info", envir = e)
  } else NULL
  stopifnot(is.list(fit), length(fit) >= 2L)

  folds <- list()
  backends <- character(2)
  for (d in 1:2) {
    backends[d] <- detect_backend(fit[[d]])
    folds[[d]] <- extract_scalar_chains(fit[[d]])
    check_expected_params(folds[[d]], row, backends[d])
  }
  rm(fit, e); invisible(gc(verbose = FALSE))
  cat(sprintf("  extracted setting %d in %.1f s: fold1 %d x %d, fold2 %d x %d\n",
              setting, proc.time()[3] - t0,
              nrow(folds[[1]]), ncol(folds[[1]]), nrow(folds[[2]]), ncol(folds[[2]])))

  list(
    setting = setting,
    backend_expected = row$backend,
    backend_detected = backends,
    folds = folds,
    runtime_info = runtime_info,
    source_file = src,
    source_mtime = file.mtime(src),
    extracted_at = Sys.time()
  )
}

load_setting_slim <- function(setting, paths = conv_paths()) {
  cf <- cache_file(setting, paths)
  if (file.exists(cf)) return(readRDS(cf))
  slim <- extract_setting_from_rdata(setting, paths)
  saveRDS(slim, cf)
  slim
}

# --------------------------------------------------------------------------
# Multichain rerun specs (04/05). Chain 1 is the production fold-1 chain
# (seed setting*100 + 1, default inits); chains 2-4 are new overdispersed
# runs. All values verified inside the samplers' prior/proposal support:
# gamma in (1e-4, 0.9999), rho in (0, 5), nu in (0, 10], AR2 phi pairs in
# the stationarity triangle. Additionally (rho.init, nu.init) must keep the
# 400-site Matern correlation numerically invertible: imputeY does a plain
# chol(cor) with no eigen fallback, and rho=4.5 with nu=2.5 makes cor
# singular (all correlations ~1) and kills the run at initialization.
# --------------------------------------------------------------------------

conv_multichain_settings <- c(55L, 111L)

conv_chain_seed <- function(setting, chain) setting * 1000L + 100L + chain

conv_chain_inits <- list(
  "2" = list(tau.init = 1, tau.alpha.init = 5, tau.beta.init = 2,
             rho.init = 0.3, nu.init = 0.15, gamma.init = 0.15,
             lambda.init = -3, z.init = 1, phi = c(0.8, -0.5)),
  "3" = list(tau.init = 0.005, tau.alpha.init = 9, tau.beta.init = 5,
             rho.init = 3, nu.init = 0.8, gamma.init = 0.85,
             lambda.init = 6, z.init = 3, phi = c(-0.4, 0.2)),
  "4" = list(tau.init = 0.2, tau.alpha.init = 1, tau.beta.init = 0.5,
             rho.init = 2, nu.init = 1, gamma.init = 0.35,
             lambda.init = 0, z.init = 0.1, phi = c(0.3, 0.3))
)

multichain_file <- function(setting, chain, paths = conv_paths()) {
  file.path(paths$cache, sprintf("multichain-%d-fold1-chain%d.rds", setting, chain))
}
