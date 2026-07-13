# Diagnostics for E1 skew-t constant-column A/B experiments
#
# Goal: for the worst skew-t cells (largest |bias| or largest SD),
# inspect chain-level correlations:
#   corr(beta0, lambda), corr(beta0, sigma_xi), corr(beta0, tau_mean)
#
# This script is read-only with respect to the model code; it only runs MCMC.

suppressMessages({
  library(fields)
  library(emulator)
})

prop_dir <- "D:/Github/spatial-skew-t/code/R/prop"

source(file.path(prop_dir, "load_ar2.R"))
source(file.path(prop_dir, "prop_utils.R"))
source(file.path(prop_dir, "prop_basis.R"))
source(file.path(prop_dir, "prop_covariance.R"))
source(file.path(prop_dir, "prop_imputation.R"))
source(file.path(prop_dir, "prop_modules.R"))
source(file.path(prop_dir, "mcmc_prop.R"))

simdata_candidates <- c(
  file.path(prop_dir, "../../analysis/simstudy/simdata.RData"),
  file.path(prop_dir, "../analysis/simstudy/simdata.RData")
)
simdata_path <- simdata_candidates[file.exists(simdata_candidates)][1]
stopifnot(!is.na(simdata_path))
load(simdata_path)

obs <- c(rep(TRUE, nrow(y) - ntest), rep(FALSE, ntest))

prepare_obs_data <- function(setting_id, dataset_id) {
  y.d <- y[, , dataset_id, setting_id]
  list(
    y = y.d[obs, , drop = FALSE],
    x = x[obs, , , drop = FALSE],
    s = s[obs, , drop = FALSE]
  )
}

e1_seed <- function(setting_id, method_id, dataset_id, prop_k, injection_a) {
  as.integer(method_id) * 100000L +
    as.integer(prop_k) * 1000L +
    as.integer(setting_id) * 100L +
    as.integer(dataset_id) +
    as.integer(round(injection_a * 10))
}

run_diag_cell <- function(dataset_id, setting_id, injection_a, trim_constant,
                           prop_k, method_id,
                           iters = 3000L, burn = 1000L) {
  dat <- prepare_obs_data(setting_id, dataset_id)
  y_inj <- dat$y + injection_a

  set.seed(e1_seed(setting_id, method_id, dataset_id, prop_k, injection_a))

  fit <- mcmc(
    y = y_inj,
    x = dat$x,
    s = dat$s,
    method = "t",
    skew = TRUE,
    thresh.all = 0,
    thresh.quant = TRUE,
    nknots = 1L,
    iterplot = FALSE,
    iters = iters,
    burn = burn,
    update = max(500L, iters %/% 10L),
    thin = 1L,
    min.s = c(0, 0),
    max.s = c(10, 10),
    temporalw = FALSE,
    temporaltau = FALSE,
    temporalz = FALSE,
    prop_k = prop_k,
    prop_cov_update_every = 1L,
    prop_trim_constant = isTRUE(trim_constant)
  )

  beta0 <- fit$beta[, 1]
  lambda <- fit$lambda
  sigma_xi <- fit$prop$sigma_xi

  tau_mean <- if (!is.null(fit$tau)) {
    apply(fit$tau, 1, mean) # average across nknots, nt
  } else {
    NA_real_
  }

  z_tmean <- if (!is.null(fit$z)) {
    # fit$z: nsaves x nknots x nt; take mean across nknots and nt
    apply(fit$z, 1, mean)
  } else {
    NA_real_
  }

  total_location <- beta0 + lambda * z_tmean

  corr_b_l <- suppressWarnings(cor(beta0, lambda))
  corr_b_sx <- suppressWarnings(cor(beta0, sigma_xi))
  corr_b_tau <- suppressWarnings(cor(beta0, tau_mean))
  corr_b_z <- suppressWarnings(cor(beta0, z_tmean))
  corr_b_total <- suppressWarnings(cor(beta0, total_location))

  cat("\n=== DIAG ===\n")
  cat(sprintf("dataset=%d setting=%d injection_a=%.1f trim=%s\n",
              dataset_id, setting_id, injection_a,
              ifelse(trim_constant, "TRUE", "FALSE")))
  cat(sprintf("nsave=%d iters=%d burn=%d\n", length(beta0), iters, burn))
  cat(sprintf("beta0: mean=%.3f sd=%.3f min=%.3f max=%.3f\n",
              mean(beta0), sd(beta0), min(beta0), max(beta0)))
  cat(sprintf("lambda: mean=%.3f sd=%.3f min=%.3f max=%.3f\n",
              mean(lambda, na.rm = TRUE),
              sd(lambda, na.rm = TRUE),
              min(lambda, na.rm = TRUE),
              max(lambda, na.rm = TRUE)))
  cat(sprintf("sigma_xi: mean=%.3f sd=%.3f\n", mean(sigma_xi), sd(sigma_xi)))
  cat(sprintf("tau_mean: mean=%.3f sd=%.3f\n", mean(tau_mean), sd(tau_mean)))
  cat(sprintf("z_tmean: mean=%.3f sd=%.3f\n", mean(z_tmean, na.rm = TRUE), sd(z_tmean, na.rm = TRUE)))
  cat(sprintf("total_location=beta0+lambda*z_tmean: mean=%.3f sd=%.3f\n",
              mean(total_location, na.rm = TRUE), sd(total_location, na.rm = TRUE)))
  cat(sprintf("corr(beta0, lambda)=%.3f\n", corr_b_l))
  cat(sprintf("corr(beta0, sigma_xi)=%.3f\n", corr_b_sx))
  cat(sprintf("corr(beta0, tau_mean)=%.3f\n", corr_b_tau))
  cat(sprintf("corr(beta0, z_tmean)=%.3f\n", corr_b_z))
  cat(sprintf("corr(beta0, total_location)=%.3f\n", corr_b_total))

  invisible(NULL)
}

prop_k <- 20L
method_id <- 2L # skew-t K=1 in prop simstudy catalog

# Worst abs_bias cell from the replicate summary:
#   replicate=6, phase B_skewt_k1 (setting=4), trim=TRUE, injection_a=5
run_diag_cell(dataset_id = 6L, setting_id = 4L, injection_a = 5,
              trim_constant = TRUE,
              prop_k = prop_k, method_id = method_id,
              iters = 3000L, burn = 1000L)

# Worst beta0_sd cell from the replicate summary:
#   replicate=3, phase B_skewt_k1 (setting=4), trim=FALSE, injection_a=0
run_diag_cell(dataset_id = 3L, setting_id = 4L, injection_a = 0,
              trim_constant = FALSE,
              prop_k = prop_k, method_id = method_id,
              iters = 3000L, burn = 1000L)

