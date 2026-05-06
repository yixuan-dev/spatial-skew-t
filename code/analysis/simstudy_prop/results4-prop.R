#########################################################################
# results4-prop.R
#
# Analogue of ../simstudy/results4.R for the simstudy_prop fits.
#
# data setting (fixed): 4 - skew t-1 (lambda = 3)
#
# analysis methods (prop catalog):
#   1 - Gaussian
#   2 - skew t, K = 1
#   3 - t,      K = 1, threshold q(0.80)
#   4 - skew t, K = 5
#   5 - t,      K = 5, threshold q(0.80)
#
# prop basis ranks: 10, 20, 30, 40, 50, 60, 70
#
# Score & parameter-interval arrays gain a 4th index for prop_k:
#   quant.score / brier.score : [probs, dataset, method, prop_k]
#   beta.0 ... lambda          : [intervals, dataset, method, prop_k]
#
# Outputs:
#   output/results/scores4-prop.RData
#########################################################################

rm(list = ls())

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]),
                               winslash = "/", mustWork = FALSE)
  if (dir.exists(dirname(script_path))) setwd(dirname(script_path))
}

source("../../R/prop/auxfunctions.R")

resolve_simdata <- function() {
  candidates <- c("./simdata.RData", "../simstudy/simdata.RData")
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0L) stop("simdata.RData not found.", call. = FALSE)
  normalizePath(hit[1], winslash = "/", mustWork = TRUE)
}
load(resolve_simdata())

setting    <- 4L
datasets   <- 1:50
methods    <- 1:5
prop_ks    <- c(10, 20, 30, 40, 50, 60, 70)
fits_dir   <- "fits"
out_dir    <- "output/results"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file   <- file.path(out_dir, "scores4-prop.RData")

nsets   <- length(datasets)
nmeth   <- length(methods)
nks     <- length(prop_ks)

obs <- c(rep(TRUE, nrow(y) - ntest), rep(FALSE, ntest))

probs     <- c(0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995)
intervals <- c(0.01, 0.025, 0.05, 0.1, 0.9, 0.95, 0.975, 0.99)

dn_score <- list(
  quantile = as.character(probs),
  dataset  = as.character(datasets),
  method   = as.character(methods),
  prop_k   = as.character(prop_ks)
)
dn_param <- list(
  interval = as.character(intervals),
  dataset  = as.character(datasets),
  method   = as.character(methods),
  prop_k   = as.character(prop_ks)
)

quant.score <- array(NA_real_, dim = c(length(probs), nsets, nmeth, nks),
                     dimnames = dn_score)
brier.score <- array(NA_real_, dim = c(length(probs), nsets, nmeth, nks),
                     dimnames = dn_score)

beta.0    <- array(NA_real_, dim = c(length(intervals), nsets, nmeth, nks), dimnames = dn_param)
beta.1    <- array(NA_real_, dim = c(length(intervals), nsets, nmeth, nks), dimnames = dn_param)
beta.2    <- array(NA_real_, dim = c(length(intervals), nsets, nmeth, nks), dimnames = dn_param)
tau.alpha <- array(NA_real_, dim = c(length(intervals), nsets, nmeth, nks), dimnames = dn_param)
tau.beta  <- array(NA_real_, dim = c(length(intervals), nsets, nmeth, nks), dimnames = dn_param)
rho       <- array(NA_real_, dim = c(length(intervals), nsets, nmeth, nks), dimnames = dn_param)
nu        <- array(NA_real_, dim = c(length(intervals), nsets, nmeth, nks), dimnames = dn_param)
gamma     <- array(NA_real_, dim = c(length(intervals), nsets, nmeth, nks), dimnames = dn_param)
lambda    <- array(NA_real_, dim = c(length(intervals), nsets, nmeth, nks), dimnames = dn_param)

elapsed_sec <- array(NA_real_, dim = c(nsets, nmeth, nks),
                     dimnames = list(dataset = as.character(datasets),
                                     method  = as.character(methods),
                                     prop_k  = as.character(prop_ks)))

skew.methods <- c(2L, 4L)

for (di in seq_along(datasets)) {
  set <- datasets[di]
  thresholds <- quantile(y[, , set, setting], probs = probs, na.rm = TRUE)
  validate   <- y[!obs, , set, setting]

  for (mi in seq_along(methods)) {
    method <- methods[mi]
    for (ki in seq_along(prop_ks)) {
      prop_k <- prop_ks[ki]
      f <- file.path(fits_dir,
                     sprintf("%d-%d-%d-p%d.RData", setting, method, set, prop_k))
      if (!file.exists(f)) {
        cat("missing: ", f, "\n", sep = "")
        next
      }
      env <- new.env(parent = emptyenv())
      load(f, envir = env)
      fit <- env$fit.1
      if (is.null(fit)) { cat("no fit.1: ", f, "\n", sep = ""); next }

      if (!is.null(fit$yp)) {
        brier.score[, di, mi, ki] <- BrierScore(fit$yp, thresholds, validate)
        quant.score[, di, mi, ki] <- QuantScore(fit$yp, probs, validate)
      }

      if (!is.null(fit$beta) && ncol(fit$beta) >= 3L) {
        beta.0[, di, mi, ki] <- quantile(fit$beta[, 1], probs = intervals, na.rm = TRUE)
        beta.1[, di, mi, ki] <- quantile(fit$beta[, 2], probs = intervals, na.rm = TRUE)
        beta.2[, di, mi, ki] <- quantile(fit$beta[, 3], probs = intervals, na.rm = TRUE)
      }
      if (!is.null(fit$tau.alpha))
        tau.alpha[, di, mi, ki] <- quantile(fit$tau.alpha, probs = intervals, na.rm = TRUE)
      if (!is.null(fit$tau.beta))
        tau.beta[, di, mi, ki]  <- quantile(fit$tau.beta,  probs = intervals, na.rm = TRUE)
      if (!is.null(fit$rho))
        rho[, di, mi, ki]       <- quantile(fit$rho,       probs = intervals, na.rm = TRUE)
      if (!is.null(fit$nu))
        nu[, di, mi, ki]        <- quantile(fit$nu,        probs = intervals, na.rm = TRUE)
      if (!is.null(fit$gamma))
        gamma[, di, mi, ki]     <- quantile(fit$gamma,     probs = intervals, na.rm = TRUE)
      if (method %in% skew.methods && !is.null(fit$lambda))
        lambda[, di, mi, ki]    <- quantile(fit$lambda,    probs = intervals, na.rm = TRUE)

      if (exists("runtime_info", envir = env, inherits = FALSE)) {
        rt <- env$runtime_info
        if (!is.null(rt$elapsed_sec))
          elapsed_sec[di, mi, ki] <- as.numeric(rt$elapsed_sec)
      }

      rm(fit, env)
    }
    cat(sprintf("dataset %d  method %d done\n", set, method))
  }

  if (set %% 10 == 0) {
    save(quant.score, brier.score,
         beta.0, beta.1, beta.2,
         tau.alpha, tau.beta, rho, nu, gamma, lambda,
         elapsed_sec,
         probs, intervals, prop_ks, datasets, methods, setting,
         file = out_file)
    cat("  -> checkpoint saved (", out_file, ")\n", sep = "")
  }
}

save(quant.score, brier.score,
     beta.0, beta.1, beta.2,
     tau.alpha, tau.beta, rho, nu, gamma, lambda,
     elapsed_sec,
     probs, intervals, prop_ks, datasets, methods, setting,
     file = out_file)

cat("\nWrote ", out_file, "\n", sep = "")
