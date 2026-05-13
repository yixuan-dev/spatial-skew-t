#########################################################################
# scores.R - Stage 1 of the simstudy post-fit pipeline.
#
# Loads fits from <results_dir>/<setting>-<method>-<dataset>[-K<mrts_k>].RData,
# computes Brier / Quantile scores against the held-out test set in the
# loaded simdata, captures parameter quantile intervals + elapsed_sec,
# and writes a single .RData cache to:
#
#   output/results/scores<setting><suffix>.RData
#
# Mirrors code/analysis/simstudy_prop/scores-prop.R, with prop_k -> mrts_k
# and fits_dir -> results_dir. mrts_k = 0 is treated as the baseline
# (filename has no -K suffix); mrts_k > 0 maps to the MRTS variant
# `<setting>-<method>-<dataset>-K<mrts_k>.RData`.
#
# This script supersedes results1.R..results8.R, combineresults.R,
# results.R, and mrts_scores.R.
#
# Usage:
#   Rscript scores.R --setting=<id>
#                    [--data=<path>]
#                    [--methods=<spec>]   default 1:5
#                    [--datasets=<spec>]  default 1..nsets
#                    [--mrts_k=<spec>]    default = auto-detect from results_dir
#
# Examples:
#   Rscript scores.R --setting=4
#   Rscript scores.R --setting=4 --mrts_k="c(0,10,20)" --datasets="1:10"
#   Rscript scores.R --setting=1 --data=simdata_def.RData
#
# results_dir resolution:
#   1. SIMSTUDY_RESULTS_DIR env var (or legacy SIMSTUDY_MRTS_RESULTS_DIR)
#   2. else derive_results_dir(data_path, "results")
#########################################################################

rm(list = ls())

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0L) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]),
                               winslash = "/", mustWork = FALSE)
  if (dir.exists(dirname(script_path))) setwd(dirname(script_path))
}

source("./mrts_cov_helpers.R")
source("../../R/auxfunctions.R")

# ---- CLI parsing -----------------------------------------------------
cli_args <- commandArgs(trailingOnly = TRUE)
parsed <- extract_leading_flags(
  cli_args,
  c("data", "setting", "methods", "datasets", "mrts_k")
)
flags <- parsed$values

if (is.null(flags$setting) || !nzchar(flags$setting)) {
  stop("scores.R: --setting=<id> is required.", call. = FALSE)
}

# ---- load data + resolve results_dir / suffix ------------------------
data_path   <- resolve_simstudy_data_path(flags$data)
load(data_path)
data_suffix <- derive_data_suffix(data_path)

setting_ids <- parse_index_expr(flags$setting, "setting")
if (length(setting_ids) != 1L) {
  stop("scores.R: --setting must be a single integer.", call. = FALSE)
}
setting <- as.integer(setting_ids[[1]])
max_setting <- as.integer(dim(y)[4])
if (setting < 1L || setting > max_setting) {
  stop(sprintf("setting must be in 1..%d", max_setting), call. = FALSE)
}

results_dir_env <- trimws(Sys.getenv("SIMSTUDY_RESULTS_DIR", unset = ""))
if (!nzchar(results_dir_env)) {
  results_dir_env <- trimws(Sys.getenv("SIMSTUDY_MRTS_RESULTS_DIR", unset = ""))
}
results_dir <- if (nzchar(results_dir_env)) {
  results_dir_env
} else {
  derive_results_dir(data_path, "results")
}
if (!dir.exists(results_dir)) {
  stop(sprintf("results directory not found: %s", results_dir), call. = FALSE)
}

# ---- methods / datasets / mrts_k -------------------------------------
methods <- if (!is.null(flags$methods) && nzchar(flags$methods)) {
  m <- parse_index_expr(flags$methods, "methods")
  if (any(m < 1L | m > 8L)) stop("methods must be in 1..8", call. = FALSE)
  sort(unique(as.integer(m)))
} else {
  1:5
}

datasets <- if (!is.null(flags$datasets) && nzchar(flags$datasets)) {
  ds <- parse_index_expr(flags$datasets, "datasets")
  max_d <- as.integer(dim(y)[3])
  if (any(ds < 1L | ds > max_d)) {
    stop(sprintf("datasets must be in 1..%d", max_d), call. = FALSE)
  }
  as.integer(ds)
} else {
  seq_len(as.integer(dim(y)[3]))
}

# mrts_k auto-detect: scan results_dir for both baseline and -K<K> patterns
auto_detect_mrts_ks <- function(results_dir, setting) {
  files <- list.files(results_dir)
  base_pat <- sprintf("^%d-\\d+-\\d+\\.RData$", setting)
  k_pat    <- sprintf("^%d-\\d+-\\d+-K(\\d+)\\.RData$", setting)
  ks <- integer(0)
  if (any(grepl(base_pat, files))) ks <- c(ks, 0L)
  m <- regmatches(files, regexec(k_pat, files))
  k_vals <- vapply(m, function(x) if (length(x) >= 2L) as.integer(x[2]) else NA_integer_,
                   integer(1))
  k_vals <- k_vals[!is.na(k_vals)]
  ks <- c(ks, k_vals)
  sort(unique(ks))
}

mrts_ks <- if (!is.null(flags$mrts_k) && nzchar(flags$mrts_k)) {
  k <- parse_index_expr(flags$mrts_k, "mrts_k")
  sort(unique(as.integer(k)))
} else {
  detected <- auto_detect_mrts_ks(results_dir, setting)
  if (length(detected) == 0L) {
    stop(sprintf(
      "No fits matching '%d-*-*[-K<K>].RData' found in %s; pass --mrts_k=<spec> explicitly.",
      setting, results_dir
    ), call. = FALSE)
  }
  detected
}

# ---- output target ---------------------------------------------------
out_dir  <- "output/results"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, sprintf("scores%d%s.RData", setting, data_suffix))

cat(sprintf(
  "scores: setting=%d data=%s results_dir=%s out=%s\n  methods=%s datasets=%s mrts_k=%s\n",
  setting, data_path, results_dir, out_file,
  paste(methods, collapse = ","),
  paste(range(datasets), collapse = ".."),
  paste(mrts_ks, collapse = ",")
))

# ---- preallocate arrays ---------------------------------------------
nsets <- length(datasets)
nmeth <- length(methods)
nks   <- length(mrts_ks)

probs     <- c(0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995)
intervals <- c(0.01, 0.025, 0.05, 0.1, 0.9, 0.95, 0.975, 0.99)

dn_score <- list(
  quantile = as.character(probs),
  dataset  = as.character(datasets),
  method   = as.character(methods),
  mrts_k   = as.character(mrts_ks)
)
dn_param <- list(
  interval = as.character(intervals),
  dataset  = as.character(datasets),
  method   = as.character(methods),
  mrts_k   = as.character(mrts_ks)
)

quant.score <- array(NA_real_, dim = c(length(probs), nsets, nmeth, nks), dimnames = dn_score)
brier.score <- array(NA_real_, dim = c(length(probs), nsets, nmeth, nks), dimnames = dn_score)

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
                                     mrts_k  = as.character(mrts_ks)))

skew.methods <- c(2L, 4L, 7L, 8L)

obs <- c(rep(TRUE, nrow(y) - ntest), rep(FALSE, ntest))

# ---- score loop ------------------------------------------------------
save_checkpoint <- function() {
  # Save the input fits directory under the label `fits_dir` (not
  # `results_dir`) so loading the cache in tables.R does not shadow
  # tables.R's local `results_dir = "output/results"` variable.
  fits_dir <- results_dir
  save(quant.score, brier.score,
       beta.0, beta.1, beta.2,
       tau.alpha, tau.beta, rho, nu, gamma, lambda,
       elapsed_sec,
       probs, intervals, mrts_ks, datasets, methods, setting,
       data_path, data_suffix, fits_dir,
       file = out_file)
}

for (di in seq_along(datasets)) {
  set <- datasets[di]
  thresholds <- quantile(y[, , set, setting], probs = probs, na.rm = TRUE)
  validate   <- y[!obs, , set, setting]

  for (mi in seq_along(methods)) {
    method <- methods[mi]
    for (ki in seq_along(mrts_ks)) {
      mrts_k <- mrts_ks[ki]
      f <- build_simstudy_result_file(
        results_dir = results_dir,
        setting_id  = setting,
        method_id   = method,
        dataset_id  = set,
        mrts_k      = if (mrts_k == 0L) NA_integer_ else mrts_k
      )
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

  if (set %% 10L == 0L) {
    save_checkpoint()
    cat("  -> checkpoint saved (", out_file, ")\n", sep = "")
  }
}

save_checkpoint()
cat("\nWrote ", out_file, "\n", sep = "")
