#########################################################################
# scores-prop.R - Stage 1 of the simstudy_prop post-fit pipeline.
#
# Loads prop fits from <results_dir>/<setting>-<method>-<dataset>-p<K>.RData,
#   (results_dir = results/, results_def/, ... depending on the dataset)
# computes Brier / Quantile scores against the held-out test set in the
# loaded simdata, captures parameter quantile intervals + elapsed_sec,
# and writes a single .RData cache to:
#
#   output/results/scores<setting>-prop<suffix>.RData
#
# This script supersedes the per-setting results4-prop.R / results6-prop.R
# pair and the score-computation halves of analyze_mrts.R / results-prop.R.
#
# Usage:
#   Rscript scores-prop.R --setting=<id>
#                         [--data=<path>]
#                         [--methods=<spec>]   default 1:5
#                         [--datasets=<spec>]  default 1..nsets
#                         [--prop_k=<spec>]    default = auto-detect from results
#
# Examples:
#   Rscript scores-prop.R --setting=4
#   Rscript scores-prop.R --setting=1 --data=simdata_def.RData
#   Rscript scores-prop.R --setting=4 --prop_k="c(20,30)" --datasets="1:10"
#
# Filename suffix:
#   simdata.RData      -> scores<setting>-prop.RData      (from results/)
#   simdata_def.RData  -> scores<setting>-prop_def.RData  (from results_def/)
#
# results_dir resolution:
#   1. SIMSTUDY_PROP_RESULTS_DIR env var if set
#   2. else derive_prop_results_dir(data_path, "results")
#########################################################################

rm(list = ls())

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0L) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]),
                               winslash = "/", mustWork = FALSE)
  if (dir.exists(dirname(script_path))) setwd(dirname(script_path))
}

source("./prop_simstudy_helpers.R")
source("../../R/ar2/auxfunctions.R")

# ---- CLI parsing -----------------------------------------------------
cli_args <- commandArgs(trailingOnly = TRUE)
parsed <- extract_leading_flags(
  cli_args,
  c("data", "setting", "methods", "datasets", "prop_k")
)
flags <- parsed$values

if (is.null(flags$setting) || !nzchar(flags$setting)) {
  stop("scores-prop.R: --setting=<id> is required.", call. = FALSE)
}

# ---- load data + resolve results_dir / suffix ---------------------------
data_path   <- resolve_simstudy_data_path(flags$data)
load(data_path)
data_suffix <- derive_data_suffix(data_path)

setting <- parse_setting_spec(flags$setting, y)

results_dir_env <- trimws(Sys.getenv("SIMSTUDY_PROP_RESULTS_DIR", unset = ""))
results_dir <- if (nzchar(results_dir_env)) {
  results_dir_env
} else {
  derive_prop_results_dir(data_path, "fits")
}
if (!dir.exists(results_dir)) {
  stop(sprintf("results directory not found: %s", results_dir), call. = FALSE)
}

# ---- methods / datasets / prop_k -------------------------------------
methods <- if (!is.null(flags$methods) && nzchar(flags$methods)) {
  parse_prop_methods_spec(flags$methods)
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

# prop_k auto-detect: scan results_dir for "<setting>-*-*-p<K>.RData"
auto_detect_prop_ks <- function(results_dir, setting) {
  pat <- sprintf("^%d-\\d+-\\d+-[pP](\\d+)\\.RData$", setting)
  files <- list.files(results_dir, pattern = pat)
  if (length(files) == 0L) return(integer(0))
  m <- regmatches(files, regexec(pat, files))
  ks <- vapply(m, function(x) as.integer(x[2]), integer(1))
  sort(unique(ks))
}

prop_ks <- if (!is.null(flags$prop_k) && nzchar(flags$prop_k)) {
  parse_prop_k_spec(flags$prop_k)
} else {
  detected <- auto_detect_prop_ks(results_dir, setting)
  if (length(detected) == 0L) {
    stop(sprintf(
      "No results matching '%d-*-*-p<K>.RData' found in %s; pass --prop_k=<spec> explicitly.",
      setting, results_dir
    ), call. = FALSE)
  }
  detected
}

# ---- output target ---------------------------------------------------
out_dir  <- "output/results"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, sprintf("scores%d-prop%s.RData", setting, data_suffix))

cat(sprintf(
  "scores-prop: setting=%d data=%s results_dir=%s out=%s\n  methods=%s datasets=%s prop_k=%s\n",
  setting, data_path, results_dir, out_file,
  paste(methods, collapse = ","),
  paste(range(datasets), collapse = ".."),
  paste(prop_ks, collapse = ",")
))

# ---- preallocate arrays ---------------------------------------------
nsets <- length(datasets)
nmeth <- length(methods)
nks   <- length(prop_ks)

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
                                     prop_k  = as.character(prop_ks)))

skew.methods <- c(2L, 4L)

obs <- c(rep(TRUE, nrow(y) - ntest), rep(FALSE, ntest))

# ---- score loop ------------------------------------------------------
save_checkpoint <- function() {
  save(quant.score, brier.score,
       beta.0, beta.1, beta.2,
       tau.alpha, tau.beta, rho, nu, gamma, lambda,
       elapsed_sec,
       probs, intervals, prop_ks, datasets, methods, setting,
       data_path, data_suffix, results_dir,
       file = out_file)
}

for (di in seq_along(datasets)) {
  set <- datasets[di]
  thresholds <- quantile(y[, , set, setting], probs = probs, na.rm = TRUE)
  validate   <- y[!obs, , set, setting]

  for (mi in seq_along(methods)) {
    method <- methods[mi]
    for (ki in seq_along(prop_ks)) {
      prop_k <- prop_ks[ki]
      f <- file.path(results_dir,
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

  if (set %% 10L == 0L) {
    save_checkpoint()
    cat("  -> checkpoint saved (", out_file, ")\n", sep = "")
  }
}

save_checkpoint()
cat("\nWrote ", out_file, "\n", sep = "")
