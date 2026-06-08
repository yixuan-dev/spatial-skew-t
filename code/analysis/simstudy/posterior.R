#########################################################################
# posterior.R - Posterior mean and SD summary for simstudy fits.
#
# Loads fits from <results_dir>/<setting>-<method>-<dataset>[-K<mrts_k>].RData,
# computes posterior mean and SD for each scalar parameter, and writes:
#
#   output/results/posterior<setting><suffix>.RData   (per-fit arrays)
#   output/tables/posterior_summary<setting><suffix>.csv  (aggregated)
#
# Mirrors scores.R structure: same CLI interface, same path resolution,
# same file-loading pattern.
#
# Usage:
#   Rscript posterior.R --setting=<id>
#                       [--data=<path>]
#                       [--methods=<spec>]   default 1:10
#                       [--datasets=<spec>]  default = auto-scan from results_dir
#                       [--mrts_k=<spec>]    default = auto-detect from results_dir
#
# Examples:
#   Rscript posterior.R --setting=9 --methods="7:8"
#   Rscript posterior.R --setting=10 --methods="c(7,8)" --datasets="1:10"
#########################################################################

rm(list = ls())

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0L) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]),
                               winslash = "/", mustWork = FALSE)
  if (dir.exists(dirname(script_path))) setwd(dirname(script_path))
}

source("./mrts_cov_helpers.R")

# ---- CLI parsing ---------------------------------------------------------
cli_args <- commandArgs(trailingOnly = TRUE)
parsed   <- extract_leading_flags(cli_args, c("data", "setting", "methods", "datasets", "mrts_k"))
flags    <- parsed$values

if (is.null(flags$setting) || !nzchar(flags$setting))
  stop("posterior.R: --setting=<id> is required.", call. = FALSE)

# ---- resolve results_dir / data_suffix -----------------------------------
data_path   <- resolve_simstudy_data_path(flags$data)
data_suffix <- derive_data_suffix(data_path)

results_dir_env <- trimws(Sys.getenv("SIMSTUDY_RESULTS_DIR", unset = ""))
if (!nzchar(results_dir_env))
  results_dir_env <- trimws(Sys.getenv("SIMSTUDY_MRTS_RESULTS_DIR", unset = ""))
results_dir <- if (nzchar(results_dir_env)) {
  results_dir_env
} else {
  derive_results_dir(data_path, "results")
}
if (!dir.exists(results_dir))
  stop(sprintf("results directory not found: %s", results_dir), call. = FALSE)

# ---- setting -------------------------------------------------------------
setting_ids <- parse_index_expr(flags$setting, "setting")
if (length(setting_ids) != 1L)
  stop("posterior.R: --setting must be a single integer.", call. = FALSE)
setting <- as.integer(setting_ids)

# ---- methods -------------------------------------------------------------
methods <- if (!is.null(flags$methods) && nzchar(flags$methods)) {
  m <- parse_index_expr(flags$methods, "methods")
  if (any(m < 1L | m > 10L)) stop("methods must be in 1..10", call. = FALSE)
  sort(unique(as.integer(m)))
} else {
  1:10
}

# ---- mrts_k (same auto-detect as scores.R) -------------------------------
auto_detect_mrts_ks <- function(results_dir, setting) {
  files    <- list.files(results_dir)
  base_pat <- sprintf("^%d-\\d+-\\d+\\.RData$", setting)
  k_pat    <- sprintf("^%d-\\d+-\\d+-K(\\d+)\\.RData$", setting)
  ks       <- integer(0)
  if (any(grepl(base_pat, files))) ks <- c(ks, 0L)
  mm    <- regmatches(files, regexec(k_pat, files))
  kvals <- vapply(mm, function(x) if (length(x) >= 2L) as.integer(x[2]) else NA_integer_, integer(1))
  c(ks, kvals[!is.na(kvals)]) |> unique() |> sort()
}

mrts_ks <- if (!is.null(flags$mrts_k) && nzchar(flags$mrts_k)) {
  sort(unique(as.integer(parse_index_expr(flags$mrts_k, "mrts_k"))))
} else {
  detected <- auto_detect_mrts_ks(results_dir, setting)
  if (length(detected) == 0L)
    stop(sprintf("No fits for setting %d found in %s; pass --mrts_k=<spec> explicitly.",
                 setting, results_dir), call. = FALSE)
  detected
}

# ---- datasets (auto-scan from results_dir) --------------------------------
auto_detect_datasets <- function(results_dir, setting, methods, mrts_ks) {
  files <- list.files(results_dir)
  meth_re <- paste(methods, collapse = "|")
  ds <- integer(0)
  for (mk in mrts_ks) {
    pat <- if (mk == 0L) {
      sprintf("^%d-(%s)-(\\d+)\\.RData$", setting, meth_re)
    } else {
      sprintf("^%d-(%s)-(\\d+)-K%d\\.RData$", setting, meth_re, mk)
    }
    mm    <- regmatches(files, regexec(pat, files))
    found <- vapply(mm, function(x) if (length(x) >= 3L) as.integer(x[3]) else NA_integer_, integer(1))
    ds    <- c(ds, found[!is.na(found)])
  }
  sort(unique(ds))
}

datasets <- if (!is.null(flags$datasets) && nzchar(flags$datasets)) {
  as.integer(parse_index_expr(flags$datasets, "datasets"))
} else {
  found <- auto_detect_datasets(results_dir, setting, methods, mrts_ks)
  if (length(found) == 0L)
    stop(sprintf("No dataset files found for setting %d (methods %s) in %s.",
                 setting, paste(methods, collapse = ","), results_dir), call. = FALSE)
  found
}

# ---- output paths --------------------------------------------------------
for (d in c("output/results", "output/tables"))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

out_rdata <- file.path("output/results", sprintf("posterior%d%s.RData",   setting, data_suffix))
out_csv   <- file.path("output/tables",  sprintf("posterior_summary%d%s.csv", setting, data_suffix))

cat(sprintf(
  "posterior: setting=%d  results_dir=%s\n  methods=%s  datasets=%s  mrts_k=%s\n  out=%s\n",
  setting, results_dir,
  paste(methods, collapse = ","),
  paste(range(datasets), collapse = ".."),
  paste(mrts_ks, collapse = ","),
  out_rdata
))

# ---- parameter registry --------------------------------------------------
skew_methods <- c(2L, 4L, 7L, 8L, 9L, 10L)
ar2_methods  <- c(7L, 8L)
ar1_methods  <- c(9L, 10L)

params <- list(
  list(name = "beta.0",    extract = function(fit) fit$beta[, 1],    applies = NULL),
  list(name = "beta.1",    extract = function(fit) fit$beta[, 2],    applies = NULL),
  list(name = "beta.2",    extract = function(fit) fit$beta[, 3],    applies = NULL),
  list(name = "tau.alpha", extract = function(fit) fit$tau.alpha,    applies = NULL),
  list(name = "tau.beta",  extract = function(fit) fit$tau.beta,     applies = NULL),
  list(name = "rho",       extract = function(fit) fit$rho,          applies = NULL),
  list(name = "nu",        extract = function(fit) fit$nu,           applies = NULL),
  list(name = "gamma",     extract = function(fit) fit$gamma,        applies = NULL),
  list(name = "lambda",    extract = function(fit) fit$lambda,       applies = skew_methods),
  list(name = "phi1.tau",  extract = function(fit) fit$phi.tau[, 1], applies = c(ar2_methods, ar1_methods)),
  list(name = "phi2.tau",  extract = function(fit) fit$phi.tau[, 2], applies = ar2_methods),
  list(name = "phi1.z",    extract = function(fit) fit$phi.z[, 1],   applies = c(ar2_methods, ar1_methods)),
  list(name = "phi2.z",    extract = function(fit) fit$phi.z[, 2],   applies = ar2_methods),
  list(name = "phi1.w",    extract = function(fit) fit$phi.w[, 1],   applies = c(ar2_methods, ar1_methods)),
  list(name = "phi2.w",    extract = function(fit) fit$phi.w[, 2],   applies = ar2_methods)
)
param_names <- vapply(params, `[[`, character(1), "name")

# ---- preallocate result arrays -------------------------------------------
nsets <- length(datasets)
nmeth <- length(methods)
nks   <- length(mrts_ks)
dn    <- list(dataset = as.character(datasets),
              method  = as.character(methods),
              mrts_k  = as.character(mrts_ks))

make_arr <- function() array(NA_real_, dim = c(nsets, nmeth, nks), dimnames = dn)
post_mean <- setNames(lapply(params, function(p) make_arr()), param_names)
post_sd   <- setNames(lapply(params, function(p) make_arr()), param_names)

niter_detected <- NA_integer_

# ---- main loop -----------------------------------------------------------
for (di in seq_along(datasets)) {
  set <- datasets[di]
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
      if (!file.exists(f)) { cat("missing: ", f, "\n", sep = ""); next }

      env <- new.env(parent = emptyenv())
      load(f, envir = env)
      fit <- env$fit.1
      if (is.null(fit)) { cat("no fit.1: ", f, "\n", sep = ""); rm(env); next }

      if (is.na(niter_detected) && !is.null(fit$rho)) {
        niter_detected <- length(fit$rho)
        cat(sprintf("  (niter detected: %d)\n", niter_detected))
      }

      for (pi in seq_along(params)) {
        p <- params[[pi]]
        if (!is.null(p$applies) && !(method %in% p$applies)) next
        draws <- tryCatch(p$extract(fit), error = function(e) NULL)
        if (is.null(draws) || !is.numeric(draws) || length(draws) == 0L) next
        post_mean[[pi]][di, mi, ki] <- mean(draws, na.rm = TRUE)
        post_sd[[pi]][di, mi, ki]   <- sd(draws,   na.rm = TRUE)
      }

      rm(fit, env)
    }
    cat(sprintf("dataset %d  method %d done\n", set, method))
  }
}

# ---- save RData ----------------------------------------------------------
fits_dir <- results_dir
save(post_mean, post_sd, params, param_names,
     datasets, methods, mrts_ks, setting, data_suffix, niter_detected, fits_dir,
     file = out_rdata)
cat("\nWrote ", out_rdata, "\n", sep = "")

# ---- build and write summary CSV -----------------------------------------
rows <- list()
for (pi in seq_along(params)) {
  pname <- param_names[pi]
  pm    <- post_mean[[pi]]
  ps    <- post_sd[[pi]]
  for (mi in seq_along(methods)) {
    method <- methods[mi]
    for (ki in seq_along(mrts_ks)) {
      mrts_k  <- mrts_ks[ki]
      mn_vals <- pm[, mi, ki]
      sd_vals <- ps[, mi, ki]
      ok      <- is.finite(mn_vals) & is.finite(sd_vals)
      if (!any(ok)) next
      rows[[length(rows) + 1L]] <- data.frame(
        setting          = setting,
        method           = method,
        mrts_k           = mrts_k,
        parameter        = pname,
        n_datasets       = sum(ok),
        mean_post_mean   = mean(mn_vals[ok]),
        median_post_mean = median(mn_vals[ok]),
        sd_post_mean     = if (sum(ok) > 1L) sd(mn_vals[ok]) else NA_real_,
        mean_post_sd     = mean(sd_vals[ok]),
        stringsAsFactors = FALSE
      )
    }
  }
}

if (length(rows) > 0L) {
  summary_df <- do.call(rbind, rows)
  write.csv(summary_df, out_csv, row.names = FALSE)
  cat("Wrote ", out_csv, "\n", sep = "")
} else {
  cat("No results to summarise — check that result files exist.\n")
}
