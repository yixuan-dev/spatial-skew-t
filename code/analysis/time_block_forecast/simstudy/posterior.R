#########################################################################
# posterior.R -- posterior mean and SD summary for time-block fits.
#
# Mirrors code/analysis/simstudy/posterior.R, with two structural
# differences forced by this study's fit files:
#
#   - run-settings.R discards the fit object (rm(fit)) before saving, so
#     the raw parameter draws do not exist anywhere. The finest available
#     resolution is the per-block post_summary data.frame (mean, sd,
#     q025, q50, q975 per parameter) that was computed while the fit was
#     alive. This script aggregates that; it cannot recompute anything.
#   - fits are per-block (5 seams per cell), so the block axis replaces
#     the sibling's mrts_k axis: arrays are [dataset, method, block].
#
# Loads fits from <results_dir>/<setting>-<method>-<dataset>.RData and
# writes:
#
#   output/results/posterior<setting><suffix>.RData      (per-fit arrays)
#   output/tables/posterior_summary<setting><suffix>.csv (aggregated)
#
# Note the fits are ~800 MB each because of the predictive draws; the
# post_summary rows ride along at ~0 MB but load() has to swallow the
# whole file, so a full 60-cell campaign takes ~10 minutes. The same
# rows also reach the score cache as the post.summary array; this script
# reads the fits deliberately, as the authoritative source.
#
# Usage:
#   Rscript posterior.R --setting=<id>
#                       [--data=<path>]
#                       [--methods=<spec>]   default = auto-scan results_dir
#                       [--datasets=<spec>]  default = auto-scan results_dir
#
# Examples:
#   Rscript posterior.R --setting=7
#   Rscript posterior.R --setting=5 --methods="c(1,2,4)" --datasets="1:10"
#########################################################################

rm(list = ls())

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0L) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]),
    winslash = "/", mustWork = FALSE)
  if (dir.exists(dirname(script_path))) setwd(dirname(script_path))
}

source("./time_block_helpers.R")

# ---- CLI parsing ------------------------------------------------------
cli_args <- commandArgs(trailingOnly = TRUE)
prior <- tbf_take_prior_flag(cli_args,
  c("data", "setting", "methods", "datasets"))
parsed <- prior$parsed
flags <- parsed$values
prior_tag <- prior$prior_tag

if (is.null(flags$setting) || !nzchar(flags$setting)) {
  stop("posterior.R: --setting=<id> is required.", call. = FALSE)
}

data_path <- resolve_simstudy_data_path(flags$data)
data_suffix <- derive_data_suffix(data_path)

setting_ids <- parse_index_expr(flags$setting, "setting")
if (length(setting_ids) != 1L) {
  stop("posterior.R: --setting must be a single integer.", call. = FALSE)
}
setting <- as.integer(setting_ids)

results_dir <- derive_results_dir(data_path, "results", prior_tag)
if (!dir.exists(results_dir)) {
  stop(sprintf("results directory not found: %s (run run-settings.R --prior=%s first)",
    results_dir, prior_tag), call. = FALSE)
}

# ---- methods / datasets (auto-scan from results_dir) ------------------
scan_fit_files <- function(results_dir, setting) {
  files <- list.files(results_dir)
  mm <- regmatches(files,
    regexec(sprintf("^%d-(\\d+)-(\\d+)\\.RData$", setting), files))
  hits <- Filter(function(x) length(x) == 3L, mm)
  data.frame(method  = vapply(hits, function(x) as.integer(x[2]), integer(1)),
             dataset = vapply(hits, function(x) as.integer(x[3]), integer(1)))
}

found <- scan_fit_files(results_dir, setting)

methods <- if (!is.null(flags$methods) && nzchar(flags$methods)) {
  sort(unique(as.integer(parse_index_expr(flags$methods, "methods"))))
} else {
  if (nrow(found) == 0L) {
    stop(sprintf("No fits for setting %d found in %s; pass --methods explicitly.",
      setting, results_dir), call. = FALSE)
  }
  sort(unique(found$method))
}

datasets <- if (!is.null(flags$datasets) && nzchar(flags$datasets)) {
  sort(unique(as.integer(parse_index_expr(flags$datasets, "datasets"))))
} else {
  keep <- found$method %in% methods
  if (!any(keep)) {
    stop(sprintf("No fits for setting %d methods [%s] in %s.",
      setting, paste(methods, collapse = ","), results_dir), call. = FALSE)
  }
  sort(unique(found$dataset[keep]))
}

# ---- output paths -----------------------------------------------------
for (d in c("output/results", "output/tables")) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
out_rdata <- tbf_score_cache_file(setting, data_suffix, prior_tag,
  dir = "output/results", prefix = "posterior")
out_csv <- file.path("output/tables",
  sprintf("posterior_summary%d%s_%s.csv", setting, data_suffix, prior_tag))

cat(sprintf(
  "posterior: setting=%d  results_dir=%s\n  methods=%s  datasets=%s\n  out=%s\n",
  setting, results_dir,
  paste(methods, collapse = ","),
  paste(range(datasets), collapse = ".."),
  out_rdata))

# ---- main loop: harvest post_summary from each fit --------------------
# Accumulate the long-form rows first; the arrays are built afterwards,
# so nothing needs the block count or parameter set up front, and a
# parameter absent from a method (e.g. the phis under i.i.d.) is NA in
# the arrays by construction rather than by an applies-registry.
long_rows <- vector("list", length(datasets) * length(methods))
iters_detected <- NA_integer_

for (di in seq_along(datasets)) {
  set <- datasets[di]
  for (mi in seq_along(methods)) {
    method <- methods[mi]
    f <- build_tbf_result_file(results_dir, setting, method, set)
    if (!file.exists(f)) { cat("missing: ", f, "\n", sep = ""); next }

    env <- new.env(parent = emptyenv())
    load(f, envir = env)
    ps <- env$post_summary
    if (is.null(ps) || !is.data.frame(ps) || nrow(ps) == 0L) {
      cat("no post_summary: ", f, "\n", sep = "")
      rm(env); gc(FALSE); next
    }

    if (is.na(iters_detected) &&
          !is.null(env$runtime_info$control$iters)) {
      iters_detected <- as.integer(env$runtime_info$control$iters)
      cat(sprintf("  (iters detected: %d)\n", iters_detected))
    }

    long_rows[[(di - 1L) * length(methods) + mi]] <-
      ps[, c("method", "dataset", "block", "param", "mean", "sd")]

    # the fit is ~800 MB of predictive draws; drop it before the next one.
    rm(ps, env)
    gc(FALSE)
    cat(sprintf("dataset %d  method %d done\n", set, method))
  }
}

long <- do.call(rbind, Filter(Negate(is.null), long_rows))
if (is.null(long) || nrow(long) == 0L) {
  stop("No post_summary rows harvested -- check that fits exist.", call. = FALSE)
}

# ---- build [dataset, method, block] arrays ----------------------------
blocks_seen <- sort(unique(as.integer(long$block)))
param_names <- intersect(tbf_post_params(), unique(long$param))
extra <- setdiff(unique(long$param), tbf_post_params())
if (length(extra) > 0L) param_names <- c(param_names, sort(extra))

dn <- list(dataset = as.character(datasets),
           method  = as.character(methods),
           block   = as.character(blocks_seen))
make_arr <- function() {
  array(NA_real_,
    dim = c(length(datasets), length(methods), length(blocks_seen)),
    dimnames = dn)
}
post_mean <- setNames(lapply(param_names, function(p) make_arr()), param_names)
post_sd   <- setNames(lapply(param_names, function(p) make_arr()), param_names)

for (pname in param_names) {
  rows <- long[long$param == pname, , drop = FALSE]
  ii <- cbind(match(as.integer(rows$dataset), datasets),
              match(as.integer(rows$method), methods),
              match(as.integer(rows$block), blocks_seen))
  post_mean[[pname]][ii] <- rows$mean
  post_sd[[pname]][ii]   <- rows$sd
}

# ---- save RData -------------------------------------------------------
fits_dir <- results_dir
save(post_mean, post_sd, param_names,
     datasets, methods, blocks_seen, setting, data_suffix,
     iters_detected, fits_dir,
     file = out_rdata)
cat("\nWrote ", out_rdata, "\n", sep = "")

# ---- build and write summary CSV --------------------------------------
# One row per (method, block, parameter), aggregated over datasets --
# the sibling's layout with block in place of mrts_k.
rows <- list()
for (pname in param_names) {
  pm <- post_mean[[pname]]
  ps <- post_sd[[pname]]
  for (mi in seq_along(methods)) {
    for (bi in seq_along(blocks_seen)) {
      mn_vals <- pm[, mi, bi]
      sd_vals <- ps[, mi, bi]
      ok <- is.finite(mn_vals) & is.finite(sd_vals)
      if (!any(ok)) next
      rows[[length(rows) + 1L]] <- data.frame(
        setting          = setting,
        method           = methods[mi],
        block            = blocks_seen[bi],
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
  cat("No results to summarise -- check that result files exist.\n")
}
