#########################################################################
# run-settings.R - fit + forecast driver for the time-block forecast
# simulation study (Stage A of the pipeline).
#
# For each (setting, method, dataset) and each held-out time block:
#   1. fit the spatial skew-t model on the contiguous prefix y[, 1:T_o]
#      (all sites observed -- no spatial holdout);
#   2. forecast the window (T_o, T_o + H] with forecast_block(), i.e.
#      Algorithm 1 of tex/time_block_strategy/ar2_rethink.tex.
#
# The fit conditions the forecast on the posterior of the seam state,
# so the AR(2) and i.i.d. methods are no longer structurally
# exchangeable (contrast Proposition 2).
#
# Output: results/<setting>-<method>-<dataset>.RData, holding a list of
# per-block predictive samples + the validation truth + runtime info.
#
# Usage:
#   Rscript run-settings.R --setting=<id> [--data=<path>]
#                          [datasets] [workers] [methods]
# Examples:
#   Rscript run-settings.R --setting=3 "1"   1 "1:2"
#   Rscript run-settings.R --setting=3 "1:5" 4 "(1,2)"
#########################################################################

# ar2_load.R does rm(list = ls()); source it first, then re-source the
# time-block helpers it wiped. Per the project instruction we reuse the
# Morris study's loader rather than keeping a private copy.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1]),
    winslash = "/", mustWork = FALSE
  ))
  if (dir.exists(script_dir)) setwd(script_dir)
}

source("../../simstudy/ar2_load.R", chdir = TRUE)
source("./time_block_helpers.R")

# -----------------------------
# User controls
# -----------------------------
iters <- 20000
burn <- 10000
update <- 1000
thin <- 1

# ---- CLI parsing ------------------------------------------------------
args_raw <- commandArgs(trailingOnly = TRUE)
parsed <- extract_leading_flags(args_raw, c("data", "setting"))
args <- parsed$args
data_flag <- parsed$values$data
setting_flag <- parsed$values$setting

data_path <- if (!is.null(data_flag) && nzchar(data_flag)) data_flag else "./simdata.RData"
if (!file.exists(data_path)) {
  stop(sprintf("Data file not found: %s (run setup.R first)", data_path), call. = FALSE)
}
# ar2_load.R has loaded the Morris ./simdata.RData; reload the time-block
# dataset over it so y/x/s/block_* refer to this study.
load(data_path, envir = .GlobalEnv)

if (is.null(setting_flag) || !nzchar(setting_flag)) {
  stop("run-settings.R: --setting=<id> is required.", call. = FALSE)
}
setting <- as.integer(parse_index_expr(setting_flag, "setting"))
if (length(setting) != 1L || setting < 1L || setting > dim(y)[4]) {
  stop(sprintf("--setting must be a single integer in 1..%d", dim(y)[4]), call. = FALSE)
}

datasets_spec <- if (length(args) >= 1) args[1] else "1"
workers <- if (length(args) >= 2) as.integer(args[2]) else 1L
methods_spec <- if (length(args) >= 3) args[3] else "1:2"

dataset_ids <- parse_index_expr(datasets_spec, "datasets")
if (any(dataset_ids < 1 | dataset_ids > dim(y)[3])) {
  stop(sprintf("datasets must be in 1..%d", dim(y)[3]), call. = FALSE)
}
method_ids <- parse_index_expr(methods_spec, "methods")
catalog <- get_tbf_method_catalog()
if (any(!method_ids %in% catalog$method_id)) {
  stop(sprintf("methods must be in 1..%d", max(catalog$method_id)), call. = FALSE)
}
if (is.na(workers) || workers < 1) stop("workers must be a positive integer", call. = FALSE)

results_dir <- derive_results_dir(data_path, "results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

blocks <- tbf_blocks(block_seams, block_H, nt)

cat(sprintf(
  "Data=%s setting=%d datasets=%s methods=%s blocks=%d (H=%d) iters=%d\n",
  data_path, setting, paste(dataset_ids, collapse = ","),
  paste(method_ids, collapse = ","), length(blocks), block_H, iters
))

# ---- one (method, dataset) task --------------------------------------
run_method <- function(method_id, dataset_id) {
  spec <- catalog[catalog$method_id == method_id, , drop = FALSE]
  y.d <- y[, , dataset_id, setting]                   # ns x nt
  outputfile <- build_tbf_result_file(results_dir, setting, method_id, dataset_id)

  set.seed(get_tbf_seed(setting, method_id, dataset_id))
  cat(sprintf("[Dataset %d][Method %d: %s] start\n", dataset_id, method_id, spec$label[1]))
  started_at <- Sys.time()
  tic <- proc.time()

  block_out <- vector("list", length(blocks))
  for (b in seq_along(blocks)) {
    blk <- blocks[[b]]
    y.train <- y.d[, blk$train_times, drop = FALSE]
    x.train <- x[, blk$train_times, , drop = FALSE]
    x.block <- x[, blk$test_times, , drop = FALSE]

    # fit on the prefix; every site observed, so no s.pred / x.pred.
    fit <- mcmc(
      y = y.train, x = x.train, s = s,
      method = "t", skew = isTRUE(spec$skew[1]),
      thresh.all = 0, thresh.quant = TRUE,
      nknots = spec$nknots[1],
      iterplot = FALSE, iters = iters, burn = burn, update = update,
      min.s = c(0, 0), max.s = c(10, 10),
      temporalw = isTRUE(spec$temporal[1]),
      temporaltau = isTRUE(spec$temporal[1]),
      temporalz = isTRUE(spec$temporal[1]),
      ar2_w = isTRUE(spec$ar2[1]),
      ar2_tau = isTRUE(spec$ar2[1]),
      ar2_z = isTRUE(spec$ar2[1]),
      rho.upper = 15, nu.upper = 10
    )

    yhat <- forecast_block(
      fit = fit, seam = blk$seam, H = block_H,
      x_block = x.block, s = s, ar2 = isTRUE(spec$ar2[1])
    )

    block_out[[b]] <- list(
      block_id = b,
      seam = blk$seam,
      leads = blk$leads,
      test_times = blk$test_times,
      yhat = yhat,                       # iters x ns x H predictive sample
      y_val = y.d[, blk$test_times]      # ns x H validation truth
    )
    rm(fit)
    gc()
    cat(sprintf("  block %d (seam %d) forecast done\n", b, blk$seam))
  }

  elapsed_sec <- unname((proc.time() - tic)[3])
  runtime_info <- build_tbf_runtime_info(
    started_at = started_at, elapsed_sec = elapsed_sec, outputfile = outputfile,
    setting_id = setting, dataset_id = dataset_id,
    runner = "run-settings.R::run_method",
    method_id = method_id, method_key = as.character(method_id),
    control = list(iters = iters, burn = burn, update = update, thin = thin,
                   block_H = block_H, block_seams = block_seams)
  )
  forecast <- list(
    setting = setting, method_id = method_id, dataset_id = dataset_id,
    blocks = block_out
  )
  save(forecast, runtime_info, file = outputfile)
  rm(forecast, block_out, runtime_info)
  gc()
  cat(sprintf("[Dataset %d][Method %d] done -> %s\n", dataset_id, method_id, outputfile))
}

# ---- task grid --------------------------------------------------------
tasks <- expand.grid(
  dataset = dataset_ids, method = method_ids,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
workers_use <- min(workers, nrow(tasks))

if (workers_use > 1) {
  cat(sprintf("Running %d tasks across %d workers...\n", nrow(tasks), workers_use))
  cl <- parallel::makeCluster(workers_use, type = "PSOCK")
  wd <- getwd()
  parallel::clusterExport(cl, "wd", envir = environment())
  parallel::clusterEvalQ(cl, {
    setwd(wd)
    source("../../simstudy/ar2_load.R", chdir = TRUE)
    source("./time_block_helpers.R")
    NULL
  })
  parallel::clusterExport(cl, c("data_path"), envir = environment())
  parallel::clusterEvalQ(cl, {
    load(data_path, envir = .GlobalEnv)
    NULL
  })
  parallel::clusterExport(cl, c(
    "setting", "iters", "burn", "update", "thin", "results_dir",
    "catalog", "blocks", "block_H", "block_seams", "run_method"
  ), envir = environment())
  parallel::parLapply(cl, seq_len(nrow(tasks)), function(i) {
    run_method(tasks$method[i], tasks$dataset[i])
  })
  parallel::stopCluster(cl)
} else {
  for (i in seq_len(nrow(tasks))) {
    tic <- proc.time()
    run_method(tasks$method[i], tasks$dataset[i])
    cat(sprintf("elapsed (dataset %d, method %d): %.2f sec\n",
      tasks$dataset[i], tasks$method[i], (proc.time() - tic)[3]))
  }
}

cat("All specified methods and datasets finished.\n")
