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
# per-block predictive samples + the validation truth + runtime info +
# per-block fit diagnostics and posterior summaries.
#
# The last two exist because the chunked driver deletes the fit files
# after scoring (DESKTOP-61SBCCI/expA_hn_driver.ps1): anything not summarised
# while `fit` is still in memory is unrecoverable without a full refit.
# scores.R carries both into the score cache, which is the durable copy.
#
# Usage:
#   Rscript run-settings.R --setting=<id> [--data=<path>] [--hn]
#                          [--iters=<n>] [--burn=<n>]
#                          [datasets] [workers] [methods]
# Examples:
#   Rscript run-settings.R --setting=3 "1"   1 "1:2"
#   Rscript run-settings.R --setting=3 "1:5" 4 "(1,2)"
#   Rscript run-settings.R --hn --setting=5 "1:5" 6 "c(1,2,4)"
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
source("./fit_diag_utils.R")

# -----------------------------
# User controls
# -----------------------------
iters <- 20000
burn <- 10000
update <- 1000
thin <- 1

# ---- CLI parsing ------------------------------------------------------
args_raw <- commandArgs(trailingOnly = TRUE)
# extract_leading_flags() demands a value after a bare "--flag", so pull the
# valueless form of --hn out first; --hn=TRUE still goes through the parser.
hn_bare <- any(args_raw == "--hn")
args_raw <- args_raw[args_raw != "--hn"]
parsed <- extract_leading_flags(args_raw,
  c("data", "setting", "hn", "iters", "burn"))
args <- parsed$args
data_flag <- parsed$values$data
setting_flag <- parsed$values$setting

# lambda ~ HN(0, lambda.s) instead of N(0, lambda.s): removes the discrete
# sign reflection of the (beta0, lambda, z) ridge at the model level, so no
# guard or reseed loop is needed. See tex/lambda_phiz_ridge and the block-1
# counterfactual in block1_positive_control/hn_prior_experiment/.
lambda_positive <- hn_bare ||
  (!is.null(parsed$values$hn) &&
     tolower(parsed$values$hn) %in% c("1", "t", "true", "yes"))

# iters/burn overrides exist for the structural smoke test only; production
# runs must leave them at the defaults above so every arm stays comparable.
if (!is.null(parsed$values$iters)) iters <- as.integer(parsed$values$iters)
if (!is.null(parsed$values$burn)) burn <- as.integer(parsed$values$burn)
if (is.na(iters) || is.na(burn) || burn >= iters) {
  stop("--iters/--burn must be integers with burn < iters", call. = FALSE)
}

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

prior_tag <- tbf_prior_tag(lambda_positive)
results_dir <- derive_results_dir(data_path, "results", prior_tag)
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

blocks <- tbf_blocks(block_seams, block_H, nt)

# Reference scale for the run banner, computed once on the master and
# exported so every worker reports against the same number. The ground
# truth that used to sit beside it went with the B assertion to
# block1_positive_control/fit_assertions.R.
marginal_sd <- sd(y[, , , setting])

cat(sprintf(
  "Data=%s setting=%d datasets=%s methods=%s blocks=%d (H=%d) iters=%d\n",
  data_path, setting, paste(dataset_ids, collapse = ","),
  paste(method_ids, collapse = ","), length(blocks), block_H, iters
))
cat(sprintf("lambda prior = %s | marginal sd = %.3f | warn = %s\n",
  if (lambda_positive) "HN(0, 20) [--hn]" else "N(0, 20) [default]",
  marginal_sd, format(getOption("warn"))
))

# ---- one (method, dataset) task --------------------------------------
run_method <- function(method_id, dataset_id) {
  spec <- catalog[catalog$method_id == method_id, , drop = FALSE]
  y.d <- y[, , dataset_id, setting]                   # ns x nt
  outputfile <- build_tbf_result_file(results_dir, setting, method_id, dataset_id)

  # Resume: a finished cell is never redone. The file is written once, at
  # the end of all blocks, so its presence means the whole cell is done.
  if (file.exists(outputfile)) {
    cat(sprintf("[Dataset %d][Method %d] exists, skip -> %s\n",
      dataset_id, method_id, outputfile))
    return(invisible(sprintf("s%d m%d d%d skip", setting, method_id, dataset_id)))
  }

  seed_used <- get_tbf_seed(setting, method_id, dataset_id)
  set.seed(seed_used)
  # Method 2 recurses on two lags, method 4 on one, method 1 draws from the
  # stationary marginal. Both flags are passed to forecast_block() AND
  # recorded, because an AR(1) fit forecast from the i.i.d. branch is the
  # kind of silent null that looks entirely plausible in a score table.
  ar2_used <- isTRUE(spec$ar2[1])
  ar1_used <- isTRUE(spec$temporal[1]) && !ar2_used
  cat(sprintf("[Dataset %d][Method %d: %s] start\n", dataset_id, method_id, spec$label[1]))
  started_at <- Sys.time()
  tic <- proc.time()

  block_out <- vector("list", length(blocks))
  diag_rows <- vector("list", length(blocks))
  post_rows <- vector("list", length(blocks))
  for (b in seq_along(blocks)) {
    blk <- blocks[[b]]
    y.train <- y.d[, blk$train_times, drop = FALSE]
    x.train <- x[, blk$train_times, , drop = FALSE]
    x.block <- x[, blk$test_times, , drop = FALSE]

    # fit on the prefix; every site observed, so no s.pred / x.pred.
    #
    # DO NOT pass z.init. The backend default (mcmc_ar2.R:249-250) is
    #   z.init <- 0.6745 / sqrt(tau.init)
    # i.e. the MEDIAN of HalfNormal(0, 1/sqrt(tau.init)), which maps to
    # z.star = 0 -- dead centre on the Gaussian copula scale. It is finite and
    # correct. (An earlier comment here claimed the default was 0 and that
    # hn.cop(0, sig) = -Inf made it fatal; the default is NOT 0, and the
    # "fix", z.init = 0.01, was the actual bug.)
    #
    # z.init = 0.01 starts z.star at hn.cop(0.01, 1) = -2.41, a 2.4-sigma tail.
    # The MH chain never escapes: z froze near 0.2 when the model wants ~1.6.
    # A frozen z looks like a near-constant series, so phi.z was driven to a
    # near-unit root (1.11 vs a truth of 0.80), which pinned z harder still --
    # a self-reinforcing freeze. Since the data only see the product lambda*z,
    # lambda then compensated, running to -98 (truth +3), and beta0 to 36
    # (truth 10), all while reproducing the data mean so the likelihood looked
    # fine from inside the fit. forecast_block correctly redraws z at its proper
    # scale (~1.6), multiplies by the corrupted lambda, and the predictive
    # explodes: SD 104 against a marginal SD of 4.38, median -32 against a data
    # centre of +14. That is the entire "AR(2) is 4x worse than i.i.d." result.
    #
    # Verified: dropping this one argument restores lambda to +2.8/+2.9,
    # beta0 to 10.4, z to 1.68, phi.z to 0.62, and the lead-15 predictive SD to
    # 4.9 (vs a marginal 4.38) on the exact chain that previously blew up.
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
      ar2_w = ar2_used,
      ar2_tau = ar2_used,
      ar2_z = ar2_used,
      rho.upper = 15, nu.upper = 10,
      lambda.positive = lambda_positive
    )

    # ar1 = TRUE is what routes method 4 into the single-lag recursion.
    # Without it a temporal AR(1) fit silently forecasts from the i.i.d.
    # branch (time_block_helpers.R:304), i.e. no temporal signal at all.
    yhat <- forecast_block(
      fit = fit, seam = blk$seam, H = block_H,
      x_block = x.block, s = s, ar2 = ar2_used, ar1 = ar1_used
    )

    block_out[[b]] <- list(
      block_id = b,
      seam = blk$seam,
      leads = blk$leads,
      test_times = blk$test_times,
      yhat = yhat,                       # iters x ns x H predictive sample
      y_val = y.d[, blk$test_times]      # ns x H validation truth
    )

    # ---- everything below must happen while `fit` is still alive ------
    # The fit is discarded here and the whole file is deleted after
    # scoring, so an unsummarised quantity needs a full refit to recover.
    # Numeric summaries only. The A/A'/B/C pass flags predate the HN
    # prior, which removed the ridge they guarded against; this study has
    # not recorded them since commit bcc2d39 and never gated on them, so
    # since 2026-08-05 it does not compute them either -- the assertions
    # live in block1_positive_control/fit_assertions.R, the study that
    # gates on B and C.
    chk <- fit_diag_summary(fit, yhat, data_mean = mean(y.train))
    diag_rows[[b]] <- cbind(
      setting = setting, method = method_id, dataset = dataset_id,
      block = b, seam = blk$seam, seed = seed_used,
      hn = lambda_positive, ar2_used = ar2_used, ar1_used = ar1_used,
      chk
    )
    post_rows[[b]] <- cbind(
      setting = setting, method = method_id, dataset = dataset_id,
      block = b, tbf_posterior_summary(fit)
    )

    rm(fit)
    gc()
    cat(sprintf("  block %d (seam %d) forecast done  lambda %+.2f  SDmax %.2f\n",
      b, blk$seam, chk$lambda, chk$sd_lead_max))
  }

  # fit_diag and post_summary travel inside the fit file only; scores.R
  # carries them into the score cache (fit.diag / post.summary), which is
  # what survives the fit deletion. The per-cell side files this used to
  # write (output/diag/{diag,post,pred}_S-m-d.*) were a second copy of
  # the same rows plus a predictive summary nothing consumed; removed
  # 2026-08-03.
  fit_diag <- do.call(rbind, diag_rows)
  post_summary <- do.call(rbind, post_rows)

  elapsed_sec <- unname((proc.time() - tic)[3])
  runtime_info <- build_tbf_runtime_info(
    started_at = started_at, elapsed_sec = elapsed_sec, outputfile = outputfile,
    setting_id = setting, dataset_id = dataset_id,
    runner = "run-settings.R::run_method",
    method_id = method_id, method_key = as.character(method_id),
    control = list(iters = iters, burn = burn, update = update, thin = thin,
                   block_H = block_H, block_seams = block_seams,
                   lambda_positive = lambda_positive,
                   ar2_used = ar2_used, ar1_used = ar1_used,
                   seed = seed_used,
                   warn_option = getOption("warn"),
                   host = unname(Sys.info()[["nodename"]]),
                   r_version = as.character(getRversion()),
                   git_sha = tbf_git_sha())
  )
  forecast <- list(
    setting = setting, method_id = method_id, dataset_id = dataset_id,
    blocks = block_out
  )
  save(forecast, runtime_info, fit_diag, post_summary,
    file = outputfile)
  rm(forecast, block_out, runtime_info)
  gc()
  cat(sprintf("[Dataset %d][Method %d] done -> %s\n", dataset_id, method_id, outputfile))
  invisible(sprintf("s%d m%d d%d ok (%.0f min)",
    setting, method_id, dataset_id, elapsed_sec / 60))
}

# One failing cell must not take down the other workers in flight: PSOCK
# propagates an error out of parLapply and the whole batch dies with it.
# Cells that already saved survive; this one is reported and re-run by the
# driver's next inventory round.
run_method_safe <- function(method_id, dataset_id) {
  tryCatch(run_method(method_id, dataset_id), error = function(e) {
    msg <- sprintf("s%d m%d d%d FAILED: %s",
      setting, method_id, dataset_id, conditionMessage(e))
    cat(msg, "\n")
    invisible(msg)
  })
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
    source("./fit_diag_utils.R")
    NULL
  })
  parallel::clusterExport(cl, c("data_path"), envir = environment())
  parallel::clusterEvalQ(cl, {
    load(data_path, envir = .GlobalEnv)
    NULL
  })
  # PSOCK resolves free names against the WORKER's globals, so a global
  # missing from this list either errors or -- worse -- silently binds to
  # something else on the worker. Every name run_method() touches goes here.
  parallel::clusterExport(cl, c(
    "setting", "iters", "burn", "update", "thin", "results_dir",
    "catalog", "blocks", "block_H", "block_seams", "run_method", "tasks",
    "lambda_positive", "marginal_sd", "run_method_safe"
  ), envir = environment())
  out <- parallel::parLapply(cl, seq_len(nrow(tasks)), function(i) {
    run_method_safe(tasks$method[i], tasks$dataset[i])
  })
  parallel::stopCluster(cl)
  cat(paste(unlist(out), collapse = "\n"), "\n")
} else {
  for (i in seq_len(nrow(tasks))) {
    tic <- proc.time()
    run_method_safe(tasks$method[i], tasks$dataset[i])
    cat(sprintf("elapsed (dataset %d, method %d): %.2f sec\n",
      tasks$dataset[i], tasks$method[i], (proc.time() - tic)[3]))
  }
}

cat("All specified methods and datasets finished.\n")
