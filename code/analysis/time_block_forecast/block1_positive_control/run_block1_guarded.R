#########################################################################
# run_block1_guarded.R -- block-1 positive control for AR(2) vs i.i.d.
#
# Question: fitting AR(2)-generated data with the AR(2) model, does the
# time-block forecast Brier score beat the i.i.d.-in-time skew-t baseline?
#
# Design (see plan): settings {4 strong, 5 nearunit} x methods {1 iid,
# 2 AR(2)} x datasets 1..10, block 1 ONLY (train y[, 1:50], forecast
# leads 1..15, i.e. t = 51..65). Data: ../simstudy/simdata.RData
# (setup.R with setting 5 = phi (0.15, 0.80), spectral radius 0.973).
#
# GUARD (pre-registered, symmetric across methods): block 1's short prefix
# is where the reflected lambda ridge bites (tex/lambda_phiz_ridge). Each
# fit must pass B (truth recovery: lambda sign+size, beta0) AND C
# (predictive spread <= 3 x marginal SD) from check_fit_consistency();
# otherwise refit with seed + 7919*attempt, at most 3 reseeds. If still
# failing, the last fit is kept and FLAGGED (pass = FALSE): the analysis
# excludes it from the primary comparison and reports it in a sensitivity
# run. A and A' (sdz_ratio, transverse to the ridge) are recorded always.
#
#   $R = "C:\Program Files\R\R-4.5.1\bin\Rscript.exe"
#   & $R run_block1_guarded.R --settings=4,5 --datasets=1:10 --methods=1:2
#########################################################################

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1]),
    winslash = "/", mustWork = FALSE
  ))
  if (dir.exists(script_dir)) setwd(script_dir)
}

# ar2_load.R does rm(list = ls()): source it FIRST, then the helpers it wiped.
source("../../simstudy/ar2_load.R", chdir = TRUE)
source("../simstudy/time_block_helpers.R")
source("./fit_assertions.R")
load("../simstudy/simdata.RData", envir = .GlobalEnv)
options(warn = 1)

# ---- CLI ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
getflag <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
settings <- eval(parse(text = paste0("c(", getflag("settings", "4,5"), ")")))
datasets <- eval(parse(text = getflag("datasets", "1:10")))
method_ids <- eval(parse(text = getflag("methods", "1:2")))
workers <- as.integer(getflag("workers", "1"))

iters <- 20000
burn <- 10000
update <- 2000
max_attempts <- 4L                     # attempt 0 + up to 3 reseeds
blk <- tbf_blocks(block_seams, block_H, nt)[[1]]   # seam 50, leads 1..15
catalog <- get_tbf_method_catalog()

dir.create("results", showWarnings = FALSE)

# ---- one guarded cell: fit -> forecast -> guard/reseed -> save ---------
# Self-contained so it can run inside a PSOCK worker. Reads globals
# (y, x, s, blk, block_H, catalog, iters, burn, update, max_attempts,
# phi.path, setting.label); returns a one-line status string.
run_one_cell <- function(setting, m, d) {
  outfile <- sprintf("results/blk1-%d-%d-%d.RData", setting, m, d)
  if (file.exists(outfile)) return(sprintf("s%d d%d m%d: exists, skip", setting, d, m))
  marginal_sd <- sd(y[, , , setting])
  truth <- list(lambda = 3, beta0 = 10)   # B uses lambda/beta0 only
  spec <- catalog[catalog$method_id == m, , drop = FALSE]
  y.train <- y[, blk$train_times, d, setting]
  x.train <- x[, blk$train_times, , drop = FALSE]
  x.block <- x[, blk$test_times, , drop = FALSE]
  y.val <- y[, blk$test_times, d, setting]
  tic <- proc.time()

  # Early stop for REPRODUCIBLE failures: if two consecutive attempts land on
  # the same lambda (within 15%) and both fail, the failure is data-level (a
  # flat long-memory z realization under-identifies lambda in a 50-day window)
  # -- reseeding cannot fix it. Never changes which fits PASS.
  chk <- NULL
  yhat <- NULL
  lam_prev <- NA_real_
  for (attempt in 0:(max_attempts - 1L)) {
    seed <- get_tbf_seed(setting, m, d) + 7919L * attempt
    set.seed(seed)
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
      x_block = x.block, s = s,
      ar2 = isTRUE(spec$ar2[1]),
      ar1 = isTRUE(spec$temporal[1]) && !isTRUE(spec$ar2[1]))
    chk <- check_fit_consistency(fit, yhat, truth, marginal_sd,
                                 data_mean = mean(y.train))
    chk <- cbind(setting = setting, method = m, dataset = d,
                 attempt = attempt, seed = seed,
                 pass = isTRUE(chk$B_truth) && isTRUE(chk$C_spread), chk)
    rm(fit)
    gc()
    if (chk$pass) break
    if (!is.na(lam_prev) &&
        abs(chk$lambda - lam_prev) < 0.15 * max(abs(lam_prev), 0.5)) break
    lam_prev <- chk$lambda
  }

  elapsed <- unname((proc.time() - tic)[3])
  res <- list(setting = setting, method_id = m, dataset = d,
              seam = blk$seam, leads = blk$leads, test_times = blk$test_times,
              yhat = yhat, y_val = y.val, chk = chk, elapsed_sec = elapsed)
  save(res, file = outfile)
  rm(res, yhat)
  gc()
  sprintf("s%d d%d m%d saved (%.0f s, %d attempt(s), lam %+.2f pass=%s)",
          setting, d, m, elapsed, chk$attempt + 1L, chk$lambda, chk$pass)
}

# ---- build the pending-cell grid (skip existing) ----------------------
grid <- expand.grid(setting = settings, m = method_ids, d = datasets,
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
grid <- grid[order(grid$setting, grid$d, grid$m), ]
grid$file <- sprintf("results/blk1-%d-%d-%d.RData", grid$setting, grid$m, grid$d)
pending <- grid[!file.exists(grid$file), , drop = FALSE]

cat(sprintf(
  "block1 positive control: settings=%s methods=%s datasets=%s\n",
  paste(settings, collapse = ","), paste(method_ids, collapse = ","),
  paste(range(datasets), collapse = "..")
))
cat(sprintf("block 1: train t=1..%d, forecast t=%d..%d (H=%d)\n",
            blk$seam, blk$seam + 1L, blk$seam + block_H, block_H))
cat(sprintf("cells: %d total, %d pending; workers = %d\n\n",
            nrow(grid), nrow(pending), workers))

if (nrow(pending) == 0L) {
  cat("nothing to do (all cells exist)\n")
} else if (workers > 1L) {
  wk <- min(workers, nrow(pending))
  cl <- parallel::makeCluster(wk, type = "PSOCK", outfile = "")
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  wd <- getwd()
  parallel::clusterExport(cl, "wd", envir = environment())
  parallel::clusterEvalQ(cl, {
    setwd(wd)
    source("../../simstudy/ar2_load.R", chdir = TRUE)
    source("../simstudy/time_block_helpers.R")
    source("./fit_assertions.R")
    load("../simstudy/simdata.RData", envir = .GlobalEnv)
    options(warn = 1)
    blk <- tbf_blocks(block_seams, block_H, nt)[[1]]
    catalog <- get_tbf_method_catalog()
    NULL
  })
  parallel::clusterExport(cl, c("iters", "burn", "update", "max_attempts",
                                "run_one_cell"), envir = environment())
  # Pass each cell's coords as the iterated argument -- do NOT reference the
  # master-side `pending` inside the worker (PSOCK resolves free names against
  # the worker's own global env, where it is unbound).
  cells <- lapply(seq_len(nrow(pending)), function(i) {
    c(setting = pending$setting[i], m = pending$m[i], d = pending$d[i])
  })
  res <- parallel::parLapply(cl, cells, function(cell) {
    run_one_cell(cell[["setting"]], cell[["m"]], cell[["d"]])
  })
  parallel::stopCluster(cl)
  cat(paste(unlist(res), collapse = "\n"), "\n")
} else {
  for (i in seq_len(nrow(pending))) {
    cat(run_one_cell(pending$setting[i], pending$m[i], pending$d[i]), "\n")
  }
}
cat("\nall requested fits finished\n")
