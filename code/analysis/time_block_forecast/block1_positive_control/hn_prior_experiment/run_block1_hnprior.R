#########################################################################
# run_block1_hnprior.R -- HN-prior counterfactual of the block-1 guard.
#
# Question: replace the guard (assertion C + reseed) with the model-level
# fix lambda ~ HN(0, 20) (lambda.positive = TRUE): how many attempt-0
# fits are healthy?
#
# Design: IDENTICAL cells and seeds to run_block1_guarded.R -- settings
# {4 strong, 5 nearunit} x methods {1 iid, 2 AR(2)} x datasets 1..10,
# block 1 (train t=1..50, forecast t=51..65), seed =
# get_tbf_seed(setting, m, d), i.e. every cell's attempt-0 seed. ONE fit
# per cell: no guard, no reseed. Assertions A/A'/B/C are RECORDED as
# diagnostics only (C is not a gate here -- the prior replaces it).
# Head-to-head readout: the guarded study saw 13/40 attempt-0 failures;
# the HN prior removes the reflected branch (lambda < 0 impossible), so
# the open question is how many cells land healthy vs. positive-ridge
# inflated (cf. s5 d6 m2, lambda0 = +5.40) vs. under-identified (s5 d2).
#
#   $R = "C:\Program Files\R\R-4.5.1\bin\Rscript.exe"
#   & $R run_block1_hnprior.R --settings=4,5 --datasets=1:10 --methods=1:2 --workers=6
#########################################################################

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1]),
    winslash = "/", mustWork = FALSE
  ))
  if (dir.exists(script_dir)) setwd(script_dir)
}

# ar2_load.R does rm(list = ls()): source it FIRST, then the helpers it wiped.
source("../../../simstudy/ar2_load.R", chdir = TRUE)
source("../../simstudy/time_block_helpers.R")
source("../../simstudy/fit_diag_utils.R")
load("../../simstudy/simdata.RData", envir = .GlobalEnv)
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
blk <- tbf_blocks(block_seams, block_H, nt)[[1]]   # seam 50, leads 1..15
catalog <- get_tbf_method_catalog()

dir.create("results", showWarnings = FALSE)

# ---- one cell: single HN-prior fit -> forecast -> record assertions ----
run_one_cell <- function(setting, m, d) {
  outfile <- sprintf("results/hn-%d-%d-%d.RData", setting, m, d)
  if (file.exists(outfile)) return(sprintf("s%d d%d m%d: exists, skip", setting, d, m))
  marginal_sd <- sd(y[, , , setting])
  truth <- list(lambda = 3, beta0 = 10)
  spec <- catalog[catalog$method_id == m, , drop = FALSE]
  y.train <- y[, blk$train_times, d, setting]
  x.train <- x[, blk$train_times, , drop = FALSE]
  x.block <- x[, blk$test_times, , drop = FALSE]
  y.val <- y[, blk$test_times, d, setting]
  tic <- proc.time()

  seed <- get_tbf_seed(setting, m, d)          # attempt-0 seed, no reseed
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
    rho.upper = 15, nu.upper = 10,
    lambda.positive = TRUE                     # <-- the intervention
  )
  yhat <- forecast_block(
    fit = fit, seam = blk$seam, H = block_H,
    x_block = x.block, s = s,
    ar2 = isTRUE(spec$ar2[1]),
    ar1 = isTRUE(spec$temporal[1]) && !isTRUE(spec$ar2[1]))
  chk <- check_fit_consistency(fit, yhat, truth, marginal_sd,
                               data_mean = mean(y.train))
  chk <- cbind(setting = setting, method = m, dataset = d,
               attempt = 0L, seed = seed,
               pass = isTRUE(chk$B_truth) && isTRUE(chk$C_spread), chk)
  rm(fit)
  gc()

  elapsed <- unname((proc.time() - tic)[3])
  res <- list(setting = setting, method_id = m, dataset = d,
              seam = blk$seam, leads = blk$leads, test_times = blk$test_times,
              yhat = yhat, y_val = y.val, chk = chk, elapsed_sec = elapsed,
              prior = "lambda ~ HN(0, 20) via lambda.positive = TRUE")
  save(res, file = outfile)
  rm(res, yhat)
  gc()
  sprintf("s%d d%d m%d saved (%.0f s, lam %+.2f b0 %.1f sdzr %.2f B=%s C=%s)",
          setting, d, m, elapsed, chk$lambda, chk$beta0, chk$sdz_ratio,
          chk$B_truth, chk$C_spread)
}

# ---- build the pending-cell grid (skip existing) ----------------------
grid <- expand.grid(setting = settings, m = method_ids, d = datasets,
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
grid <- grid[order(grid$setting, grid$d, grid$m), ]
grid$file <- sprintf("results/hn-%d-%d-%d.RData", grid$setting, grid$m, grid$d)
pending <- grid[!file.exists(grid$file), , drop = FALSE]

cat(sprintf(
  "block1 HN-prior experiment: settings=%s methods=%s datasets=%s\n",
  paste(settings, collapse = ","), paste(method_ids, collapse = ","),
  paste(range(datasets), collapse = "..")
))
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
    source("../../../simstudy/ar2_load.R", chdir = TRUE)
    source("../../simstudy/time_block_helpers.R")
    source("../../simstudy/fit_diag_utils.R")
    load("../../simstudy/simdata.RData", envir = .GlobalEnv)
    options(warn = 1)
    blk <- tbf_blocks(block_seams, block_H, nt)[[1]]
    catalog <- get_tbf_method_catalog()
    NULL
  })
  parallel::clusterExport(cl, c("iters", "burn", "update", "run_one_cell"),
                          envir = environment())
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
