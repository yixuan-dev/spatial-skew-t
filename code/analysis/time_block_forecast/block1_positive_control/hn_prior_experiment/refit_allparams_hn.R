#########################################################################
# refit_allparams_hn.R -- fill in the parameters check_fit_consistency()
# does not record: rho, nu, gamma, beta1, beta2 (and full posterior
# summaries for the rest), under lambda ~ HN(0, 20).
#
# The main HN run saved only yhat + chk (it rm()s the fit), so rho/nu/
# gamma/beta1/beta2 recovery cannot be read off disk. This refits a
# subset -- setting {4, 5} x method 2 x datasets 1..5 -- bit-for-bit with
# the same attempt-0 seeds, and stores posterior mean/sd/quantiles of
# every scalar parameter instead of the chains.
#
# Truth (setup.R): beta = (10, 0, 0), lambda = 3, rho = 1, nu = 0.5,
# gamma = 0.9, tau.alpha = 6, tau.beta = 16 (internal shape/rate 3, 8).
#
#   & $R refit_allparams_hn.R --datasets=1:5 --workers=6
#########################################################################

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1]),
    winslash = "/", mustWork = FALSE))
  if (dir.exists(script_dir)) setwd(script_dir)
}
source("../../../simstudy/ar2_load.R", chdir = TRUE)
source("../../simstudy/time_block_helpers.R")
load("../../simstudy/simdata.RData", envir = .GlobalEnv)
options(warn = 1)

args <- commandArgs(trailingOnly = TRUE)
getflag <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
datasets <- eval(parse(text = getflag("datasets", "1:5")))
settings <- eval(parse(text = paste0("c(", getflag("settings", "4,5"), ")")))
method_ids <- eval(parse(text = paste0("c(", getflag("methods", "2"), ")")))
workers <- as.integer(getflag("workers", "1"))

iters <- 20000
burn <- 10000
update <- 20001
blk <- tbf_blocks(block_seams, block_H, nt)[[1]]
catalog <- get_tbf_method_catalog()
dir.create("results_allparams", showWarnings = FALSE)

summarise_chain <- function(v) {
  v <- as.numeric(v)
  q <- quantile(v, c(0.025, 0.5, 0.975), na.rm = TRUE)
  c(mean = mean(v, na.rm = TRUE), sd = sd(v, na.rm = TRUE),
    q025 = q[[1]], q50 = q[[2]], q975 = q[[3]])
}

run_one_cell <- function(setting, m, d) {
  outfile <- sprintf("results_allparams/ap-%d-%d-%d.RData", setting, m, d)
  if (file.exists(outfile)) {
    return(sprintf("s%d d%d m%d: exists, skip", setting, d, m))
  }
  spec <- catalog[catalog$method_id == m, , drop = FALSE]
  y.train <- y[, blk$train_times, d, setting]
  x.train <- x[, blk$train_times, , drop = FALSE]
  tic <- proc.time()

  set.seed(get_tbf_seed(setting, m, d))
  fit <- mcmc(
    y = y.train, x = x.train, s = s,
    method = "t", skew = isTRUE(spec$skew[1]),
    thresh.all = 0, thresh.quant = TRUE, nknots = spec$nknots[1],
    iterplot = FALSE, iters = iters, burn = burn, update = update,
    min.s = c(0, 0), max.s = c(10, 10),
    temporalw = isTRUE(spec$temporal[1]),
    temporaltau = isTRUE(spec$temporal[1]),
    temporalz = isTRUE(spec$temporal[1]),
    ar2_w = isTRUE(spec$ar2[1]), ar2_tau = isTRUE(spec$ar2[1]),
    ar2_z = isTRUE(spec$ar2[1]),
    rho.upper = 15, nu.upper = 10,
    lambda.positive = TRUE
  )

  post <- list(
    beta0 = summarise_chain(fit$beta[, 1]),
    beta1 = summarise_chain(fit$beta[, 2]),
    beta2 = summarise_chain(fit$beta[, 3]),
    lambda = summarise_chain(fit$lambda),
    rho = summarise_chain(fit$rho),
    nu = summarise_chain(fit$nu),
    gamma = summarise_chain(fit$gamma),
    tau.alpha = summarise_chain(fit$tau.alpha / 2),
    tau.beta = summarise_chain(fit$tau.beta / 2)
  )
  add_phi <- function(p, nm) {
    if (is.null(p)) return(invisible())
    pm <- as.matrix(p)
    post[[paste0(nm, "1")]] <<- summarise_chain(pm[, 1])
    if (ncol(pm) >= 2) post[[paste0(nm, "2")]] <<- summarise_chain(pm[, 2])
  }
  add_phi(fit$phi.z, "phi.z")
  add_phi(fit$phi.tau, "phi.tau")
  add_phi(fit$phi.w, "phi.w")

  elapsed <- unname((proc.time() - tic)[3])
  res <- list(setting = setting, method_id = m, dataset = d,
              post = post, elapsed_sec = elapsed,
              prior = "lambda ~ HN(0, 20)")
  save(res, file = outfile)
  rm(fit, res)
  gc()
  sprintf("s%d d%d m%d done (%.0f s)", setting, d, m, elapsed)
}

grid <- expand.grid(setting = settings, m = method_ids, d = datasets,
                    KEEP.OUT.ATTRS = FALSE)
grid <- grid[order(grid$setting, grid$d, grid$m), ]
grid$file <- sprintf("results_allparams/ap-%d-%d-%d.RData",
                     grid$setting, grid$m, grid$d)
pending <- grid[!file.exists(grid$file), , drop = FALSE]
cat(sprintf("all-params refit: %d total, %d pending, workers = %d\n\n",
            nrow(grid), nrow(pending), workers))

if (nrow(pending) == 0L) {
  cat("nothing to do\n")
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
    load("../../simstudy/simdata.RData", envir = .GlobalEnv)
    options(warn = 1)
    blk <- tbf_blocks(block_seams, block_H, nt)[[1]]
    catalog <- get_tbf_method_catalog()
    NULL
  })
  parallel::clusterExport(cl, c("iters", "burn", "update", "run_one_cell",
                               "summarise_chain"), envir = environment())
  cells <- lapply(seq_len(nrow(pending)), function(i) {
    c(setting = pending$setting[i], m = pending$m[i], d = pending$d[i])
  })
  out <- parallel::parLapply(cl, cells, function(cell) {
    run_one_cell(cell[["setting"]], cell[["m"]], cell[["d"]])
  })
  parallel::stopCluster(cl)
  cat(paste(unlist(out), collapse = "\n"), "\n")
} else {
  for (i in seq_len(nrow(pending))) {
    cat(run_one_cell(pending$setting[i], pending$m[i], pending$d[i]), "\n")
  }
}
cat("\nall-params refit finished\n")
