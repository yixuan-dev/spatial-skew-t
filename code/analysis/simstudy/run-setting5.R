#########################################################################
# Single-dataset replication script
# Data setting: 5 (Skew-t, K = 5, lambda = 3)
# Run all 6 analysis methods for one dataset index.
#########################################################################

# source("./package_load.R", chdir = TRUE)
source("./ar2_load.R", chdir = TRUE)

# -----------------------------
# User controls
# -----------------------------
setting <- 5
iters <- 20000
burn <- 10000
update <- 1000
thin <- 1

# Parse command-line arguments
# Usage: Rscript run-setting5.R [datasets] [workers] [ms_threads] [methods]
# datasets/methods examples: "1:5", "20:28", "c(1,2,6)", "(1,2,6)"
parse_index_expr <- function(expr_str, arg_name) {
  expr <- trimws(expr_str)
  if (is.na(expr) || expr == "") {
    stop(sprintf("%s specification cannot be empty", arg_name))
  }

  # Allow shorthand like (1,2,6) by converting to c(1,2,6)
  if (grepl("^\\(.*\\)$", expr)) {
    expr <- paste0("c", expr)
  }

  values_raw <- tryCatch(
    eval(parse(text = expr), envir = baseenv()),
    error = function(e) {
      stop(sprintf(
        "Invalid %s expression '%s'. Use forms like 1:5, c(1,2,6), (1,2,6)",
        arg_name, expr_str
      ))
    }
  )

  if (!is.numeric(values_raw) || length(values_raw) == 0) {
    stop(sprintf("%s must evaluate to a non-empty numeric vector", arg_name))
  }

  values_int <- as.integer(values_raw)
  if (any(is.na(values_int)) || any(values_raw != values_int)) {
    stop(sprintf("%s must contain integers only", arg_name))
  }

  sort(unique(values_int))
}

parse_dataset_spec <- function(dataset_str) {
  dataset_ids <- parse_index_expr(dataset_str, "datasets")
  if (any(dataset_ids < 1 | dataset_ids > 50)) {
    stop("datasets must be integers in 1..50")
  }
  dataset_ids
}

parse_methods_spec <- function(methods_str) {
  method_ids <- parse_index_expr(methods_str, "methods")
  if (any(method_ids < 1 | method_ids > 6)) {
    stop("methods must be integers in 1..6")
  }
  list(mcmc = method_ids[method_ids <= 5], maxstable = 6 %in% method_ids)
}

args <- commandArgs(trailingOnly = TRUE)
datasets_spec <- if (length(args) >= 1) args[1] else "1"
workers <- if (length(args) >= 2) as.integer(args[2]) else 1
ms_threads <- if (length(args) >= 3) as.integer(args[3]) else 2
methods_spec <- if (length(args) >= 4) args[4] else "1:6"

if (is.na(workers) || workers < 1) {
  stop("workers must be a positive integer")
}
if (is.na(ms_threads) || ms_threads < 1) {
  stop("ms_threads must be a positive integer")
}

dataset_ids <- parse_dataset_spec(datasets_spec)
methods_to_run <- parse_methods_spec(methods_spec)

# keep output path consistent with results*.R scripts
if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}

run_method_mcmc <- function(method_id, dataset_id) {
  obs <- c(rep(TRUE, 100), rep(FALSE, 44))
  y.d <- y[, , dataset_id, setting]
  y.o <- y.d[obs, ]
  x.o <- x[obs, , ]
  s.o <- s[obs, ]
  x.p <- x[!obs, , ]
  s.p <- s[!obs, ]

  # Seed rule (methods 1-5): method * 1000 + setting * 100 + dataset
  set.seed(method_id * 1000 + setting * 100 + dataset_id)
  outputfile <- sprintf("results/%d-%d-%d.RData", setting, method_id, dataset_id)

  cat(sprintf("[Dataset %d][Method %d] start\n", dataset_id, method_id))

  if (method_id == 1) {
    cat("[Method 1] Gaussian start\n")
    fit.1 <- mcmc(
      y = y.o, x = x.o, s = s.o, s.pred = s.p, x.pred = x.p,
      method = "gaussian", skew = FALSE, thresh.all = 0,
      thresh.quant = TRUE, nknots = 1, iterplot = FALSE, iters = iters,
      burn = burn, update = update, min.s = c(0, 0), max.s = c(10, 10),
      temporalw = FALSE, temporaltau = FALSE, temporalz = FALSE,
      rho.upper = 15, nu.upper = 10
    )
  } else if (method_id == 2) {
    cat("[Method 2] Skew-t, K=1 start\n")
    fit.1 <- mcmc(
      y = y.o, x = x.o, s = s.o, s.pred = s.p, x.pred = x.p,
      method = "t", skew = TRUE, thresh.all = 0,
      thresh.quant = TRUE, nknots = 1, iterplot = FALSE, iters = iters,
      burn = burn, update = update, min.s = c(0, 0), max.s = c(10, 10),
      temporalw = FALSE, temporaltau = FALSE, temporalz = FALSE,
      rho.upper = 15, nu.upper = 10
    )
  } else if (method_id == 3) {
    cat("[Method 3] t, K=1, T=q(0.80) start\n")
    fit.1 <- tryCatch(
      mcmc(
        y = y.o, x = x.o, s = s.o, s.pred = s.p, x.pred = x.p,
        method = "t", skew = FALSE, thresh.all = 0.80,
        thresh.quant = TRUE, nknots = 1, iterplot = FALSE, iters = iters,
        burn = burn, update = update, min.s = c(0, 0), max.s = c(10, 10),
        temporalw = FALSE, temporaltau = FALSE, temporalz = FALSE,
        rho.upper = 15, nu.upper = 10
      ),
      error = function(e) {
        mcmc(
          y = y.o, x = x.o, s = s.o, s.pred = s.p, x.pred = x.p,
          method = "t", skew = FALSE, thresh.all = 0.80,
          thresh.quant = TRUE, nknots = 1, iterplot = FALSE, iters = iters,
          burn = burn, update = update, min.s = c(0, 0), max.s = c(10, 10),
          temporalw = FALSE, temporaltau = FALSE, temporalz = FALSE,
          rho.upper = 15, nu.upper = 10, cov.model = "exponential"
        )
      }
    )
  } else if (method_id == 4) {
    cat("[Method 4] Skew-t, K=5 start\n")
    fit.1 <- mcmc(
      y = y.o, x = x.o, s = s.o, s.pred = s.p, x.pred = x.p,
      method = "t", skew = TRUE, thresh.all = 0,
      thresh.quant = TRUE, nknots = 5, iterplot = FALSE, iters = iters,
      burn = burn, update = update, min.s = c(0, 0), max.s = c(10, 10),
      temporalw = FALSE, temporaltau = FALSE, temporalz = FALSE,
      rho.upper = 15, nu.upper = 10
    )
  } else if (method_id == 5) {
    cat("[Method 5] t, K=5, T=q(0.80) start\n")
    fit.1 <- tryCatch(
      mcmc(
        y = y.o, x = x.o, s = s.o, s.pred = s.p, x.pred = x.p,
        method = "t", skew = FALSE, thresh.all = 0.80,
        thresh.quant = TRUE, nknots = 5, iterplot = FALSE, iters = iters,
        burn = burn, update = update, min.s = c(0, 0), max.s = c(10, 10),
        temporalw = FALSE, temporaltau = FALSE, temporalz = FALSE,
        rho.upper = 15, nu.upper = 10
      ),
      error = function(e) {
        mcmc(
          y = y.o, x = x.o, s = s.o, s.pred = s.p, x.pred = x.p,
          method = "t", skew = FALSE, thresh.all = 0.80,
          thresh.quant = TRUE, nknots = 5, iterplot = FALSE, iters = iters,
          burn = burn, update = update, min.s = c(0, 0), max.s = c(10, 10),
          temporalw = FALSE, temporaltau = FALSE, temporalz = FALSE,
          rho.upper = 15, nu.upper = 10, cov.model = "exponential"
        )
      }
    )
  } else {
    stop("Unsupported method_id")
  }

  save(fit.1, file = outputfile)
  rm(fit.1)
  gc()
  cat(sprintf("[Dataset %d][Method %d] done -> %s\n", dataset_id, method_id, outputfile))
}

run_method_maxstable <- function(dataset_id) {
  obs <- c(rep(TRUE, 100), rep(FALSE, 44))
  y.d <- y[, , dataset_id, setting]
  y.o <- y.d[obs, ]
  x.o <- x[obs, , ]
  s.o <- s[obs, ]
  x.p <- x[!obs, , ]
  s.p <- s[!obs, ]

  method_id <- 6
  # Seed rule (method 6): setting * 100 + dataset
  set.seed(setting * 100 + dataset_id)
  outputfile <- sprintf("results/%d-%d-%d.RData", setting, method_id, dataset_id)

  knots.x <- seq(1, 9, length = 12)
  knots <- expand.grid(knots.x, knots.x)

  y.ms <- t(y.o)
  thresh <- quantile(y.ms, probs = 0.80, na.rm = TRUE)

  cat(sprintf("[Dataset %d][Method 6] Max-stable, T=q(0.80) start\n", dataset_id))
  fit.1 <- maxstable(
    y = y.ms, x = x.o, s = s.o, sp = s.p, xp = x.p, thresh = thresh,
    knots = knots, iters = iters, burn = burn, update = update,
    threads = ms_threads, thin = thin
  )

  save(fit.1, file = outputfile)
  rm(fit.1)
  gc()
  cat(sprintf("[Dataset %d][Method 6] done -> %s\n", dataset_id, outputfile))
}

cat(sprintf(
  "Run setting=%d, datasets=%s, iters=%d, burn=%d, workers=%d, ms_threads=%d\n",
  setting, paste(dataset_ids, collapse = ","), iters, burn, workers, ms_threads
))
selected_methods <- c(methods_to_run$mcmc, if (methods_to_run$maxstable) 6)
cat(sprintf("Methods to run: %s\n", paste(selected_methods, collapse = ",")))

mcmc_methods <- methods_to_run$mcmc
run_maxstable <- methods_to_run$maxstable
mcmc_tasks <- if (length(mcmc_methods) > 0) {
  expand.grid(
    dataset = dataset_ids,
    method = mcmc_methods,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
} else {
  data.frame(dataset = integer(0), method = integer(0))
}
workers_use <- min(workers, nrow(mcmc_tasks))

if (workers_use > 1) {
  cat(sprintf("Running methods 1-5 in parallel with %d workers across %d tasks...\n", workers_use, nrow(mcmc_tasks)))
  cl <- parallel::makeCluster(workers_use, type = "PSOCK")

  wd <- getwd()
  parallel::clusterExport(cl, varlist = c("wd"), envir = environment())
  parallel::clusterEvalQ(cl, {
    setwd(wd)
    # source("./package_load.R", chdir = TRUE)
    source("./ar2_load.R", chdir = TRUE)
    NULL
  })
  parallel::clusterExport(
    cl,
    varlist = c("setting", "iters", "burn", "update", "thin", "mcmc_tasks", "run_method_mcmc"),
    envir = environment()
  )

  elapsed <- parallel::parLapply(cl, seq_len(nrow(mcmc_tasks)), function(i) {
    d <- mcmc_tasks$dataset[i]
    m <- mcmc_tasks$method[i]
    tic <- proc.time()
    run_method_mcmc(m, d)
    unname((proc.time() - tic)[3])
  })
  for (i in seq_len(nrow(mcmc_tasks))) {
    cat(sprintf(
      "elapsed (dataset %d, method %d): %.2f sec\n",
      mcmc_tasks$dataset[i], mcmc_tasks$method[i], elapsed[[i]]
    ))
  }

  parallel::stopCluster(cl)
} else {
  for (i in seq_len(nrow(mcmc_tasks))) {
    d <- mcmc_tasks$dataset[i]
    m <- mcmc_tasks$method[i]
    tic <- proc.time()
    run_method_mcmc(m, d)
    toc <- proc.time()
    cat(sprintf("elapsed (dataset %d, method %d): %.2f sec\n", d, m, (toc - tic)[3]))
  }
}

if (run_maxstable) {
  for (d in dataset_ids) {
    tic <- proc.time()
    run_method_maxstable(d)
    toc <- proc.time()
    cat(sprintf("elapsed (dataset %d, method 6): %.2f sec\n", d, (toc - tic)[3]))
  }
}

cat("All specified methods and datasets finished.\n")
