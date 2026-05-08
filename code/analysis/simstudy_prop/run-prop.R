#########################################################################
# Prop simulation launcher
# Keeps the simstudy run-settings.R CLI shape:
#   Rscript run-prop.R [--data=<path>|--data <path>] \
#                      [--setting=<id>|--setting <id>] \
#                      [datasets] [workers] [ms_threads] [methods] [mrts_k]
# In this launcher, the final [mrts_k] slot is reused as the prop basis rank.
#
# --data defaults to ./simdata.RData (with fallback to ../simstudy/simdata.RData).
# When --data is supplied with a non-default file (e.g. simdata_def.RData),
# fits_dir is derived from the basename: simdata_def.RData -> fits_def/
# (unless SIMSTUDY_PROP_FITS_DIR is set, which always wins).
#########################################################################

script_dir_override <- trimws(Sys.getenv("SIMSTUDY_PROP_SCRIPT_DIR", unset = ""))
if (nzchar(script_dir_override) && dir.exists(script_dir_override)) {
  setwd(normalizePath(script_dir_override, winslash = "/", mustWork = TRUE))
} else {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_arg) > 0) {
    script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = FALSE)
    script_dir <- dirname(script_path)
    if (dir.exists(script_dir)) {
      setwd(script_dir)
    }
  }
}

source("./prop_load.R", chdir = TRUE)
source("./prop_simstudy_helpers.R", chdir = TRUE)

default_setting <- 5L
iters <- as.integer(Sys.getenv("SIMSTUDY_PROP_ITERS", unset = "20000"))
burn <- as.integer(Sys.getenv("SIMSTUDY_PROP_BURN", unset = "10000"))
update <- as.integer(Sys.getenv("SIMSTUDY_PROP_UPDATE", unset = "1000"))
thin <- as.integer(Sys.getenv("SIMSTUDY_PROP_THIN", unset = "1"))
prop_cov_update_every <- as.integer(Sys.getenv("SIMSTUDY_PROP_COV_UPDATE_EVERY", unset = "1"))

args_raw <- commandArgs(trailingOnly = TRUE)
parsed_args <- extract_leading_flags(args_raw, c("data", "setting"))
args <- parsed_args$args
data_flag_value <- parsed_args$values$data
setting_flag_value <- parsed_args$values$setting

# Default = whatever prop_load.R already loaded. If --data was passed,
# resolve it (with fallback to ../simstudy/<basename>) and reload over
# the default-loaded objects so y/x/s/etc. reflect the requested dataset.
default_data_path <- resolve_simstudy_prop_data()
data_path <- resolve_simstudy_prop_data_path(data_flag_value)

if (!identical(data_path, default_data_path)) {
  load(data_path, envir = .GlobalEnv)
  cat(sprintf("Loaded dataset from %s\n", data_path))
}

datasets_spec <- if (length(args) >= 1L) args[1] else "1"
workers <- if (length(args) >= 2L) as.integer(args[2]) else 1L
ms_threads <- if (length(args) >= 3L) as.integer(args[3]) else 1L
methods_spec <- if (length(args) >= 4L) args[4] else "1"
prop_k_spec <- if (length(args) >= 5L) args[5] else ""

setting <- if (!is.null(setting_flag_value)) {
  parse_setting_spec(setting_flag_value, y)
} else {
  parse_setting_spec(as.character(default_setting), y)
}

if (is.na(workers) || workers < 1L) {
  stop("workers must be a positive integer", call. = FALSE)
}
if (is.na(ms_threads) || ms_threads < 1L) {
  stop("ms_threads must be a positive integer", call. = FALSE)
}
if (is.na(iters) || is.na(burn) || is.na(update) || is.na(thin) ||
    iters < 2L || burn < 0L || burn >= iters || update < 1L || thin < 1L) {
  stop("Invalid MCMC control values. Check SIMSTUDY_PROP_{ITERS,BURN,UPDATE,THIN}.", call. = FALSE)
}
if (is.na(prop_cov_update_every) || prop_cov_update_every < 1L) {
  stop("SIMSTUDY_PROP_COV_UPDATE_EVERY must be a positive integer.", call. = FALSE)
}

dataset_ids <- parse_dataset_spec(datasets_spec)
method_ids <- parse_prop_methods_spec(methods_spec)
prop_k_values <- parse_prop_k_spec(prop_k_spec)

method_plan <- build_prop_run_plan(method_ids = method_ids, prop_k_values = prop_k_values)
unsupported <- method_plan[!method_plan$prop_supported, , drop = FALSE]
if (nrow(unsupported) > 0L) {
  stop(
    sprintf(
      "The prop backend does not support the requested method keys: %s",
      paste(unsupported$method_key, collapse = ", ")
    ),
    call. = FALSE
  )
}

fits_dir_env <- trimws(Sys.getenv("SIMSTUDY_PROP_FITS_DIR", unset = ""))
if (nzchar(fits_dir_env)) {
  fits_dir <- fits_dir_env
} else {
  fits_dir <- derive_prop_results_dir(data_path, default_dir = "fits")
}
if (!dir.exists(fits_dir)) {
  dir.create(fits_dir, recursive = TRUE, showWarnings = FALSE)
}

write.csv(
  method_plan,
  file.path(fits_dir, sprintf("run-plan-setting-%d.csv", setting)),
  row.names = FALSE
)

cat(sprintf(
  "Prop run: data=%s fits_dir=%s setting=%d datasets=%s iters=%d burn=%d workers=%d ms_threads=%d cov_update_every=%d\n",
  data_path,
  fits_dir,
  setting,
  paste(dataset_ids, collapse = ","),
  iters,
  burn,
  workers,
  ms_threads,
  prop_cov_update_every
))
cat(sprintf("Prop method keys: %s\n", paste(method_plan$method_key, collapse = ", ")))

run_method_prop_mcmc <- function(plan_id, dataset_id) {
  spec <- method_plan[method_plan$plan_id == plan_id, , drop = FALSE]
  if (nrow(spec) != 1L) {
    stop(sprintf("Expected exactly one plan row for plan_id=%s", plan_id), call. = FALSE)
  }

  obs <- c(rep(TRUE, nrow(y) - ntest), rep(FALSE, ntest))
  y.d <- y[, , dataset_id, setting]
  y.o <- y.d[obs, , drop = FALSE]
  x.o <- x[obs, , , drop = FALSE]
  s.o <- s[obs, , drop = FALSE]
  x.p <- x[!obs, , , drop = FALSE]
  s.p <- s[!obs, , drop = FALSE]

  set.seed(get_prop_seed(setting, spec$method_id[1], dataset_id, spec$prop_k[1]))
  outputfile <- build_prop_result_file(
    results_dir = fits_dir,
    setting_id = setting,
    method_id = spec$method_id[1],
    dataset_id = dataset_id,
    prop_k = spec$prop_k[1]
  )

  cat(sprintf("[Dataset %d][Method %s] start\n", dataset_id, spec$method_key[1]))
  started_at <- Sys.time()
  tic <- proc.time()
  fit.1 <- mcmc(
    y = y.o,
    x = x.o,
    s = s.o,
    s.pred = s.p,
    x.pred = x.p,
    method = spec$method[1],
    skew = isTRUE(spec$skew[1]),
    thresh.all = spec$thresh_all[1],
    thresh.quant = isTRUE(spec$thresh_quant[1]),
    nknots = spec$nknots[1],
    iterplot = FALSE,
    iters = iters,
    burn = burn,
    update = update,
    thin = thin,
    min.s = c(0, 0),
    max.s = c(10, 10),
    temporalw = FALSE,
    temporaltau = FALSE,
    temporalz = FALSE,
    prop_k = spec$prop_k[1],
    prop_cov_update_every = prop_cov_update_every
  )
  elapsed_sec <- unname((proc.time() - tic)[3])
  analysis_spec <- spec
  prop_meta <- fit.1$prop
  runtime_info <- build_simstudy_runtime_info(
    started_at = started_at,
    elapsed_sec = elapsed_sec,
    outputfile = outputfile,
    setting_id = setting,
    dataset_id = dataset_id,
    runner = "simstudy_prop/run-prop.R::run_method_prop_mcmc",
    method_id = spec$method_id[1],
    method_key = spec$method_key[1],
    prop_k = spec$prop_k[1],
    backend = "prop-simstudy",
    control = list(
      iters = as.integer(iters),
      burn = as.integer(burn),
      update = as.integer(update),
      thin = as.integer(thin),
      cov_update_every = as.integer(prop_cov_update_every)
    )
  )

  save(fit.1, analysis_spec, prop_meta, runtime_info, file = outputfile)
  rm(fit.1, analysis_spec, prop_meta, runtime_info)
  gc()
  cat(sprintf("[Dataset %d][Method %s] done -> %s\n", dataset_id, spec$method_key[1], outputfile))
}

mcmc_tasks <- expand.grid(
  dataset = dataset_ids,
  plan_id = method_plan$plan_id,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

workers_use <- min(workers, nrow(mcmc_tasks))
if (workers_use > 1L) {
  cat(sprintf(
    "Running prop tasks in parallel with %d workers across %d tasks...\n",
    workers_use,
    nrow(mcmc_tasks)
  ))
  cl <- parallel::makeCluster(workers_use, type = "PSOCK")
  wd <- getwd()
  parallel::clusterExport(cl, varlist = c("wd"), envir = environment())
  parallel::clusterEvalQ(cl, {
    setwd(wd)
    source("./prop_load.R", chdir = TRUE)
    source("./prop_simstudy_helpers.R", chdir = TRUE)
    NULL
  })
  # prop_load.R doesn't rm() globals, but it does load default simdata
  # into the worker. Re-export the resolved data_path / default_data_path
  # and reload the alternate dataset on each worker if needed.
  parallel::clusterExport(
    cl,
    varlist = c("data_path", "default_data_path"),
    envir = environment()
  )
  parallel::clusterEvalQ(cl, {
    if (!identical(data_path, default_data_path)) {
      load(data_path, envir = .GlobalEnv)
    }
    NULL
  })
  parallel::clusterExport(
    cl,
    varlist = c(
      "setting", "iters", "burn", "update", "thin", "fits_dir",
      "method_plan", "mcmc_tasks", "prop_cov_update_every", "run_method_prop_mcmc"
    ),
    envir = environment()
  )

  elapsed <- parallel::parLapply(cl, seq_len(nrow(mcmc_tasks)), function(i) {
    d <- mcmc_tasks$dataset[i]
    pid <- mcmc_tasks$plan_id[i]
    tic <- proc.time()
    run_method_prop_mcmc(pid, d)
    unname((proc.time() - tic)[3])
  })
  for (i in seq_len(nrow(mcmc_tasks))) {
    method_key <- method_plan$method_key[match(mcmc_tasks$plan_id[i], method_plan$plan_id)]
    cat(sprintf(
      "elapsed (dataset %d, method %s): %.2f sec\n",
      mcmc_tasks$dataset[i],
      method_key,
      elapsed[[i]]
    ))
  }
  parallel::stopCluster(cl)
} else {
  for (i in seq_len(nrow(mcmc_tasks))) {
    d <- mcmc_tasks$dataset[i]
    pid <- mcmc_tasks$plan_id[i]
    method_key <- method_plan$method_key[match(pid, method_plan$plan_id)]
    tic <- proc.time()
    run_method_prop_mcmc(pid, d)
    toc <- proc.time()
    cat(sprintf("elapsed (dataset %d, method %s): %.2f sec\n", d, method_key, (toc - tic)[3]))
  }
}

cat("Prop tasks finished.\n")
