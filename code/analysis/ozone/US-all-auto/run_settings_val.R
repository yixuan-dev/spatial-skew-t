# run_settings_val.R
#
# Fits each candidate setting ONCE on train_sites and generates posterior
# predictive draws at val_sites.  No cross-validation — one MCMC run per setting.
#
# Output per setting: fits/val-<N>.RData containing:
#   fit          - MCMC result; fit$yp = draws x n_val x ntime
#   runtime_info - metadata (timing, model spec, run mode)
#
# How to run (PowerShell):
#   cd D:\Github\spatial-skew-t\code\analysis\ozone\US-all-auto
#   $env:US_ALL_AUTOSELECT_SETTINGS = '111,112'
#   Rscript run_settings_val.R
#
# Prerequisites:
#   Rscript us-all-setup-auto.R   (produces us-all-setup-auto.RData)
#
# Environment variables:
#   US_ALL_AUTOSELECT_SETTINGS   Range tokens e.g. "111,112,204:206"  [required]
#   US_ALL_VAL_RUN_MODE          "dev" or "prod"                       [default "prod"]
#   US_ALL_MCMC_BACKEND          "legacy" or "ar2"                     [default "ar2"]
#   US_ALL_VAL_RESULTS_DIR       Override fits/ output directory       [default: auto]

rm(list = ls())

# ---- Locate this script and load shared helpers ----

.this_script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", arg[1]), winslash = "/", mustWork = FALSE))
  } else {
    normalizePath(".", winslash = "/", mustWork = FALSE)
  }
})
setwd(.this_script_dir)

.us_all_auto_root <- normalizePath(.this_script_dir, winslash = "/", mustWork = FALSE)
.us_all_root <- normalizePath(file.path(.us_all_auto_root, "../US-all"),
  winslash = "/", mustWork = FALSE
)

source(file.path(.us_all_auto_root, "utils.R"))

# ---- Parse environment variables ----

RUN_MODE <- tolower(trimws(Sys.getenv("US_ALL_VAL_RUN_MODE", unset = "prod")))
if (!RUN_MODE %in% c("dev", "prod")) {
  stop("US_ALL_VAL_RUN_MODE must be 'dev' or 'prod'.", call. = FALSE)
}

mcmc_ctrl <- switch(RUN_MODE,
  dev  = list(iters = 2000, burn = 1000, update = 200),
  prod = list(iters = 30000, burn = 25000, update = 500)
)

backend <- tolower(trimws(Sys.getenv("US_ALL_MCMC_BACKEND", unset = "ar2")))
if (!backend %in% c("legacy", "ar2")) {
  stop("US_ALL_MCMC_BACKEND must be 'legacy' or 'ar2'.", call. = FALSE)
}

requested_ids <- expand_range_tokens(
  trimws(Sys.getenv("US_ALL_AUTOSELECT_SETTINGS", unset = ""))
)
if (length(requested_ids) == 0L) {
  stop("US_ALL_AUTOSELECT_SETTINGS must be set. Example: '111,112,204:206'", call. = FALSE)
}

fits_dir_override <- trimws(Sys.getenv("US_ALL_VAL_RESULTS_DIR", unset = ""))
fits_dir <- if (nzchar(fits_dir_override)) {
  normalizePath(fits_dir_override, winslash = "/", mustWork = FALSE)
} else {
  file.path(.us_all_auto_root, "fits")
}
dir.create(fits_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load packages and MCMC functions ----
# package_load.R reads `backend` and sources only the matching set of functions,
# avoiding name conflicts between ar2/ and legacy model files.

source(file.path(.us_all_auto_root, "package_load.R"))

if (!exists("mcmc", envir = .GlobalEnv, inherits = TRUE))
  stop(sprintf("MCMC function not available for backend '%s'.", backend), call. = FALSE)
run_mcmc <- get("mcmc", envir = .GlobalEnv, inherits = TRUE)
run_maxstable <- if (exists("maxstable", envir = .GlobalEnv, inherits = TRUE)) {
  get("maxstable", envir = .GlobalEnv, inherits = TRUE)
} else {
  NULL
}

# ---- Load setup (train/val partition + dataset) ----

setup_path <- file.path(.us_all_auto_root, "us-all-setup-auto.RData")
if (!file.exists(setup_path)) {
  stop("Setup file not found: ", setup_path, "\nRun:  Rscript us-all-setup-auto.R", call. = FALSE)
}

setup_env <- new.env(parent = emptyenv())
load(setup_path, envir = setup_env)

for (obj in c("Y", "X", "S", "split.lst", "beta.init", "tau.init")) {
  if (!exists(obj, envir = setup_env, inherits = FALSE)) {
    stop("us-all-setup-auto.RData is missing: ", obj, "\nRe-run us-all-setup-auto.R.", call. = FALSE)
  }
}

Y_data <- get("Y", envir = setup_env)
X_data <- get("X", envir = setup_env)
S_data <- get("S", envir = setup_env)
split.lst <- get("split.lst", envir = setup_env)
beta_init <- get("beta.init", envir = setup_env)
tau_init <- get("tau.init", envir = setup_env)

train_sites <- split.lst$train
val_sites <- split.lst$val

# ---- Load settings grid ----

settings_path <- file.path(.us_all_auto_root, "settings-auto.csv")
if (!file.exists(settings_path)) settings_path <- file.path(.us_all_root, "settings.csv")
if (!file.exists(settings_path)) stop("No settings file found. Run settings-auto.R first.", call. = FALSE)

settings <- read.csv(settings_path, stringsAsFactors = FALSE)
if (!"setting_num" %in% names(settings)) {
  settings$setting_num <- suppressWarnings(as.integer(settings$setting))
}

unknown <- setdiff(requested_ids, settings$setting_num[!is.na(settings$setting_num)])
if (length(unknown) > 0L) {
  message("Warning: settings not in settings file: ", paste(unknown, collapse = ", "))
}

cat("=== run_settings_val.R ===\n")
cat(sprintf("Backend   : %s\n", backend))
cat(sprintf("Run mode  : %s  (iters=%d, burn=%d)\n", RUN_MODE, mcmc_ctrl$iters, mcmc_ctrl$burn))
cat(sprintf("Train     : %d sites | Val: %d sites\n", length(train_sites), length(val_sites)))
cat(sprintf("Settings  : %s\n", paste(requested_ids, collapse = ", ")))
cat(sprintf("Fits dir  : %s\n\n", fits_dir))

# ---- Run one setting ----

run_one_setting <- function(setting_id) {
  row <- settings[!is.na(settings$setting_num) & settings$setting_num == setting_id, ,
    drop = FALSE
  ]
  if (nrow(row) == 0L) {
    message("Setting ", setting_id, " not found — skipped.")
    return(invisible(NULL))
  }
  row <- row[1L, ]

  model_spec <- resolve_model(row$method)
  threshold <- suppressWarnings(as.numeric(row$thresh))
  nknots <- suppressWarnings(as.integer(row$knots))
  use_cmaq <- truthy(row$CMAQ)
  temporal <- truthy(row$TS)
  mrts_k <- parse_optional_int(row$mrts)
  use_mrts <- !is.na(mrts_k) && mrts_k > 0

  if (is.na(threshold) || is.na(nknots)) {
    message("Setting ", setting_id, ": invalid thresh/knots — skipped.")
    return(invisible(NULL))
  }

  cat("\n==============================\n")
  cat(
    "Setting:", setting_id, "| Method:", row$method, "| K:", nknots,
    "| Threshold:", threshold, "\n"
  )
  cat(
    "CMAQ:", if (use_cmaq) "yes" else "no",
    "| TS:", if (temporal) "yes" else "no",
    "| MRTS:", if (use_mrts) paste0("k=", mrts_k) else "no", "\n"
  )

  if (model_spec$is_maxstable && use_mrts) {
    stop("MRTS not supported for max-stable settings.", call. = FALSE)
  }

  # Subset to train / val
  y_tr <- Y_data[train_sites, ]
  x_tr <- X_data[train_sites, , , drop = FALSE]
  x_vl <- X_data[val_sites, , , drop = FALSE]
  if (!use_cmaq) {
    x_tr <- x_tr[, , 1, drop = FALSE]
    x_vl <- x_vl[, , 1, drop = FALSE]
  }
  S_tr <- S_data[train_sites, ]
  S_vl <- S_data[val_sites, ]

  # Append MRTS columns
  if (use_mrts) {
    mrts_cov <- build_mrts_covariates(S_tr, S_vl, nt = dim(x_tr)[2], k = mrts_k)
    p_base <- dim(x_tr)[3]
    p_mrts <- dim(mrts_cov$train)[3]
    x_tr_ext <- array(NA_real_, dim = c(dim(x_tr)[1], dim(x_tr)[2], p_base + p_mrts))
    x_vl_ext <- array(NA_real_, dim = c(dim(x_vl)[1], dim(x_vl)[2], p_base + p_mrts))
    x_tr_ext[, , seq_len(p_base)] <- x_tr
    x_vl_ext[, , seq_len(p_base)] <- x_vl
    for (j in seq_len(p_mrts)) {
      x_tr_ext[, , p_base + j] <- mrts_cov$train[, , j]
      x_vl_ext[, , p_base + j] <- mrts_cov$pred[, , j]
    }
    x_tr <- x_tr_ext
    x_vl <- x_vl_ext
    cat(
      "MRTS:", mrts_cov$source,
      "| cols orig/kept/drop:", mrts_cov$original_cols, "/",
      mrts_cov$kept_cols, "/", mrts_cov$dropped_cols, "\n"
    )
  }

  outputfile <- file.path(fits_dir, sprintf("val-%d.RData", setting_id))
  set.seed(setting_id * 100 + 1L)
  started_at <- Sys.time()
  tic <- proc.time()

  if (model_spec$is_maxstable) {
    if (is.null(run_maxstable)) stop("maxstable() not loaded.", call. = FALSE)
    fit <- run_maxstable(
      y = t(y_tr), x = t(x_tr[, , 2]), s = S_tr, sp = S_vl,
      xp = t(x_vl[, , 2]), thresh = threshold, knots = S_tr,
      iters = mcmc_ctrl$iters, burn = mcmc_ctrl$burn,
      update = mcmc_ctrl$update, thin = 1
    )
  } else {
    call_args <- list(
      y = y_tr, s = S_tr, x = x_tr, x.pred = x_vl, s.pred = S_vl,
      method = model_spec$method, skew = model_spec$skew,
      keep.knots = FALSE, thresh.all = threshold, thresh.quant = FALSE,
      nknots = nknots, iters = mcmc_ctrl$iters, burn = mcmc_ctrl$burn,
      update = mcmc_ctrl$update, iterplot = FALSE,
      beta.init = beta_init, tau.init = pick_tau_init(model_spec$method, tau_init),
      gamma.init = 0.5, rho.init = 1, rho.upper = 5,
      nu.init = 0.5, nu.upper = 10, min.s = c(-2.25, -1.55), max.s = c(2.35, 1.30)
    )
    if (temporal) {
      call_args$temporaltau <- TRUE
      call_args$temporalw <- TRUE
      call_args$temporalz <- TRUE
    }
    if (backend == "ar2") {
      call_args$ar2_tau <- isTRUE(temporal)
      call_args$ar2_w <- isTRUE(temporal)
      call_args$ar2_z <- isTRUE(temporal)
    }
    fit <- do.call(run_mcmc, call_args)
  }

  elapsed <- (proc.time() - tic)[3]
  runtime_info <- list(
    schema_version = 1L, runner = "run_settings_val.R", output_file = outputfile,
    setting_id = as.integer(setting_id), method = as.character(model_spec$method),
    skew = isTRUE(model_spec$skew), nknots = as.integer(nknots),
    threshold = unname(as.numeric(threshold)), backend = backend, run_mode = RUN_MODE,
    use_cmaq = isTRUE(use_cmaq), temporal = isTRUE(temporal),
    use_mrts = isTRUE(use_mrts), mrts_k = if (is.na(mrts_k)) NA_integer_ else as.integer(mrts_k),
    n_train = length(train_sites), n_val = length(val_sites),
    started_at_utc = format(as.POSIXct(started_at, tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    elapsed_sec = unname(as.numeric(elapsed)),
    mcmc_control = list(
      iters = mcmc_ctrl$iters, burn = mcmc_ctrl$burn,
      update = mcmc_ctrl$update, thin = 1L
    )
  )

  save(fit, runtime_info, file = outputfile)
  cat("Saved:", outputfile, sprintf("(%.1f sec)\n", elapsed))
}

# ---- Main loop ----

for (sid in requested_ids) run_one_setting(sid)
cat("\nCompleted:", paste(requested_ids, collapse = ", "), "\n")
