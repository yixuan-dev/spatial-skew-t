# run_settings_val.R
#
# Fits each candidate setting ONCE on train_sites and generates posterior
# predictive draws at val_sites.  No cross-validation — one MCMC run per setting.
#
# Output per setting: fits/<setting>.RData containing:
#   fit          - MCMC result object; fit$yp = draws x n_val x ntime
#   runtime_info - metadata list (timing, model spec, run mode)
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
#   US_ALL_AUTOSELECT_SETTINGS   Range tokens, e.g. "111,112,204:206"  [required]
#   US_ALL_VAL_RUN_MODE          "dev" or "prod"                        [default "prod"]
#   US_ALL_MCMC_BACKEND          "legacy" or "ar2"                      [default "legacy"]
#   US_ALL_VAL_RESULTS_DIR       Override fits/ output directory        [default: auto]

rm(list = ls())

library(compiler)
enableJIT(3)

# ---- Locate this script ----

.this_script_dir <- local({
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", script_arg[1]),
                          winslash = "/", mustWork = FALSE))
  } else {
    normalizePath(".", winslash = "/", mustWork = FALSE)
  }
})
setwd(.this_script_dir)

.us_all_auto_root <- normalizePath(.this_script_dir, winslash = "/", mustWork = FALSE)
.us_all_root      <- normalizePath(file.path(.us_all_auto_root, "../US-all"),
                                   winslash = "/", mustWork = FALSE)

# ---- Helpers ----

truthy <- function(x) tolower(trimws(as.character(x))) %in% c("1", "true", "yes", "y")

parse_optional_int <- function(x) {
  txt <- trimws(as.character(x[1]))
  if (!nzchar(txt)) return(NA_integer_)
  val <- suppressWarnings(as.numeric(txt))
  if (!is.finite(val)) return(NA_integer_)
  as.integer(round(val))
}

expand_range_tokens <- function(raw_str) {
  tokens <- unlist(strsplit(raw_str, ",", fixed = TRUE))
  out <- integer(0L)
  for (tok in tokens) {
    tok <- trimws(tok)
    if (!nzchar(tok)) next
    if (grepl("^[0-9]+:[0-9]+$", tok)) {
      parts <- as.integer(strsplit(tok, ":", fixed = TRUE)[[1]])
      if (length(parts) == 2L && !anyNA(parts)) {
        step <- if (parts[1L] <= parts[2L]) 1L else -1L
        out  <- c(out, seq.int(parts[1L], parts[2L], by = step))
      }
    } else {
      val <- suppressWarnings(as.integer(tok))
      if (!is.na(val)) out <- c(out, val)
    }
  }
  sort(unique(out))
}

# ---- Parse environment variables ----

RUN_MODE <- tolower(trimws(Sys.getenv("US_ALL_VAL_RUN_MODE", unset = "prod")))
if (!RUN_MODE %in% c("dev", "prod")) {
  stop("US_ALL_VAL_RUN_MODE must be 'dev' or 'prod'.", call. = FALSE)
}

mcmc_ctrl <- switch(RUN_MODE,
  dev  = list(iters = 2000,  burn = 1000,  update = 200),
  prod = list(iters = 30000, burn = 25000, update = 500)
)

backend <- tolower(trimws(Sys.getenv("US_ALL_MCMC_BACKEND", unset = "legacy")))
if (!backend %in% c("legacy", "ar2")) {
  stop("US_ALL_MCMC_BACKEND must be 'legacy' or 'ar2'.", call. = FALSE)
}

settings_raw <- trimws(Sys.getenv("US_ALL_AUTOSELECT_SETTINGS", unset = ""))
if (!nzchar(settings_raw)) {
  stop(
    "US_ALL_AUTOSELECT_SETTINGS must be set.\n",
    "Example:  $env:US_ALL_AUTOSELECT_SETTINGS = '111,112,204:206'",
    call. = FALSE
  )
}
requested_ids <- expand_range_tokens(settings_raw)
if (length(requested_ids) == 0L) {
  stop("No valid settings parsed from US_ALL_AUTOSELECT_SETTINGS='", settings_raw, "'",
       call. = FALSE)
}

fits_dir_override <- trimws(Sys.getenv("US_ALL_VAL_RESULTS_DIR", unset = ""))
fits_dir <- if (nzchar(fits_dir_override)) {
  normalizePath(fits_dir_override, winslash = "/", mustWork = FALSE)
} else {
  file.path(.us_all_auto_root, "fits")
}
dir.create(fits_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load backend (MCMC functions + shared dataset) ----

backend_env <- new.env(parent = globalenv())
if (backend == "ar2") {
  sys.source(file.path(.us_all_root, "ar2_load.R"),     envir = backend_env, chdir = TRUE)
} else {
  sys.source(file.path(.us_all_root, "package_load.R"), envir = backend_env, chdir = TRUE)
}

# ---- Load setup (train/val partition + dataset) ----

setup_path <- file.path(.us_all_auto_root, "us-all-setup-auto.RData")
if (!file.exists(setup_path)) {
  stop("Setup file not found: ", setup_path,
       "\nRun:  Rscript us-all-setup-auto.R", call. = FALSE)
}
setup_env <- new.env(parent = emptyenv())
load(setup_path, envir = setup_env)

required <- c("Y", "X", "S", "split.lst", "beta.init", "tau.init")
missing  <- required[!vapply(required, exists, logical(1), envir = setup_env, inherits = FALSE)]
if (length(missing) > 0L) {
  stop("us-all-setup-auto.RData is missing: ", paste(missing, collapse = ", "),
       "\nRe-run us-all-setup-auto.R.", call. = FALSE)
}

Y_data      <- get("Y",         envir = setup_env)
X_data      <- get("X",         envir = setup_env)
S_data      <- get("S",         envir = setup_env)
split.lst   <- get("split.lst", envir = setup_env)
beta_init   <- get("beta.init", envir = setup_env)
tau_init    <- get("tau.init",  envir = setup_env)

train_sites <- split.lst$train
val_sites   <- split.lst$val

cat("=== run_settings_val.R ===\n")
cat(sprintf("Backend   : %s\n", backend))
cat(sprintf("Run mode  : %s  (iters=%d, burn=%d)\n",
    RUN_MODE, mcmc_ctrl$iters, mcmc_ctrl$burn))
cat(sprintf("Train     : %d sites\n", length(train_sites)))
cat(sprintf("Val       : %d sites\n", length(val_sites)))
cat(sprintf("Settings  : %s\n", paste(requested_ids, collapse = ", ")))
cat(sprintf("Fits dir  : %s\n\n", fits_dir))

# ---- Resolve MCMC function ----

resolve_fn <- function(fname) {
  if (exists(fname, envir = backend_env, mode = "function", inherits = FALSE))
    return(get(fname, envir = backend_env, inherits = FALSE))
  if (exists(fname, mode = "function", inherits = TRUE))
    return(get(fname, inherits = TRUE))
  NULL
}

run_mcmc    <- if (backend == "ar2") resolve_fn("mcmc_ar2") else resolve_fn("mcmc")
run_maxstable <- resolve_fn("maxstable")

if (is.null(run_mcmc)) {
  stop(sprintf("MCMC function not available for backend '%s'.", backend), call. = FALSE)
}

# ---- Load settings grid ----

settings_path <- file.path(.us_all_auto_root, "settings-auto.csv")
if (!file.exists(settings_path))
  settings_path <- file.path(.us_all_root, "settings.csv")
if (!file.exists(settings_path))
  stop("No settings file found. Run settings-auto.R first.", call. = FALSE)

settings <- read.csv(settings_path, stringsAsFactors = FALSE)
if (!"setting_num" %in% names(settings))
  settings$setting_num <- suppressWarnings(as.integer(settings$setting))

unknown <- setdiff(requested_ids, settings$setting_num[!is.na(settings$setting_num)])
if (length(unknown) > 0L)
  message("Warning: settings not in settings file: ", paste(unknown, collapse = ", "))

# ---- MRTS helper (identical logic to us-all-run.R) ----

build_mrts_covariates <- function(S_train, S_pred, nt, k) {
  if (is.na(k) || k <= 0) stop("MRTS k must be a positive integer.", call. = FALSE)

  safe_mrts_call <- function(fn, label) {
    train_obj <- tryCatch(fn(S_train, k = k), error = function(e) NULL)
    pred_obj  <- tryCatch(fn(S_train, k = k, x = S_pred), error = function(e) NULL)
    if (!is.null(train_obj) && !is.null(pred_obj)) {
      return(list(train = as.matrix(train_obj), pred = as.matrix(pred_obj),
                  source = label))
    }
    all_obj <- tryCatch(fn(rbind(S_train, S_pred), k = k), error = function(e) NULL)
    if (is.null(all_obj)) return(NULL)
    all_mat <- as.matrix(all_obj)
    if (nrow(all_mat) != nrow(S_train) + nrow(S_pred)) return(NULL)
    n_tr <- nrow(S_train)
    list(train = all_mat[seq_len(n_tr), , drop = FALSE],
         pred  = all_mat[(n_tr + 1):nrow(all_mat), , drop = FALSE],
         source = paste0(label, " [combined fallback]"))
  }

  basis_obj <- NULL
  if (requireNamespace("autoFRK", quietly = TRUE) &&
      exists("mrts", where = asNamespace("autoFRK"), inherits = FALSE)) {
    basis_obj <- safe_mrts_call(autoFRK::mrts, "autoFRK::mrts")
  }
  if (is.null(basis_obj) && exists("mrts", mode = "function"))
    basis_obj <- safe_mrts_call(mrts, "mrts()")
  if (is.null(basis_obj))
    stop("Unable to build MRTS covariates (mrts function unavailable or failed).", call. = FALSE)

  col_sd <- apply(basis_obj$train, 2, sd)
  keep   <- which(is.finite(col_sd) & col_sd > 1e-10)
  if (length(keep) == 0L)
    stop("All MRTS columns are near-constant after construction.", call. = FALSE)

  basis_train <- basis_obj$train[, keep, drop = FALSE]
  basis_pred  <- basis_obj$pred[, keep, drop = FALSE]

  mrts_train <- array(0, dim = c(nrow(S_train), nt, ncol(basis_train)))
  mrts_pred  <- array(0, dim = c(nrow(S_pred),  nt, ncol(basis_pred)))
  for (t in seq_len(nt)) {
    mrts_train[, t, ] <- basis_train
    mrts_pred[, t, ]  <- basis_pred
  }

  list(train = mrts_train, pred = mrts_pred, source = basis_obj$source,
       original_cols = ncol(basis_obj$train),
       kept_cols     = ncol(basis_train),
       dropped_cols  = ncol(basis_obj$train) - ncol(basis_train))
}

# ---- Resolve model spec ----

resolve_model <- function(method_raw) {
  switch(tolower(trimws(method_raw)),
    "skew-t"     = list(method = "t",          skew = TRUE,  is_maxstable = FALSE),
    "t"          = list(method = "t",          skew = FALSE, is_maxstable = FALSE),
    "gaussian"   = list(method = "gaussian",   skew = FALSE, is_maxstable = FALSE),
    "max-stable" = list(method = "max-stable", skew = FALSE, is_maxstable = TRUE),
    stop(sprintf("Unsupported method: %s", method_raw), call. = FALSE)
  )
}

pick_tau_init <- function(method_name) {
  if (tolower(method_name) == "gaussian") tau_init else 0.05
}

# ---- Run one setting ----

run_one_setting <- function(setting_id) {
  row <- settings[!is.na(settings$setting_num) & settings$setting_num == setting_id,
                  , drop = FALSE]
  if (nrow(row) == 0L) {
    message("Setting ", setting_id, " not found in settings file — skipped.")
    return(invisible(NULL))
  }
  row <- row[1L, ]

  model_spec <- resolve_model(row$method)
  threshold  <- suppressWarnings(as.numeric(row$thresh))
  nknots     <- suppressWarnings(as.integer(row$knots))
  use_cmaq   <- truthy(row$CMAQ)
  temporal   <- truthy(row$TS)
  mrts_k     <- parse_optional_int(row$mrts)
  use_mrts   <- !is.na(mrts_k) && mrts_k > 0

  if (is.na(threshold) || is.na(nknots)) {
    message("Setting ", setting_id, ": invalid thresh/knots — skipped.")
    return(invisible(NULL))
  }

  cat("\n==============================\n")
  cat("Setting:", setting_id, "| Backend:", backend, "\n")
  cat("RUN_MODE:", RUN_MODE,
      "| iters:", mcmc_ctrl$iters, "| burn:", mcmc_ctrl$burn, "\n")
  cat("Method:", row$method, "| K:", nknots, "| Threshold:", threshold, "\n")
  cat(
    "CMAQ:", if (use_cmaq) "yes" else "no",
    "| TS:", if (temporal) "yes" else "no",
    "| MRTS:", if (use_mrts) paste0("k=", mrts_k) else "no",
    "\n"
  )

  if (model_spec$is_maxstable && use_mrts)
    stop("MRTS not supported for max-stable settings.", call. = FALSE)

  # ---- Subset data ----

  y_tr <- Y_data[train_sites, ]

  x_tr <- X_data[train_sites, , , drop = FALSE]
  x_vl <- X_data[val_sites,   , , drop = FALSE]

  if (!use_cmaq) {
    x_tr <- x_tr[, , 1, drop = FALSE]
    x_vl <- x_vl[, , 1, drop = FALSE]
  }

  S_tr <- S_data[train_sites, ]
  S_vl <- S_data[val_sites,   ]

  # ---- MRTS covariates ----

  if (use_mrts) {
    mrts_cov <- build_mrts_covariates(
      S_train = S_tr, S_pred = S_vl,
      nt = dim(x_tr)[2], k = mrts_k
    )
    p_base <- dim(x_tr)[3]
    p_mrts <- dim(mrts_cov$train)[3]

    x_tr_ext <- array(NA_real_, dim = c(dim(x_tr)[1], dim(x_tr)[2], p_base + p_mrts))
    x_vl_ext <- array(NA_real_, dim = c(dim(x_vl)[1], dim(x_vl)[2], p_base + p_mrts))

    x_tr_ext[, , seq_len(p_base)] <- x_tr
    x_vl_ext[, , seq_len(p_base)] <- x_vl
    for (j in seq_len(p_mrts)) {
      x_tr_ext[, , p_base + j] <- mrts_cov$train[, , j]
      x_vl_ext[, , p_base + j] <- mrts_cov$pred[,  , j]
    }

    x_tr <- x_tr_ext
    x_vl <- x_vl_ext

    cat("MRTS source:", mrts_cov$source,
        "| cols(orig/kept/drop):",
        mrts_cov$original_cols, "/", mrts_cov$kept_cols, "/", mrts_cov$dropped_cols,
        "| total covariates:", dim(x_tr)[3], "\n")
  }

  # ---- MCMC call ----

  outputfile <- file.path(fits_dir, sprintf("val-%d.RData", setting_id))
  set.seed(setting_id * 100 + 1L)

  started_at <- Sys.time()
  tic        <- proc.time()

  if (model_spec$is_maxstable) {
    if (is.null(run_maxstable))
      stop("maxstable() is not loaded in current backend.", call. = FALSE)
    fit <- run_maxstable(
      y      = t(y_tr),
      x      = t(x_tr[, , 2]),
      s      = S_tr,
      sp     = S_vl,
      xp     = t(x_vl[, , 2]),
      thresh = threshold,
      knots  = S_tr,
      iters  = mcmc_ctrl$iters,
      burn   = mcmc_ctrl$burn,
      update = mcmc_ctrl$update,
      thin   = 1
    )
  } else {
    call_args <- list(
      y          = y_tr,
      s          = S_tr,
      x          = x_tr,
      x.pred     = x_vl,
      s.pred     = S_vl,
      method     = model_spec$method,
      skew       = model_spec$skew,
      keep.knots = FALSE,
      thresh.all   = threshold,
      thresh.quant = FALSE,
      nknots     = nknots,
      iters      = mcmc_ctrl$iters,
      burn       = mcmc_ctrl$burn,
      update     = mcmc_ctrl$update,
      iterplot   = FALSE,
      beta.init  = beta_init,
      tau.init   = pick_tau_init(model_spec$method),
      gamma.init = 0.5,
      rho.init   = 1,
      rho.upper  = 5,
      nu.init    = 0.5,
      nu.upper   = 10,
      min.s      = c(-2.25, -1.55),
      max.s      = c(2.35,  1.30)
    )

    if (temporal) {
      call_args$temporaltau <- TRUE
      call_args$temporalw   <- TRUE
      call_args$temporalz   <- TRUE
    }

    if (backend == "ar2") {
      call_args$ar2_tau <- isTRUE(temporal)
      call_args$ar2_w   <- isTRUE(temporal)
      call_args$ar2_z   <- isTRUE(temporal)
    }

    fit <- do.call(run_mcmc, call_args)
  }

  elapsed <- (proc.time() - tic)[3]

  runtime_info <- list(
    schema_version = 1L,
    runner         = "run_settings_val.R",
    output_file    = outputfile,
    setting_id     = as.integer(setting_id),
    method         = as.character(model_spec$method),
    skew           = isTRUE(model_spec$skew),
    nknots         = as.integer(nknots),
    threshold      = unname(as.numeric(threshold)),
    backend        = backend,
    run_mode       = RUN_MODE,
    use_cmaq       = isTRUE(use_cmaq),
    temporal       = isTRUE(temporal),
    use_mrts       = isTRUE(use_mrts),
    mrts_k         = if (is.na(mrts_k)) NA_integer_ else as.integer(mrts_k),
    n_train        = length(train_sites),
    n_val          = length(val_sites),
    started_at_utc = format(as.POSIXct(started_at, tz = "UTC"),
                            "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    elapsed_sec    = unname(as.numeric(elapsed)),
    mcmc_control   = list(iters = mcmc_ctrl$iters, burn = mcmc_ctrl$burn,
                          update = mcmc_ctrl$update, thin = 1L)
  )

  save(fit, runtime_info, file = outputfile)
  cat("Saved:", outputfile, sprintf("(%.1f sec)\n", elapsed))
}

# ---- Main loop ----

for (sid in requested_ids) {
  run_one_setting(sid)
}

cat("\nAll requested settings completed:", paste(requested_ids, collapse = ", "), "\n")
