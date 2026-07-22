rm(list = ls())

library(compiler)
enableJIT(3)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = FALSE)
  script_dir <- dirname(script_path)
  if (dir.exists(script_dir)) {
    setwd(script_dir)
  }
}

truthy <- function(x) {
  tolower(trimws(as.character(x))) %in% c("1", "true", "yes", "y")
}

# Tiny runtime switch for MCMC lengths.
# Use "dev" for quick checks and "prod" for full experiments.
RUN_MODE <- tolower(Sys.getenv("US_ALL_RUN_MODE", unset = "prod"))
if (!RUN_MODE %in% c("dev", "prod")) {
  stop("RUN_MODE must be one of: dev, prod", call. = FALSE)
}

mcmc_ctrl <- switch(RUN_MODE,
  dev = list(iters = 2000, burn = 1000, update = 200),
  prod = list(iters = 30000, burn = 25000, update = 500)
)

# US_ALL_RESULTS_DIR is documented in US-ALL-RUN-GUIDE.md; without this
# override, a caller who sets it to protect results/ silently overwrites
# results/ anyway. Relative paths resolve against this script's directory.
results_dir <- Sys.getenv("US_ALL_RESULTS_DIR", unset = "results")
if (!grepl("^([A-Za-z]:|/|\\\\)", results_dir)) {
  results_dir <- file.path(getwd(), results_dir)
}
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  cat("Created results directory:", results_dir, "\n")
}

backend <- tolower(Sys.getenv("US_ALL_MCMC_BACKEND", unset = "legacy"))
if (!backend %in% c("legacy", "ar2")) {
  stop("US_ALL_MCMC_BACKEND must be one of: legacy, ar2", call. = FALSE)
}

backend_env <- new.env(parent = globalenv())
if (backend == "ar2") {
  sys.source("./ar2_load.R", envir = backend_env, chdir = TRUE)
} else {
  sys.source("./package_load.R", envir = backend_env, chdir = TRUE)
}

settings <- read.csv("settings.csv", stringsAsFactors = FALSE)
settings$setting_chr <- trimws(as.character(settings$setting))

required_objects <- c("Y", "X", "S", "cv.lst", "beta.init", "tau.init")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1), envir = backend_env, inherits = FALSE)]
if (length(missing_objects) > 0) {
  stop(sprintf("Missing required objects after load step: %s", paste(missing_objects, collapse = ", ")), call. = FALSE)
}

Y_data <- get("Y", envir = backend_env, inherits = FALSE)
X_data <- get("X", envir = backend_env, inherits = FALSE)
S_data <- get("S", envir = backend_env, inherits = FALSE)
cv_folds <- get("cv.lst", envir = backend_env, inherits = FALSE)
beta_init_default <- get("beta.init", envir = backend_env, inherits = FALSE)
tau_init_default <- get("tau.init", envir = backend_env, inherits = FALSE)

resolve_backend_fn <- function(fname) {
  if (exists(fname, envir = backend_env, mode = "function", inherits = FALSE)) {
    return(get(fname, envir = backend_env, inherits = FALSE))
  }
  if (exists(fname, envir = .GlobalEnv, mode = "function", inherits = TRUE)) {
    return(get(fname, envir = .GlobalEnv, inherits = TRUE))
  }
  NULL
}

run_mcmc <- if (backend == "ar2") resolve_backend_fn("mcmc_ar2") else resolve_backend_fn("mcmc")
if (is.null(run_mcmc)) {
  stop(sprintf("Required MCMC function not available for backend '%s'.", backend), call. = FALSE)
}

run_maxstable <- resolve_backend_fn("maxstable")

# Helper function to expand range expressions like "1:124" -> c("1", "2", ..., "124")
expand_ranges <- function(tokens) {
  result <- c()
  for (token in tokens) {
    if (grepl("^\\d+:\\d+$", token)) {
      # Match pattern like "1:124"
      parts <- as.integer(strsplit(token, ":", fixed = TRUE)[[1]])
      if (length(parts) == 2 && parts[1] <= parts[2]) {
        result <- c(result, as.character(seq(parts[1], parts[2])))
      } else {
        # Invalid range, keep original
        result <- c(result, token)
      }
    } else {
      result <- c(result, token)
    }
  }
  return(result)
}

args <- commandArgs(trailingOnly = TRUE)
env_settings <- trimws(Sys.getenv("US_ALL_SETTINGS", unset = ""))

if (nzchar(env_settings)) {
  setting_tokens <- unlist(strsplit(env_settings, ",", fixed = TRUE))
  setting_tokens <- trimws(setting_tokens)
  requested_settings <- expand_ranges(setting_tokens)
} else if (length(args) > 0) {
  args <- trimws(args)
  requested_settings <- expand_ranges(args)
} else {
  stop(
    "No settings provided. Use US_ALL_SETTINGS=114 or pass args like: 114 115",
    call. = FALSE
  )
}

requested_settings <- requested_settings[nzchar(requested_settings)]
if (length(requested_settings) == 0) {
  stop("No valid setting token found.", call. = FALSE)
}

if (any(!requested_settings %in% settings$setting_chr)) {
  unknown <- requested_settings[!requested_settings %in% settings$setting_chr]
  stop(sprintf("Unknown setting(s): %s", paste(unique(unknown), collapse = ", ")), call. = FALSE)
}

resolve_model <- function(method_raw) {
  method_lc <- tolower(trimws(method_raw))
  if (method_lc == "skew-t") {
    return(list(method = "t", skew = TRUE, is_maxstable = FALSE))
  }
  if (method_lc == "t") {
    return(list(method = "t", skew = FALSE, is_maxstable = FALSE))
  }
  if (method_lc == "gaussian") {
    return(list(method = "gaussian", skew = FALSE, is_maxstable = FALSE))
  }
  if (method_lc == "max-stable") {
    return(list(method = "max-stable", skew = FALSE, is_maxstable = TRUE))
  }
  stop(sprintf("Unsupported method in settings.csv: %s", method_raw), call. = FALSE)
}

pick_tau_init <- function(method_name) {
  if (tolower(method_name) == "gaussian") {
    return(tau_init_default)
  }
  0.05
}

parse_optional_int <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(NA_integer_)
  }

  txt <- trimws(as.character(x[1]))
  if (!nzchar(txt)) {
    return(NA_integer_)
  }

  val <- suppressWarnings(as.numeric(txt))
  if (!is.finite(val)) {
    return(NA_integer_)
  }

  as.integer(round(val))
}

format_runtime_timestamp <- function(timestamp) {
  if (length(timestamp) == 0 || anyNA(timestamp)) {
    return(NA_character_)
  }

  format(as.POSIXct(timestamp, tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

build_us_all_runtime_info <- function(
    started_at,
    elapsed_sec,
    outputfile,
    setting_token,
    setting_id,
    model_spec,
    nknots,
    threshold,
    backend,
    run_mode,
    use_cmaq,
    temporal,
    use_mrts,
    mrts_k,
    fold_elapsed_sec
) {
  finished_at <- started_at + as.numeric(elapsed_sec)
  completed_folds <- sum(is.finite(fold_elapsed_sec))

  list(
    schema_version = 1L,
    runner = "us-all-run.R",
    output_file = outputfile,
    setting_token = setting_token,
    setting_id = as.integer(setting_id),
    method = as.character(model_spec$method),
    skew = isTRUE(model_spec$skew),
    is_maxstable = isTRUE(model_spec$is_maxstable),
    nknots = as.integer(nknots),
    threshold = unname(as.numeric(threshold)),
    backend = backend,
    run_mode = run_mode,
    use_cmaq = isTRUE(use_cmaq),
    temporal = isTRUE(temporal),
    use_mrts = isTRUE(use_mrts),
    mrts_k = if (is.na(mrts_k)) NA_integer_ else as.integer(mrts_k),
    started_at_utc = format_runtime_timestamp(started_at),
    finished_at_utc = format_runtime_timestamp(finished_at),
    elapsed_sec = unname(as.numeric(elapsed_sec)),
    completed_folds = completed_folds,
    fold_elapsed_sec = as.numeric(fold_elapsed_sec),
    average_fold_elapsed_sec = if (completed_folds > 0) mean(fold_elapsed_sec[is.finite(fold_elapsed_sec)]) else NA_real_,
    mcmc_control = list(
      iters = as.integer(mcmc_ctrl$iters),
      burn = as.integer(mcmc_ctrl$burn),
      update = as.integer(mcmc_ctrl$update),
      thin = 1L
    )
  )
}

# build_mrts_covariates() lives in mrts_basis.R, shared with the full-data
# map pipeline (us-all-full-204.R / predict-cincy-204.R).
source("mrts_basis.R")

run_one_setting <- function(setting_token) {
  row <- settings[settings$setting_chr == setting_token, , drop = FALSE][1, ]

  numeric_prefix <- sub("^([0-9]+).*$", "\\1", setting_token)
  setting_seed_base <- suppressWarnings(as.integer(numeric_prefix))
  if (is.na(setting_seed_base)) {
    stop(sprintf("Setting must begin with digits for seeding: %s", setting_token), call. = FALSE)
  }

  model_spec <- resolve_model(row$method)
  threshold <- suppressWarnings(as.numeric(row$thresh))
  nknots <- suppressWarnings(as.integer(row$knots))
  use_cmaq <- truthy(row$CMAQ)
  temporal <- truthy(row$TS)
  ar2_label <- truthy(row$ar2)
  mrts_k <- parse_optional_int(row$mrts)
  use_mrts <- !is.na(mrts_k) && mrts_k > 0

  if (is.na(threshold) || is.na(nknots)) {
    stop(sprintf("Invalid thresh/knots for setting %s", setting_token), call. = FALSE)
  }

  cat("\n==============================\n")
  cat("Setting:", setting_token, "| Backend:", backend, "\n")
  cat("RUN_MODE:", RUN_MODE, "| iters:", mcmc_ctrl$iters, "| burn:", mcmc_ctrl$burn, "| update:", mcmc_ctrl$update, "\n")
  cat("Method:", row$method, "| K:", nknots, "| Threshold:", threshold, "\n")
  cat(
    "CMAQ:", if (use_cmaq) "yes" else "no",
    "| TS:", if (temporal) "yes" else "no",
    "| AR2 row:", if (ar2_label) "yes" else "no",
    "| MRTS covariates:", if (use_mrts) paste0("k=", mrts_k) else "no",
    "\n"
  )

  if (model_spec$is_maxstable) {
    if (is.null(run_maxstable)) {
      stop("maxstable() is not loaded in current backend.", call. = FALSE)
    }
    if (!use_cmaq) {
      stop("max-stable requires CMAQ covariate (x[, , 2]).", call. = FALSE)
    }
    cat("Special setting lane: max-stable (setting 2 style)\n")
  }

  if (model_spec$is_maxstable && use_mrts) {
    stop("MRTS covariates are not currently supported for max-stable settings.", call. = FALSE)
  }

  alias_match <- regmatches(setting_token, regexec("^([0-9]+)a$", setting_token))[[1]]
  if (length(alias_match) == 2) {
    outputfile <- file.path(results_dir, sprintf("us-all-%s-a.RData", alias_match[2]))
  } else {
    outputfile <- file.path(results_dir, sprintf("us-all-%s.RData", setting_token))
  }

  fit <- vector(mode = "list", length = 2)

  runtime_started_at <- Sys.time()
  start <- proc.time()
  fold_elapsed_sec <- rep(NA_real_, length(fit))

  for (val in 1:2) {
    set.seed(setting_seed_base * 100 + val)

    cat("CV", val, "started\n")
    val.idx <- cv_folds[[val]]

    y.o <- Y_data[-val.idx, ]

    x.o <- X_data[-val.idx, , , drop = FALSE]
    x.p <- X_data[val.idx, , , drop = FALSE]

    if (use_cmaq) {
      X.o <- x.o
      X.p <- x.p
    } else {
      X.o <- x.o[, , 1, drop = FALSE]
      X.p <- x.p[, , 1, drop = FALSE]
    }

    S.o <- S_data[-val.idx, ]
    S.p <- S_data[val.idx, ]

    if (use_mrts) {
      mrts_cov <- build_mrts_covariates(
        S_train = S.o,
        S_pred = S.p,
        nt = dim(X.o)[2],
        k = mrts_k
      )

      p_base <- dim(X.o)[3]
      p_mrts <- dim(mrts_cov$train)[3]

      X.o.ext <- array(NA_real_, dim = c(dim(X.o)[1], dim(X.o)[2], p_base + p_mrts))
      X.p.ext <- array(NA_real_, dim = c(dim(X.p)[1], dim(X.p)[2], p_base + p_mrts))

      X.o.ext[, , seq_len(p_base)] <- X.o
      X.p.ext[, , seq_len(p_base)] <- X.p

      for (j in seq_len(p_mrts)) {
        X.o.ext[, , p_base + j] <- mrts_cov$train[, , j]
        X.p.ext[, , p_base + j] <- mrts_cov$pred[, , j]
      }

      X.o <- X.o.ext
      X.p <- X.p.ext

      cat(
        "CV", val,
        "| MRTS source:", mrts_cov$source,
        "| cols(original/kept/dropped)=",
        mrts_cov$original_cols, "/", mrts_cov$kept_cols, "/", mrts_cov$dropped_cols,
        "| total covariates:", dim(X.o)[3],
        "\n"
      )
    }

    tic.set <- proc.time()

    if (model_spec$is_maxstable) {
      fit[[val]] <- run_maxstable(
        y = t(y.o),
        x = t(X.o[, , 2]),
        s = S.o,
        sp = S.p,
        xp = t(X.p[, , 2]),
        thresh = threshold,
        knots = S.o,
        iters = mcmc_ctrl$iters,
        burn = mcmc_ctrl$burn,
        update = mcmc_ctrl$update,
        thin = 1
      )
    } else {
      call_common <- list(
        y = y.o,
        s = S.o,
        x = X.o,
        x.pred = X.p,
        s.pred = S.p,
        method = model_spec$method,
        skew = model_spec$skew,
        keep.knots = FALSE,
        thresh.all = threshold,
        thresh.quant = FALSE,
        nknots = nknots,
        iters = mcmc_ctrl$iters,
        burn = mcmc_ctrl$burn,
        update = mcmc_ctrl$update,
        iterplot = FALSE,
        beta.init = beta_init_default,
        tau.init = pick_tau_init(model_spec$method),
        gamma.init = 0.5,
        rho.init = 1,
        rho.upper = 5,
        nu.init = 0.5,
        nu.upper = 10,
        min.s = c(-2.25, -1.55),
        max.s = c(2.35, 1.30)
      )

      if (temporal) {
        call_common$temporaltau <- TRUE
        call_common$temporalw <- TRUE
        call_common$temporalz <- TRUE
      }

      if (backend == "ar2") {
        call_common$ar2_tau <- isTRUE(temporal)
        call_common$ar2_w <- isTRUE(temporal)
        call_common$ar2_z <- isTRUE(temporal)
      }

      fit[[val]] <- do.call(run_mcmc, call_common)
    }

    toc.set <- proc.time()
    time.set <- (toc.set - tic.set)[3]
    fold_elapsed_sec[val] <- unname(as.numeric(time.set))

    elap.time.val <- (proc.time() - start)[3]
    avg.time.val <- elap.time.val / val
    runtime_info <- build_us_all_runtime_info(
      started_at = runtime_started_at,
      elapsed_sec = elap.time.val,
      outputfile = outputfile,
      setting_token = setting_token,
      setting_id = setting_seed_base,
      model_spec = model_spec,
      nknots = nknots,
      threshold = threshold,
      backend = backend,
      run_mode = RUN_MODE,
      use_cmaq = use_cmaq,
      temporal = temporal,
      use_mrts = use_mrts,
      mrts_k = mrts_k,
      fold_elapsed_sec = fold_elapsed_sec
    )

    cat("CV", val, "finished. Fold sec:", round(time.set, 2), "| Avg sec/dataset:", round(avg.time.val, 2), "\n")
    save(fit, runtime_info, file = outputfile)
  }

  cat("Saved:", outputfile, "\n")
}

for (setting_token in requested_settings) {
  run_one_setting(setting_token)
}

cat("\nAll requested settings completed:\n")
cat(paste(requested_settings, collapse = ", "), "\n")
