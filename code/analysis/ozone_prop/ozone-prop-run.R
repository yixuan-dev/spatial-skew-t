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

# RUN_MODE: dev for smoke tests, prod for full experiments.
RUN_MODE <- tolower(Sys.getenv("OZONE_PROP_RUN_MODE", unset = "prod"))
if (!RUN_MODE %in% c("dev", "prod")) {
  stop("OZONE_PROP_RUN_MODE must be one of: dev, prod", call. = FALSE)
}

mcmc_ctrl <- switch(RUN_MODE,
  dev  = list(iters = 2000, burn = 1000, update = 200),
  prod = list(iters = 30000, burn = 25000, update = 500)
)

results_dir <- Sys.getenv("OZONE_PROP_RESULTS_DIR", unset = "results")
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  cat("Created results directory:", results_dir, "\n")
}

prop_cov_update_every <- as.integer(Sys.getenv("OZONE_PROP_COV_UPDATE_EVERY", unset = "1"))
if (is.na(prop_cov_update_every) || prop_cov_update_every < 1L) {
  stop("OZONE_PROP_COV_UPDATE_EVERY must be a positive integer.", call. = FALSE)
}

# Load prop backend + ozone setup data into a controlled env.
backend_env <- new.env(parent = globalenv())
sys.source("./prop_load.R", envir = backend_env, chdir = TRUE)

required_objects <- c("Y", "X", "S", "cv.lst", "beta.init", "tau.init")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1), envir = backend_env, inherits = FALSE)]
if (length(missing_objects) > 0) {
  stop(sprintf("Missing required objects after prop_load: %s", paste(missing_objects, collapse = ", ")), call. = FALSE)
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

run_mcmc <- resolve_backend_fn("mcmc")
if (is.null(run_mcmc)) {
  stop("prop mcmc() function not available after prop_load.", call. = FALSE)
}

# settings file: defaults to settings.csv in the working dir.
settings_path <- Sys.getenv("OZONE_PROP_SETTINGS_FILE", unset = "settings_prop.csv")
if (!file.exists(settings_path)) {
  stop(sprintf("settings file not found: %s", settings_path), call. = FALSE)
}
settings <- read.csv(settings_path, stringsAsFactors = FALSE)
settings$setting_chr <- trimws(as.character(settings$setting))

expand_ranges <- function(tokens) {
  result <- c()
  for (token in tokens) {
    if (grepl("^\\d+:\\d+$", token)) {
      parts <- as.integer(strsplit(token, ":", fixed = TRUE)[[1]])
      if (length(parts) == 2 && parts[1] <= parts[2]) {
        result <- c(result, as.character(seq(parts[1], parts[2])))
      } else {
        result <- c(result, token)
      }
    } else {
      result <- c(result, token)
    }
  }
  result
}

args <- commandArgs(trailingOnly = TRUE)
env_settings <- trimws(Sys.getenv("OZONE_PROP_SETTINGS", unset = ""))

if (nzchar(env_settings)) {
  setting_tokens <- trimws(unlist(strsplit(env_settings, ",", fixed = TRUE)))
  requested_settings <- expand_ranges(setting_tokens)
} else if (length(args) > 0) {
  requested_settings <- expand_ranges(trimws(args))
} else {
  stop("No settings provided. Pass args like: 1 2 3 or 1:5, or set OZONE_PROP_SETTINGS.",
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
    return(list(method = "t", skew = TRUE))
  }
  if (method_lc == "t") {
    return(list(method = "t", skew = FALSE))
  }
  if (method_lc == "gaussian") {
    return(list(method = "gaussian", skew = FALSE))
  }
  stop(sprintf("Prop backend does not support method: %s (allowed: gaussian, t, skew-t)", method_raw),
    call. = FALSE
  )
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

build_runtime_info <- function(started_at, elapsed_sec, outputfile,
                               setting_token, setting_id, model_spec,
                               nknots, threshold, run_mode, use_cmaq,
                               prop_k, fold_elapsed_sec) {
  finished_at <- started_at + as.numeric(elapsed_sec)
  completed_folds <- sum(is.finite(fold_elapsed_sec))
  list(
    schema_version = 1L,
    runner = "ozone-prop-run.R",
    output_file = outputfile,
    setting_token = setting_token,
    setting_id = as.integer(setting_id),
    method = as.character(model_spec$method),
    skew = isTRUE(model_spec$skew),
    nknots = as.integer(nknots),
    threshold = unname(as.numeric(threshold)),
    backend = "prop",
    run_mode = run_mode,
    use_cmaq = isTRUE(use_cmaq),
    prop_k = as.integer(prop_k),
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
      thin = 1L,
      cov_update_every = as.integer(prop_cov_update_every)
    )
  )
}

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
  prop_k <- parse_optional_int(row$mrts)

  if (is.na(threshold) || is.na(nknots)) {
    stop(sprintf("Invalid thresh/knots for setting %s", setting_token), call. = FALSE)
  }
  if (is.na(prop_k) || prop_k <= 0) {
    stop(sprintf(
      "Setting %s requires a positive integer in the `mrts` column (used as prop_k).",
      setting_token
    ), call. = FALSE)
  }
  if (temporal) {
    stop(sprintf(
      "Setting %s has TS=yes; prop backend does not support temporal blocks.",
      setting_token
    ), call. = FALSE)
  }
  if (ar2_label) {
    stop(sprintf(
      "Setting %s has ar2=yes; prop backend does not support AR2 blocks.",
      setting_token
    ), call. = FALSE)
  }
  if (model_spec$method == "gaussian" && (model_spec$skew || nknots != 1L)) {
    stop(sprintf(
      "Setting %s: gaussian prop only supports skew=FALSE, knots=1.",
      setting_token
    ), call. = FALSE)
  }

  cat("\n==============================\n")
  cat("Setting:", setting_token, "| Backend: prop\n")
  cat(
    "RUN_MODE:", RUN_MODE, "| iters:", mcmc_ctrl$iters,
    "| burn:", mcmc_ctrl$burn, "| update:", mcmc_ctrl$update, "\n"
  )
  cat(
    "Method:", row$method, "| K:", nknots, "| Threshold:", threshold,
    "| CMAQ:", if (use_cmaq) "yes" else "no",
    "| prop_k:", prop_k,
    "| cov_update_every:", prop_cov_update_every, "\n"
  )

  alias_match <- regmatches(setting_token, regexec("^([0-9]+)a$", setting_token))[[1]]
  if (length(alias_match) == 2) {
    outputfile <- file.path(results_dir, sprintf("ozone-prop-%s-a.RData", alias_match[2]))
  } else {
    outputfile <- file.path(results_dir, sprintf("ozone-prop-%s.RData", setting_token))
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

    tic.set <- proc.time()
    fit[[val]] <- run_mcmc(
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
      max.s = c(2.35, 1.30),
      prop_k = as.integer(prop_k),
      prop_cov_update_every = as.integer(prop_cov_update_every)
    )

    toc.set <- proc.time()
    time.set <- (toc.set - tic.set)[3]
    fold_elapsed_sec[val] <- unname(as.numeric(time.set))

    elap.time.val <- (proc.time() - start)[3]
    avg.time.val <- elap.time.val / val
    runtime_info <- build_runtime_info(
      started_at = runtime_started_at,
      elapsed_sec = elap.time.val,
      outputfile = outputfile,
      setting_token = setting_token,
      setting_id = setting_seed_base,
      model_spec = model_spec,
      nknots = nknots,
      threshold = threshold,
      run_mode = RUN_MODE,
      use_cmaq = use_cmaq,
      prop_k = prop_k,
      fold_elapsed_sec = fold_elapsed_sec
    )

    cat(
      "CV", val, "finished. Fold sec:", round(time.set, 2),
      "| Avg sec/fold:", round(avg.time.val, 2), "\n"
    )
    save(fit, runtime_info, file = outputfile)
  }

  cat("Saved:", outputfile, "\n")
}

for (setting_token in requested_settings) {
  run_one_setting(setting_token)
}

cat("\nAll requested settings completed:\n")
cat(paste(requested_settings, collapse = ", "), "\n")
