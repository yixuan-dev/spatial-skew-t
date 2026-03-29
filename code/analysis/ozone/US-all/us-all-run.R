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

mcmc_ctrl <- switch(
  RUN_MODE,
  dev = list(iters = 2000, burn = 1000, update = 200),
  prod = list(iters = 30000, burn = 25000, update = 500)
)

results_dir <- trimws(Sys.getenv("US_ALL_RESULTS_DIR", unset = "results_new"))
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  cat("Created results directory:", results_dir, "\n")
}

backend <- tolower(Sys.getenv("US_ALL_MCMC_BACKEND", unset = "legacy"))
if (!backend %in% c("legacy", "ar2")) {
  stop("US_ALL_MCMC_BACKEND must be one of: legacy, ar2", call. = FALSE)
}

if (backend == "ar2") {
  source("./ar2_load.R", chdir = TRUE)
} else {
  source("./package_load.R", chdir = TRUE)
}

settings <- read.csv("settings.csv", stringsAsFactors = FALSE)
settings$setting_chr <- trimws(as.character(settings$setting))

required_objects <- c("Y", "X", "S", "cv.lst", "beta.init", "tau.init")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1), inherits = TRUE)]
if (length(missing_objects) > 0) {
  stop(sprintf("Missing required objects after load step: %s", paste(missing_objects, collapse = ", ")), call. = FALSE)
}

Y_data <- Y
X_data <- X
S_data <- S
cv_folds <- cv.lst
beta_init_default <- beta.init
tau_init_default <- tau.init

run_mcmc <- if (backend == "ar2") mcmc_ar2 else mcmc
run_maxstable <- if (exists("maxstable", mode = "function")) maxstable else NULL

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

  if (is.na(threshold) || is.na(nknots)) {
    stop(sprintf("Invalid thresh/knots for setting %s", setting_token), call. = FALSE)
  }

  cat("\n==============================\n")
  cat("Setting:", setting_token, "| Backend:", backend, "\n")
  cat("RUN_MODE:", RUN_MODE, "| iters:", mcmc_ctrl$iters, "| burn:", mcmc_ctrl$burn, "| update:", mcmc_ctrl$update, "\n")
  cat("Method:", row$method, "| K:", nknots, "| Threshold:", threshold, "\n")
  cat("CMAQ:", if (use_cmaq) "yes" else "no", "| TS:", if (temporal) "yes" else "no", "| AR2 row:", if (ar2_label) "yes" else "no", "\n")

  if (model_spec$is_maxstable) {
    if (is.null(run_maxstable)) {
      stop("maxstable() is not loaded in current backend.", call. = FALSE)
    }
    if (!use_cmaq) {
      stop("max-stable requires CMAQ covariate (x[, , 2]).", call. = FALSE)
    }
    cat("Special setting lane: max-stable (setting 2 style)\n")
  }

  outputfile <- file.path(results_dir, sprintf("us-all-%s.RData", setting_token))
  fit <- vector(mode = "list", length = 2)

  start <- proc.time()

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

    elap.time.val <- (proc.time() - start)[3]
    avg.time.val <- elap.time.val / val

    cat("CV", val, "finished. Fold sec:", round(time.set, 2), "| Avg sec/dataset:", round(avg.time.val, 2), "\n")
    save(fit, file = outputfile)
  }

  cat("Saved:", outputfile, "\n")
}

for (setting_token in requested_settings) {
  run_one_setting(setting_token)
}

cat("\nAll requested settings completed:\n")
cat(paste(requested_settings, collapse = ", "), "\n")
