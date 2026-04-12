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

results_dir <- trimws(Sys.getenv("US_ALL_RESULTS_DIR", unset = "results_new"))
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  cat("Created results directory:", results_dir, "\n")
}

# Optional compatibility output:
# when enabled, setting tokens like "5a" will also be saved as "us-all-5-a.RData"
write_legacy_a_alias <- truthy(Sys.getenv("US_ALL_WRITE_LEGACY_A_ALIAS", unset = "false"))
if (write_legacy_a_alias) {
  cat("Legacy -a alias output: enabled\n")
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

build_mrts_covariates <- function(S_train, S_pred, nt, k) {
  if (is.na(k) || k <= 0) {
    stop("MRTS k must be a positive integer.", call. = FALSE)
  }

  safe_mrts_call <- function(fn, label) {
    train_obj <- tryCatch(fn(S_train, k = k), error = function(e) NULL)
    pred_obj <- tryCatch(fn(S_train, k = k, x = S_pred), error = function(e) NULL)

    if (!is.null(train_obj) && !is.null(pred_obj)) {
      train_mat <- as.matrix(train_obj)
      pred_mat <- as.matrix(pred_obj)
      return(list(train = train_mat, pred = pred_mat, source = label))
    }

    all_obj <- tryCatch(fn(rbind(S_train, S_pred), k = k), error = function(e) NULL)
    if (is.null(all_obj)) {
      return(NULL)
    }

    all_mat <- as.matrix(all_obj)
    if (nrow(all_mat) != nrow(S_train) + nrow(S_pred)) {
      return(NULL)
    }

    n_train <- nrow(S_train)
    list(
      train = all_mat[seq_len(n_train), , drop = FALSE],
      pred = all_mat[(n_train + 1):nrow(all_mat), , drop = FALSE],
      source = paste0(label, " [combined fallback]")
    )
  }

  basis_obj <- NULL
  if (requireNamespace("autoFRK", quietly = TRUE) && exists("mrts", where = asNamespace("autoFRK"), inherits = FALSE)) {
    basis_obj <- safe_mrts_call(autoFRK::mrts, "autoFRK::mrts")
  }

  if (is.null(basis_obj) && exists("mrts", mode = "function")) {
    basis_obj <- safe_mrts_call(mrts, "mrts()")
  }

  if (is.null(basis_obj)) {
    stop("Unable to build MRTS covariates for this setting (mrts function unavailable or failed).", call. = FALSE)
  }

  basis_train <- as.matrix(basis_obj$train)
  basis_pred <- as.matrix(basis_obj$pred)

  if (nrow(basis_train) != nrow(S_train) || nrow(basis_pred) != nrow(S_pred)) {
    stop("MRTS basis row count does not match train/pred locations.", call. = FALSE)
  }

  if (ncol(basis_train) != ncol(basis_pred)) {
    stop("MRTS basis column mismatch between train and pred matrices.", call. = FALSE)
  }

  original_cols <- ncol(basis_train)
  col_sd <- apply(basis_train, 2, sd)
  keep <- which(is.finite(col_sd) & col_sd > 1e-10)

  if (length(keep) == 0) {
    stop("All MRTS columns are near-constant after construction.", call. = FALSE)
  }

  basis_train_kept <- basis_train[, keep, drop = FALSE]
  basis_pred_kept <- basis_pred[, keep, drop = FALSE]

  mrts_train <- array(0, dim = c(nrow(S_train), nt, ncol(basis_train_kept)))
  mrts_pred <- array(0, dim = c(nrow(S_pred), nt, ncol(basis_pred_kept)))

  for (t in seq_len(nt)) {
    mrts_train[, t, ] <- basis_train_kept
    mrts_pred[, t, ] <- basis_pred_kept
  }

  list(
    train = mrts_train,
    pred = mrts_pred,
    source = basis_obj$source,
    original_cols = original_cols,
    kept_cols = ncol(basis_train_kept),
    dropped_cols = original_cols - ncol(basis_train_kept)
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

  outputfile <- file.path(results_dir, sprintf("us-all-%s.RData", setting_token))
  legacy_a_outputfile <- NULL
  if (write_legacy_a_alias) {
    alias_match <- regmatches(setting_token, regexec("^([0-9]+)a$", setting_token))[[1]]
    if (length(alias_match) == 2) {
      legacy_a_outputfile <- file.path(results_dir, sprintf("us-all-%s-a.RData", alias_match[2]))
    }
  }

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

    elap.time.val <- (proc.time() - start)[3]
    avg.time.val <- elap.time.val / val

    cat("CV", val, "finished. Fold sec:", round(time.set, 2), "| Avg sec/dataset:", round(avg.time.val, 2), "\n")
    save(fit, file = outputfile)
    if (!is.null(legacy_a_outputfile)) {
      save(fit, file = legacy_a_outputfile)
    }
  }

  cat("Saved:", outputfile, "\n")
  if (!is.null(legacy_a_outputfile)) {
    cat("Saved legacy alias:", legacy_a_outputfile, "\n")
  }
}

for (setting_token in requested_settings) {
  run_one_setting(setting_token)
}

cat("\nAll requested settings completed:\n")
cat(paste(requested_settings, collapse = ", "), "\n")
