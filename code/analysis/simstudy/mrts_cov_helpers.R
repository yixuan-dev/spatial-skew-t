set_working_dir_to_script <- function() {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_arg) == 0) {
    return(invisible(NULL))
  }

  script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = FALSE)
  script_dir <- dirname(script_path)
  if (dir.exists(script_dir)) {
    setwd(script_dir)
  }

  invisible(NULL)
}

parse_index_expr <- function(expr_str, arg_name) {
  expr <- trimws(expr_str)
  if (is.na(expr) || !nzchar(expr)) {
    stop(sprintf("%s specification cannot be empty", arg_name), call. = FALSE)
  }

  if (grepl("^\\(.*\\)$", expr)) {
    expr <- paste0("c", expr)
  }

  values_raw <- tryCatch(
    eval(parse(text = expr), envir = baseenv()),
    error = function(e) {
      stop(sprintf(
        "Invalid %s expression '%s'. Use forms like 1:5, c(1,2,6), (1,2,6)",
        arg_name, expr_str
      ), call. = FALSE)
    }
  )

  if (!is.numeric(values_raw) || length(values_raw) == 0) {
    stop(sprintf("%s must evaluate to a non-empty numeric vector", arg_name), call. = FALSE)
  }

  values_int <- as.integer(values_raw)
  if (any(is.na(values_int)) || any(values_raw != values_int)) {
    stop(sprintf("%s must contain integers only", arg_name), call. = FALSE)
  }

  sort(unique(values_int))
}

parse_env_indices <- function(env_name, default_ids, upper_bound) {
  raw <- trimws(Sys.getenv(env_name, unset = ""))
  if (!nzchar(raw)) {
    return(default_ids)
  }

  ids <- parse_index_expr(raw, env_name)
  if (any(ids < 1 | ids > upper_bound)) {
    stop(sprintf("%s must stay within 1..%d", env_name, upper_bound), call. = FALSE)
  }
  ids
}

parse_positive_int_expr <- function(expr_str, arg_name) {
  values <- parse_index_expr(expr_str, arg_name)
  if (any(values < 1)) {
    stop(sprintf("%s must contain positive integers only", arg_name), call. = FALSE)
  }
  values
}

resolve_mrts_k_values <- function(mrts_k_values = NULL, default_ids = c(5L, 10L, 15L)) {
  if (is.null(mrts_k_values)) {
    raw <- trimws(Sys.getenv("SIMSTUDY_MRTS_K", unset = ""))
    if (!nzchar(raw)) {
      return(as.integer(default_ids))
    }
    mrts_k_values <- parse_positive_int_expr(raw, "SIMSTUDY_MRTS_K")
  }

  values <- sort(unique(as.integer(mrts_k_values)))
  if (length(values) == 0 || any(is.na(values)) || any(values < 1)) {
    stop("MRTS k must be a non-empty set of positive integers.", call. = FALSE)
  }
  values
}

format_mrts_tag <- function(mrts_k) {
  sprintf("mrts%d", as.integer(mrts_k))
}

format_mrts_file_tag <- function(mrts_k) {
  sprintf("K%d", as.integer(mrts_k))
}

build_simstudy_result_file <- function(results_dir, setting_id, method_id, dataset_id, mrts_k = NA_integer_) {
  tag <- if (is.na(mrts_k)) "" else paste0("-", format_mrts_file_tag(mrts_k))
  file.path(
    results_dir,
    sprintf("%d-%d-%d%s.RData", as.integer(setting_id), as.integer(method_id), as.integer(dataset_id), tag)
  )
}

format_runtime_timestamp <- function(timestamp) {
  if (length(timestamp) == 0 || anyNA(timestamp)) {
    return(NA_character_)
  }

  format(as.POSIXct(timestamp, tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

build_simstudy_runtime_info <- function(
    started_at,
    elapsed_sec,
    outputfile,
    setting_id,
    dataset_id,
    runner,
    method_id = NA_integer_,
    method_key = NA_character_,
    mrts_k = NA_integer_,
    backend = NA_character_,
    control = NULL
) {
  finished_at <- started_at + as.numeric(elapsed_sec)

  list(
    schema_version = 1L,
    runner = as.character(runner),
    output_file = outputfile,
    setting_id = as.integer(setting_id),
    dataset_id = as.integer(dataset_id),
    method_id = if (is.na(method_id)) NA_integer_ else as.integer(method_id),
    method_key = if (is.na(method_key)) NA_character_ else as.character(method_key),
    mrts_k = if (is.na(mrts_k)) NA_integer_ else as.integer(mrts_k),
    backend = if (length(backend) == 0 || is.na(backend) || !nzchar(backend)) NA_character_ else as.character(backend),
    started_at_utc = format_runtime_timestamp(started_at),
    finished_at_utc = format_runtime_timestamp(finished_at),
    elapsed_sec = unname(as.numeric(elapsed_sec)),
    control = control
  )
}

get_simstudy_seed <- function(setting_id, method_id, dataset_id, mrts_k = NA_integer_) {
  if (is.na(mrts_k) && as.integer(method_id) == 6L) {
    return(as.integer(setting_id) * 100L + as.integer(dataset_id))
  }

  as.integer(method_id) * 1000L + as.integer(setting_id) * 100L + as.integer(dataset_id)
}

get_simstudy_method_catalog <- function(include_maxstable = TRUE) {
  catalog <- data.frame(
    method_id = 1:6,
    method_key = as.character(1:6),
    family = c(
      "gaussian",
      "skew_t_k1_thresh0",
      "t_k1_thresh_q80",
      "skew_t_k5_thresh0",
      "t_k5_thresh_q80",
      "max_stable_thresh_q80"
    ),
    label = c(
      "Gaussian",
      "Skew-t, K=1",
      "t, K=1, q(0.80)",
      "Skew-t, K=5",
      "t, K=5, q(0.80)",
      "Max-stable, q(0.80)"
    ),
    runner = c(rep("mcmc", 5), "maxstable"),
    method = c("gaussian", "t", "t", "t", "t", NA_character_),
    skew = c(FALSE, TRUE, FALSE, TRUE, FALSE, NA),
    thresh_all = c(0, 0, 0.80, 0, 0.80, NA),
    thresh_quant = c(TRUE, TRUE, TRUE, TRUE, TRUE, NA),
    nknots = c(1L, 1L, 1L, 5L, 5L, NA_integer_),
    use_exponential_fallback = c(FALSE, FALSE, TRUE, FALSE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )

  if (!isTRUE(include_maxstable)) {
    catalog <- catalog[catalog$method_id <= 5L, , drop = FALSE]
  }

  rownames(catalog) <- NULL
  catalog
}

resolve_mrts_method_ids <- function(method_ids = NULL, default_ids = 1:5) {
  if (is.null(method_ids)) {
    raw <- trimws(Sys.getenv("SIMSTUDY_MRTS_METHODS", unset = ""))
    if (!nzchar(raw)) {
      return(as.integer(default_ids))
    }
    method_ids <- parse_positive_int_expr(raw, "SIMSTUDY_MRTS_METHODS")
  }

  ids <- sort(unique(as.integer(method_ids)))
  if (length(ids) == 0 || any(is.na(ids)) || any(ids < 1L | ids > 5L)) {
    stop("MRTS methods must be a non-empty set of integers in 1..5.", call. = FALSE)
  }

  ids
}

resolve_simstudy_backend <- function() {
  backend <- tolower(trimws(Sys.getenv("SIMSTUDY_MCMC_BACKEND", unset = "legacy")))
  if (!backend %in% c("legacy", "ar2")) {
    stop("SIMSTUDY_MCMC_BACKEND must be one of: legacy, ar2", call. = FALSE)
  }
  backend
}

load_simstudy_backend <- function(backend = resolve_simstudy_backend()) {
  loader_file <- if (backend == "ar2") "./ar2_load.R" else "./package_load.R"
  backend_env <- new.env(parent = globalenv())
  sys.source(loader_file, envir = backend_env, chdir = TRUE)
  backend_env
}

resolve_backend_fn <- function(backend_env, fn_name) {
  if (exists(fn_name, envir = backend_env, mode = "function", inherits = FALSE)) {
    return(get(fn_name, envir = backend_env, inherits = FALSE))
  }
  if (exists(fn_name, envir = .GlobalEnv, mode = "function", inherits = TRUE)) {
    return(get(fn_name, envir = .GlobalEnv, inherits = TRUE))
  }
  NULL
}

build_mrts_covariates <- function(S_train, S_pred, nt, k) {
  if (is.na(k) || k <= 0) {
    stop("MRTS k must be a positive integer.", call. = FALSE)
  }

  safe_mrts_call <- function(fn, label) {
    train_obj <- tryCatch(fn(S_train, k = k), error = function(e) NULL)
    pred_obj <- tryCatch(fn(S_train, k = k, x = S_pred), error = function(e) NULL)

    if (!is.null(train_obj) && !is.null(pred_obj)) {
      return(list(train = as.matrix(train_obj), pred = as.matrix(pred_obj), source = label))
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
    stop("Unable to build MRTS covariates for this simulation run.", call. = FALSE)
  }

  basis_train <- as.matrix(basis_obj$train)
  basis_pred <- as.matrix(basis_obj$pred)

  if (nrow(basis_train) != nrow(S_train) || nrow(basis_pred) != nrow(S_pred)) {
    stop("MRTS basis row count mismatch.", call. = FALSE)
  }
  if (ncol(basis_train) != ncol(basis_pred)) {
    stop("MRTS basis column mismatch.", call. = FALSE)
  }

  original_cols <- ncol(basis_train)
  col_sd <- apply(basis_train, 2, sd)
  keep <- which(is.finite(col_sd) & col_sd > 1e-10)
  if (length(keep) == 0) {
    stop("All MRTS columns are near-constant after construction.", call. = FALSE)
  }

  basis_train <- basis_train[, keep, drop = FALSE]
  basis_pred <- basis_pred[, keep, drop = FALSE]

  mrts_train <- array(0, dim = c(nrow(S_train), nt, ncol(basis_train)))
  mrts_pred <- array(0, dim = c(nrow(S_pred), nt, ncol(basis_pred)))

  for (tt in seq_len(nt)) {
    mrts_train[, tt, ] <- basis_train
    mrts_pred[, tt, ] <- basis_pred
  }

  list(
    train = mrts_train,
    pred = mrts_pred,
    source = basis_obj$source,
    original_cols = original_cols,
    kept_cols = ncol(basis_train),
    dropped_cols = original_cols - ncol(basis_train)
  )
}

append_mrts_covariates <- function(X_train, X_pred, mrts_cov) {
  p_base <- dim(X_train)[3]
  p_mrts <- dim(mrts_cov$train)[3]

  X_train_ext <- array(NA_real_, dim = c(dim(X_train)[1], dim(X_train)[2], p_base + p_mrts))
  X_pred_ext <- array(NA_real_, dim = c(dim(X_pred)[1], dim(X_pred)[2], p_base + p_mrts))

  X_train_ext[, , seq_len(p_base)] <- X_train
  X_pred_ext[, , seq_len(p_base)] <- X_pred

  for (jj in seq_len(p_mrts)) {
    X_train_ext[, , p_base + jj] <- mrts_cov$train[, , jj]
    X_pred_ext[, , p_base + jj] <- mrts_cov$pred[, , jj]
  }

  list(train = X_train_ext, pred = X_pred_ext)
}

get_mrts_baseline_plan <- function(method_ids = NULL) {
  method_ids <- resolve_mrts_method_ids(method_ids)
  catalog <- get_simstudy_method_catalog(include_maxstable = FALSE)
  baseline <- catalog[catalog$method_id %in% method_ids, c(
    "method_id", "method_key", "family", "label", "method", "skew",
    "thresh_all", "thresh_quant", "nknots", "use_exponential_fallback"
  ), drop = FALSE]

  baseline$analysis_id <- baseline$method_id
  baseline$label <- paste0(baseline$label, ", no MRTS")
  baseline$mrts_k <- NA_integer_
  baseline$output_tag <- NA_character_

  baseline <- baseline[, c(
    "analysis_id", "method_id", "method_key", "family", "label", "method",
    "skew", "thresh_all", "thresh_quant", "nknots", "mrts_k",
    "output_tag", "use_exponential_fallback"
  ), drop = FALSE]
  baseline[order(baseline$method_id), , drop = FALSE]
}

get_mrts_analysis_plan <- function(mrts_k_values = NULL, method_ids = NULL) {
  mrts_k_values <- resolve_mrts_k_values(mrts_k_values)
  baseline <- get_mrts_baseline_plan(method_ids = method_ids)
  rows <- vector("list", length = nrow(baseline) * length(mrts_k_values))
  row_id <- 1L

  for (ii in seq_len(nrow(baseline))) {
    base_spec <- baseline[ii, , drop = FALSE]
    for (kk in mrts_k_values) {
      rows[[row_id]] <- data.frame(
        analysis_id = NA_integer_,
        method_id = base_spec$method_id[1],
        method_key = paste0(base_spec$method_id[1], "+", format_mrts_tag(kk)),
        family = base_spec$family[1],
        label = sprintf("%s, MRTS k=%d", sub(", no MRTS$", "", base_spec$label[1]), kk),
        method = base_spec$method[1],
        skew = isTRUE(base_spec$skew[1]),
        thresh_all = base_spec$thresh_all[1],
        thresh_quant = isTRUE(base_spec$thresh_quant[1]),
        nknots = base_spec$nknots[1],
        mrts_k = as.integer(kk),
        output_tag = format_mrts_tag(kk),
        baseline_analysis = base_spec$analysis_id[1],
        baseline_method_id = base_spec$method_id[1],
        baseline_method_key = base_spec$method_key[1],
        baseline_label = base_spec$label[1],
        use_exponential_fallback = isTRUE(base_spec$use_exponential_fallback[1]),
        stringsAsFactors = FALSE
      )
      row_id <- row_id + 1L
    }
  }

  plan <- do.call(rbind, rows)
  plan$analysis_id <- seq.int(
    from = max(get_simstudy_method_catalog(include_maxstable = TRUE)$method_id) + 1L,
    length.out = nrow(plan)
  )
  plan[order(plan$method_id, plan$mrts_k), , drop = FALSE]
}

run_mrts_pot_setting <- function(
    setting_id,
    dataset_ids = NULL,
    results_dir = trimws(Sys.getenv("SIMSTUDY_MRTS_RESULTS_DIR", unset = "results")),
    mrts_k_values = NULL,
    method_ids = NULL
) {
  backend <- resolve_simstudy_backend()
  backend_env <- load_simstudy_backend(backend = backend)
  run_mcmc <- resolve_backend_fn(backend_env, "mcmc")
  if (is.null(run_mcmc)) {
    stop(sprintf("MCMC function not available for backend '%s'.", backend), call. = FALSE)
  }

  y_data <- get("y", envir = backend_env, inherits = FALSE)
  x_data <- get("x", envir = backend_env, inherits = FALSE)
  s_data <- get("s", envir = backend_env, inherits = FALSE)
  ntest <- get("ntest", envir = backend_env, inherits = FALSE)

  if (is.null(dataset_ids)) {
    dataset_ids <- parse_env_indices("SIMSTUDY_MRTS_DATASETS", seq_len(dim(y_data)[3]), dim(y_data)[3])
  }

  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  }

  obs <- c(rep(TRUE, nrow(y_data) - ntest), rep(FALSE, ntest))
  method_ids <- resolve_mrts_method_ids(method_ids)
  mrts_k_values <- resolve_mrts_k_values(mrts_k_values)
  plan <- get_mrts_analysis_plan(mrts_k_values = mrts_k_values, method_ids = method_ids)
  write.csv(plan, "mrts_cov_analysis_plan.csv", row.names = FALSE)

  cat("Setting:", setting_id, "| backend:", backend, "\n")
  cat("Datasets:", paste(dataset_ids, collapse = ", "), "\n")
  cat("Methods:", paste(method_ids, collapse = ", "), "\n")
  cat("Planned MRTS analyses:", paste(plan$method_key, collapse = ", "), "\n")

  for (dataset_id in dataset_ids) {
    cat("start dataset", dataset_id, "\n")
    y_d <- y_data[, , dataset_id, setting_id]
    y_o <- y_d[obs, ]
    x_o <- x_data[obs, , , drop = FALSE]
    s_o <- s_data[obs, , drop = FALSE]
    x_p <- x_data[!obs, , , drop = FALSE]
    s_p <- s_data[!obs, , drop = FALSE]

    mrts_cache <- list()
    for (k in sort(unique(plan$mrts_k))) {
      mrts_cache[[as.character(k)]] <- build_mrts_covariates(
        S_train = s_o,
        S_pred = s_p,
        nt = dim(x_o)[2],
        k = k
      )
    }

    for (ii in seq_len(nrow(plan))) {
      spec <- plan[ii, , drop = FALSE]
      analysis_id <- spec$analysis_id[1]
      method_id <- spec$method_id[1]
      mrts_k <- spec$mrts_k[1]
      outputfile <- build_simstudy_result_file(
        results_dir = results_dir,
        setting_id = setting_id,
        method_id = method_id,
        dataset_id = dataset_id,
        mrts_k = mrts_k
      )

      set.seed(get_simstudy_seed(
        setting_id = setting_id,
        method_id = method_id,
        dataset_id = dataset_id,
        mrts_k = mrts_k
      ))

      extended_x <- append_mrts_covariates(
        X_train = x_o,
        X_pred = x_p,
        mrts_cov = mrts_cache[[as.character(mrts_k)]]
      )

      cat(
        "  start:", spec$label[1],
        "| MRTS source:", mrts_cache[[as.character(mrts_k)]]$source,
        "| cols(original/kept/dropped)=",
        mrts_cache[[as.character(mrts_k)]]$original_cols, "/",
        mrts_cache[[as.character(mrts_k)]]$kept_cols, "/",
        mrts_cache[[as.character(mrts_k)]]$dropped_cols,
        "\n"
      )

      started_at <- Sys.time()
      tic <- proc.time()
      if (isTRUE(spec$use_exponential_fallback[1])) {
        fit.1 <- tryCatch(
          run_mcmc(
            y = y_o,
            x = extended_x$train,
            s = s_o,
            s.pred = s_p,
            x.pred = extended_x$pred,
            method = spec$method[1],
            skew = isTRUE(spec$skew[1]),
            thresh.all = spec$thresh_all[1],
            thresh.quant = isTRUE(spec$thresh_quant[1]),
            nknots = spec$nknots[1],
            iterplot = FALSE,
            iters = 20000,
            burn = 10000,
            update = 1000,
            min.s = c(0, 0),
            max.s = c(10, 10),
            temporalw = FALSE,
            temporaltau = FALSE,
            temporalz = FALSE,
            rho.upper = 15,
            nu.upper = 10
          ),
          error = function(e) {
            run_mcmc(
              y = y_o,
              x = extended_x$train,
              s = s_o,
              s.pred = s_p,
              x.pred = extended_x$pred,
              method = spec$method[1],
              skew = isTRUE(spec$skew[1]),
              thresh.all = spec$thresh_all[1],
              thresh.quant = isTRUE(spec$thresh_quant[1]),
              nknots = spec$nknots[1],
              iterplot = FALSE,
              iters = 20000,
              burn = 10000,
              update = 1000,
              min.s = c(0, 0),
              max.s = c(10, 10),
              temporalw = FALSE,
              temporaltau = FALSE,
              temporalz = FALSE,
              rho.upper = 15,
              nu.upper = 10,
              cov.model = "exponential"
            )
          }
        )
      } else {
        fit.1 <- run_mcmc(
          y = y_o,
          x = extended_x$train,
          s = s_o,
          s.pred = s_p,
          x.pred = extended_x$pred,
          method = spec$method[1],
          skew = isTRUE(spec$skew[1]),
          thresh.all = spec$thresh_all[1],
          thresh.quant = isTRUE(spec$thresh_quant[1]),
          nknots = spec$nknots[1],
          iterplot = FALSE,
          iters = 20000,
          burn = 10000,
          update = 1000,
          min.s = c(0, 0),
          max.s = c(10, 10),
          temporalw = FALSE,
          temporaltau = FALSE,
          temporalz = FALSE,
          rho.upper = 15,
          nu.upper = 10
        )
      }

      elapsed_sec <- unname((proc.time() - tic)[3])
      cat("  end:", spec$label[1], "| elapsed:", round(elapsed_sec, 2), "sec\n")

      analysis_spec <- spec
      mrts_meta <- mrts_cache[[as.character(mrts_k)]]
      runtime_info <- build_simstudy_runtime_info(
        started_at = started_at,
        elapsed_sec = elapsed_sec,
        outputfile = outputfile,
        setting_id = setting_id,
        dataset_id = dataset_id,
        runner = "mrts_cov_helpers.R::run_mrts_pot_setting",
        method_id = method_id,
        method_key = spec$method_key[1],
        mrts_k = mrts_k,
        backend = backend,
        control = list(
          iters = 20000L,
          burn = 10000L,
          update = 1000L,
          thin = 1L
        )
      )
      save(fit.1, analysis_spec, mrts_meta, runtime_info, file = outputfile)

      rm(fit.1, analysis_spec, mrts_meta, runtime_info)
      gc()
    }
  }

  invisible(plan)
}

# ===================================================================
# Helpers for the scores.R / tables.R post-fit pipeline.
# Mirrors the analogous helpers in
# code/analysis/simstudy_prop/prop_simstudy_helpers.R.
# ===================================================================

# Derive a fits/results directory from the dataset filename:
#   simdata.RData      -> default_dir (e.g. "results")
#   simdata_def.RData  -> "<default_dir>_def"
#   other.RData        -> "<default_dir>_<other>"
derive_results_dir <- function(data_path, default_dir = "results") {
  base <- tools::file_path_sans_ext(basename(data_path))
  if (identical(base, "simdata")) {
    return(default_dir)
  }
  suffix <- sub("^simdata", "", base)
  if (!nzchar(suffix)) {
    return(paste0(default_dir, "_", base))
  }
  paste0(default_dir, suffix)
}

# Suffix appended to output filenames, derived from the dataset filename.
# Matches the tail of derive_results_dir():
#   simdata.RData      -> ""
#   simdata_def.RData  -> "_def"
#   other.RData        -> "_other"
derive_data_suffix <- function(data_path) {
  base <- tools::file_path_sans_ext(basename(data_path))
  if (identical(base, "simdata")) return("")
  s <- sub("^simdata", "", base)
  if (!nzchar(s)) paste0("_", base) else s
}

# Resolve a dataset path with sibling-fallback to ../simstudy/<basename>.
# data_arg = NULL or empty -> default lookup (./simdata.RData).
resolve_simstudy_data_path <- function(data_arg = NULL) {
  if (!is.null(data_arg) && nzchar(data_arg)) {
    if (file.exists(data_arg)) {
      return(normalizePath(data_arg, winslash = "/", mustWork = TRUE))
    }
    sibling <- file.path("..", "simstudy", basename(data_arg))
    if (file.exists(sibling)) {
      return(normalizePath(sibling, winslash = "/", mustWork = TRUE))
    }
    stop(
      sprintf("Data file not found: %s (also tried %s)", data_arg, sibling),
      call. = FALSE
    )
  }
  candidates <- c("./simdata.RData", "../simstudy/simdata.RData")
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0L) {
    stop(
      "Unable to locate simdata.RData (looked in ./ and ../simstudy/).",
      call. = FALSE
    )
  }
  normalizePath(existing[1], winslash = "/", mustWork = TRUE)
}

# Generic leading --key=value (or --key value) parser. Pulls one or more
# named flags from the front of args, stops at the first positional.
extract_leading_flags <- function(args, flag_names) {
  values <- setNames(vector("list", length(flag_names)), flag_names)
  while (length(args) > 0L && startsWith(args[1], "--")) {
    first <- args[1]
    matched <- FALSE
    for (name in flag_names) {
      flag_eq <- paste0("--", name, "=")
      flag_word <- paste0("--", name)
      if (startsWith(first, flag_eq)) {
        if (!is.null(values[[name]])) {
          stop(sprintf("Specify --%s only once", name), call. = FALSE)
        }
        values[[name]] <- sub(paste0("^", flag_eq), "", first)
        args <- if (length(args) > 1L) args[-1] else character(0)
        matched <- TRUE
        break
      } else if (first == flag_word) {
        if (!is.null(values[[name]])) {
          stop(sprintf("Specify --%s only once", name), call. = FALSE)
        }
        if (length(args) < 2L) {
          stop(sprintf("Missing value after --%s", name), call. = FALSE)
        }
        values[[name]] <- args[2]
        args <- if (length(args) > 2L) args[-c(1, 2)] else character(0)
        matched <- TRUE
        break
      }
    }
    if (!matched) {
      stop(
        sprintf(
          "Unknown leading flag: %s. Place known flags (%s) before positional args.",
          first,
          paste0("--", flag_names, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }
  list(args = args, values = values)
}
