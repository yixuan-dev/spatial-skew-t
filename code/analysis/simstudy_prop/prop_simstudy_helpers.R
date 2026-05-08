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
        "Invalid %s expression '%s'. Use forms like 1:5, c(1,2), (1,2).",
        arg_name, expr_str
      ), call. = FALSE)
    }
  )

  if (!is.numeric(values_raw) || length(values_raw) == 0L) {
    stop(sprintf("%s must evaluate to a non-empty numeric vector", arg_name), call. = FALSE)
  }

  values_int <- as.integer(values_raw)
  if (any(is.na(values_int)) || any(values_raw != values_int)) {
    stop(sprintf("%s must contain integers only", arg_name), call. = FALSE)
  }

  sort(unique(values_int))
}

parse_positive_int_expr <- function(expr_str, arg_name) {
  values <- parse_index_expr(expr_str, arg_name)
  if (any(values < 1L)) {
    stop(sprintf("%s must contain positive integers only", arg_name), call. = FALSE)
  }
  values
}

extract_setting_arg <- function(args) {
  if (length(args) == 0L) {
    return(list(args = character(0), setting = NULL))
  }

  has_setting_anywhere <- any(grepl("^--setting=", args)) || any(args == "--setting")
  setting_value <- NULL
  remaining_args <- args

  if (grepl("^--setting=", args[1])) {
    setting_value <- sub("^--setting=", "", args[1])
    remaining_args <- if (length(args) > 1L) args[-1] else character(0)
  } else if (args[1] == "--setting") {
    if (length(args) < 2L) {
      stop("Missing value after --setting", call. = FALSE)
    }
    setting_value <- args[2]
    remaining_args <- if (length(args) > 2L) args[-c(1, 2)] else character(0)
  } else if (has_setting_anywhere) {
    stop("--setting must be the first argument after script name", call. = FALSE)
  }

  if (any(grepl("^--setting=", remaining_args)) || any(remaining_args == "--setting")) {
    stop("Specify --setting only once and place it first", call. = FALSE)
  }

  list(args = remaining_args, setting = setting_value)
}

# Generalized leading-flag parser: pulls one or more --key=value (or
# --key value) flags from the front of args, in any order. Stops at the
# first positional argument. Used by run-prop.R for --data and --setting.
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

# Derive a results directory from the dataset filename:
#   simdata.RData      -> default_dir (e.g. "fits")
#   simdata_def.RData  -> "<default_dir>_def"
#   other.RData        -> "<default_dir>_<other>"
derive_prop_results_dir <- function(data_path, default_dir = "fits") {
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

parse_setting_spec <- function(setting_str, y_array) {
  setting_ids <- parse_index_expr(setting_str, "setting")
  if (length(setting_ids) != 1L) {
    stop("setting must evaluate to a single integer", call. = FALSE)
  }

  setting_id <- as.integer(setting_ids[[1]])
  max_setting <- dim(y_array)[4]
  if (setting_id < 1L || setting_id > max_setting) {
    stop(sprintf("setting must be in 1..%d", max_setting), call. = FALSE)
  }

  setting_id
}

parse_dataset_spec <- function(dataset_str) {
  dataset_ids <- parse_index_expr(dataset_str, "datasets")
  if (any(dataset_ids < 1L | dataset_ids > 50L)) {
    stop("datasets must be integers in 1..50", call. = FALSE)
  }
  dataset_ids
}

parse_prop_methods_spec <- function(methods_str) {
  method_ids <- parse_index_expr(methods_str, "methods")
  if (any(method_ids < 1L | method_ids > 5L)) {
    stop("prop methods must be integers in 1..5", call. = FALSE)
  }
  sort(unique(as.integer(method_ids)))
}

parse_prop_k_spec <- function(prop_str, default_k = 10L) {
  expr <- trimws(prop_str)
  if (!nzchar(expr)) {
    return(as.integer(default_k))
  }
  parse_positive_int_expr(expr, "mrts_k")
}

format_prop_tag <- function(prop_k) {
  sprintf("p%d", as.integer(prop_k))
}

format_prop_file_tag <- function(prop_k) {
  format_prop_tag(prop_k)
}

format_prop_legacy_file_tag <- function(prop_k) {
  sprintf("P%d", as.integer(prop_k))
}

build_prop_result_file <- function(results_dir, setting_id, method_id, dataset_id, prop_k) {
  file.path(
    results_dir,
    sprintf(
      "%d-%d-%d-%s.RData",
      as.integer(setting_id),
      as.integer(method_id),
      as.integer(dataset_id),
      format_prop_file_tag(prop_k)
    )
  )
}

build_prop_legacy_result_file <- function(results_dir, setting_id, method_id, dataset_id, prop_k) {
  file.path(
    results_dir,
    sprintf(
      "%d-%d-%d-%s.RData",
      as.integer(setting_id),
      as.integer(method_id),
      as.integer(dataset_id),
      format_prop_legacy_file_tag(prop_k)
    )
  )
}

build_prop_result_candidates <- function(results_dir, setting_id, method_id, dataset_id, prop_k) {
  unique(c(
    build_prop_result_file(results_dir, setting_id, method_id, dataset_id, prop_k),
    build_prop_legacy_result_file(results_dir, setting_id, method_id, dataset_id, prop_k)
  ))
}

format_runtime_timestamp <- function(timestamp) {
  if (length(timestamp) == 0L || anyNA(timestamp)) {
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
    method_id,
    method_key,
    prop_k,
    backend = "prop-simstudy",
    control = NULL) {
  finished_at <- started_at + as.numeric(elapsed_sec)
  list(
    schema_version = 1L,
    runner = as.character(runner),
    output_file = outputfile,
    setting_id = as.integer(setting_id),
    dataset_id = as.integer(dataset_id),
    method_id = as.integer(method_id),
    method_key = as.character(method_key),
    prop_k = as.integer(prop_k),
    backend = as.character(backend),
    started_at_utc = format_runtime_timestamp(started_at),
    finished_at_utc = format_runtime_timestamp(finished_at),
    elapsed_sec = unname(as.numeric(elapsed_sec)),
    control = control
  )
}

get_prop_seed <- function(setting_id, method_id, dataset_id, prop_k) {
  as.integer(method_id) * 100000L +
    as.integer(prop_k) * 1000L +
    as.integer(setting_id) * 100L +
    as.integer(dataset_id)
}

get_prop_method_catalog <- function() {
  data.frame(
    method_id = 1:5,
    method_key = as.character(1:5),
    family = c(
      "gaussian",
      "skew_t_k1_thresh0",
      "t_k1_thresh_q80",
      "skew_t_k5_thresh0",
      "t_k5_thresh_q80"
    ),
    label = c(
      "Gaussian",
      "Skew-t, K=1",
      "t, K=1, q(0.80)",
      "Skew-t, K=5",
      "t, K=5, q(0.80)"
    ),
    method = c("gaussian", "t", "t", "t", "t"),
    skew = c(FALSE, TRUE, FALSE, TRUE, FALSE),
    thresh_all = c(0, 0, 0.80, 0, 0.80),
    thresh_quant = c(TRUE, TRUE, TRUE, TRUE, TRUE),
    nknots = c(1L, 1L, 1L, 5L, 5L),
    prop_supported = rep(TRUE, 5L),
    stringsAsFactors = FALSE
  )
}

build_prop_run_plan <- function(method_ids, prop_k_values) {
  catalog <- get_prop_method_catalog()
  base <- catalog[catalog$method_id %in% method_ids, , drop = FALSE]
  combos <- expand.grid(
    method_id = base$method_id,
    prop_k = as.integer(prop_k_values),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  plan <- merge(combos, base, by = "method_id", sort = FALSE)
  plan$output_tag <- vapply(plan$prop_k, format_prop_tag, FUN.VALUE = character(1))
  plan$method_key <- paste0(plan$method_id, "+", plan$output_tag)
  plan$label <- sprintf("%s + prop(k=%d)", plan$label, plan$prop_k)
  plan$plan_id <- seq_len(nrow(plan))
  plan
}

normalize_yes_no_flag <- function(x, default = FALSE) {
  value <- tolower(trimws(as.character(x)))
  if (!nzchar(value) || is.na(value)) {
    return(isTRUE(default))
  }
  value %in% c("yes", "y", "true", "1")
}

load_settings_prop_manifest <- function(path = "settings_prop.csv") {
  if (!file.exists(path)) {
    stop(sprintf("settings_prop manifest not found: %s", path), call. = FALSE)
  }

  manifest <- read.csv(path, stringsAsFactors = FALSE)
  required_cols <- c(
    "run_id", "status", "enabled", "priority", "runner_script",
    "setting", "datasets", "methods", "prop_k", "workers", "ms_threads",
    "iters", "burn", "update", "thin", "prop_cov_update_every",
    "compare_after"
  )
  missing_cols <- setdiff(required_cols, names(manifest))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf(
        "settings_prop.csv is missing required columns: %s",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  manifest$row_id <- seq_len(nrow(manifest))
  manifest$enabled_flag <- vapply(manifest$enabled, normalize_yes_no_flag, logical(1))
  manifest$compare_after_flag <- vapply(manifest$compare_after, normalize_yes_no_flag, logical(1))
  manifest
}

parse_manifest_rows_spec <- function(row_str, n_rows) {
  row_ids <- parse_index_expr(row_str, "rows")
  if (any(row_ids < 1L | row_ids > n_rows)) {
    stop(sprintf("rows must be integers in 1..%d", n_rows), call. = FALSE)
  }
  as.integer(row_ids)
}

parse_run_id_spec <- function(run_id_str) {
  raw <- trimws(run_id_str)
  if (!nzchar(raw)) {
    stop("run_id specification cannot be empty", call. = FALSE)
  }

  ids <- trimws(unlist(strsplit(raw, ",", fixed = TRUE)))
  ids <- ids[nzchar(ids)]
  if (length(ids) == 0L) {
    stop("run_id specification cannot be empty", call. = FALSE)
  }

  unique(ids)
}

resolve_rscript_binary <- function() {
  candidates <- c(
    file.path(R.home("bin"), "Rscript.exe"),
    file.path(R.home("bin"), "Rscript"),
    Sys.which("Rscript")
  )
  existing <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (length(existing) == 0L) {
    stop("Unable to locate an Rscript binary for batch execution.", call. = FALSE)
  }
  normalizePath(existing[1], winslash = "/", mustWork = TRUE)
}

quote_rscript_string <- function(x) {
  encodeString(as.character(x), quote = "\"")
}

run_rscript_source_with_env <- function(rscript_bin, script, args = character(0), env = character(0)) {
  script_path <- normalizePath(script, winslash = "/", mustWork = TRUE)
  script_dir <- dirname(script_path)
  wrapper_path <- tempfile(pattern = "simstudy_prop_batch_", fileext = ".R")
  on.exit(unlink(wrapper_path), add = TRUE)

  env_full <- env
  env_full[["SIMSTUDY_PROP_SCRIPT_DIR"]] <- script_dir

  wrapper_lines <- character(0)
  if (length(env_full) > 0L) {
    env_names <- names(env_full)
    if (is.null(env_names) || any(!nzchar(env_names))) {
      stop("env must be a named character vector.", call. = FALSE)
    }
    env_expr <- paste(
      sprintf(
        "%s=%s",
        env_names,
        vapply(env_full, quote_rscript_string, FUN.VALUE = character(1))
      ),
      collapse = ", "
    )
    wrapper_lines <- c(wrapper_lines, sprintf("Sys.setenv(%s)", env_expr))
  }
  wrapper_lines <- c(wrapper_lines, sprintf("source(%s, chdir=TRUE)", quote_rscript_string(script_path)))
  writeLines(wrapper_lines, con = wrapper_path, useBytes = TRUE)

  call_args <- c(normalizePath(wrapper_path, winslash = "/", mustWork = TRUE))
  if (length(args) > 0L) {
    call_args <- c(call_args, as.character(args))
  }

  system2(rscript_bin, args = call_args)
}
