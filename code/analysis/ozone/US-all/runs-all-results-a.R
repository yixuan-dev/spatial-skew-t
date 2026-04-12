rm(list = ls())

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = FALSE)
  script_dir <- dirname(script_path)
  if (dir.exists(script_dir)) {
    setwd(script_dir)
  }
}

source("mrts_cov_helpers.R")

resolve_rscript_path <- function() {
  sys_r <- Sys.which("Rscript")
  local_r <- file.path(
    R.home("bin"),
    ifelse(.Platform$OS.type == "windows", "Rscript.exe", "Rscript")
  )

  candidates <- c(sys_r, local_r)
  candidates <- candidates[nzchar(candidates)]
  candidates <- unique(normalizePath(candidates, winslash = "/", mustWork = FALSE))
  hit <- candidates[file.exists(candidates)]

  if (length(hit) == 0) {
    stop("Cannot find Rscript executable.")
  }

  hit[1]
}

get_arg_value <- function(args, key, default = NA_character_) {
  prefix <- paste0("--", key, "=")
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (length(hit) > 0) {
    return(sub(prefix, "", hit[1]))
  }
  idx <- which(args == paste0("--", key))
  if (length(idx) > 0 && idx[1] < length(args)) {
    return(args[idx[1] + 1])
  }
  default
}

args <- commandArgs(trailingOnly = TRUE)
run_now <- "--run" %in% args
force_run <- "--force" %in% args

settings <- read.csv("settings.csv", stringsAsFactors = FALSE)
target_pairs <- build_mrts_target_pairs(settings)

baseline_settings <- sort(unique(target_pairs$baseline_setting))
mrts_settings <- sort(unique(target_pairs$mrts_setting))
target_settings <- sort(unique(c(1L, baseline_settings, mrts_settings)))

run_mode <- tolower(get_arg_value(args, "run-mode", default = Sys.getenv("US_ALL_RUN_MODE", unset = "dev")))
backend <- tolower(get_arg_value(args, "backend", default = Sys.getenv("US_ALL_MCMC_BACKEND", unset = "legacy")))
results_dir <- get_arg_value(args, "results-dir", default = Sys.getenv("US_ALL_RESULTS_DIR", unset = "results_mrts_cov_dev"))
if (!nzchar(results_dir)) {
  results_dir <- "results_mrts_cov_dev"
}

if (!run_mode %in% c("dev", "prod")) {
  stop("--run-mode must be one of: dev, prod")
}
if (!backend %in% c("legacy", "ar2")) {
  stop("--backend must be one of: legacy, ar2")
}

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
}

runner_file <- normalizePath("us-all-run.R", winslash = "/", mustWork = TRUE)
rscript_bin <- resolve_rscript_path()

settings_lookup <- settings
settings_lookup$setting_num <- suppressWarnings(as.integer(settings_lookup$setting))

build_role_table <- function(target_settings, baseline_settings, mrts_settings, settings_lookup) {
  role <- rep("other", length(target_settings))
  role[target_settings == 1] <- "gaussian_reference"
  role[target_settings %in% baseline_settings] <- "baseline_no_mrts"
  role[target_settings %in% mrts_settings] <- "mrts_candidate"

  idx <- match(target_settings, settings_lookup$setting_num)
  mrts_value <- if ("mrts" %in% names(settings_lookup)) settings_lookup$mrts[idx] else NA
  method <- if ("method" %in% names(settings_lookup)) settings_lookup$method[idx] else NA
  thresh <- if ("thresh" %in% names(settings_lookup)) settings_lookup$thresh[idx] else NA

  data.frame(
    setting = target_settings,
    role = role,
    method = method,
    thresh = thresh,
    mrts = mrts_value,
    stringsAsFactors = FALSE
  )
}

role_table <- build_role_table(target_settings, baseline_settings, mrts_settings, settings_lookup)
write.csv(target_pairs, "runs-all-results-a-target-pairs.csv", row.names = FALSE)
write.csv(role_table, "runs-all-results-a-target-settings.csv", row.names = FALSE)

cat("\n=== MRTS comparison target plan ===\n")
cat("Baseline settings:", paste(baseline_settings, collapse = ", "), "\n")
cat("MRTS settings:", paste(mrts_settings, collapse = ", "), "\n")
cat("All target settings:", paste(target_settings, collapse = ", "), "\n")
cat("results_dir:", results_dir, "\n")
cat("run_mode:", run_mode, "| backend:", backend, "\n")
cat("run_now:", run_now, "| force_run:", force_run, "\n")

run_one_setting <- function(setting_id) {
  result_file <- file.path(results_dir, sprintf("us-all-%d.RData", setting_id))
  start_time <- Sys.time()

  if (file.exists(result_file) && !force_run) {
    cat(sprintf("[%s] setting=%d | status=already_present\n", format(Sys.time(), "%H:%M:%S"), setting_id))
    return(data.frame(
      setting = setting_id,
      run_status = "already_present",
      result_exists = TRUE,
      elapsed_sec = NA_real_,
      checked_at = format(start_time, "%Y-%m-%d %H:%M:%S"),
      error_message = "",
      stringsAsFactors = FALSE
    ))
  }

  if (!run_now) {
    cat(sprintf("[%s] setting=%d | status=dry_missing\n", format(Sys.time(), "%H:%M:%S"), setting_id))
    return(data.frame(
      setting = setting_id,
      run_status = "dry_missing",
      result_exists = file.exists(result_file),
      elapsed_sec = NA_real_,
      checked_at = format(start_time, "%Y-%m-%d %H:%M:%S"),
      error_message = "",
      stringsAsFactors = FALSE
    ))
  }

  old_env <- Sys.getenv(c("US_ALL_RUN_MODE", "US_ALL_MCMC_BACKEND", "US_ALL_RESULTS_DIR"), unset = NA)
  on.exit({
    for (i in seq_along(old_env)) {
      nm <- names(old_env)[i]
      val <- old_env[[i]]
      if (is.na(val)) {
        Sys.unsetenv(nm)
      } else {
        do.call(Sys.setenv, setNames(list(val), nm))
      }
    }
  }, add = TRUE)

  Sys.setenv(US_ALL_RUN_MODE = run_mode)
  Sys.setenv(US_ALL_MCMC_BACKEND = backend)
  Sys.setenv(US_ALL_RESULTS_DIR = results_dir)

  tic <- proc.time()[3]
  cmd_out <- character(0)
  cmd_err <- ""
  exit_status <- 0L

  tryCatch(
    {
      cmd_out <- system2(
        command = rscript_bin,
        args = c(runner_file, as.character(setting_id)),
        stdout = TRUE,
        stderr = TRUE
      )
      st <- attr(cmd_out, "status")
      if (!is.null(st)) {
        exit_status <- as.integer(st)
      }
    },
    error = function(e) {
      cmd_err <<- conditionMessage(e)
    }
  )

  elapsed <- proc.time()[3] - tic
  result_exists <- file.exists(result_file)

  if (identical(cmd_err, "") && exit_status == 0L && result_exists) {
    run_status <- "ok"
    error_message <- ""
  } else {
    run_status <- if (result_exists) "error_result_exists" else "error"
    if (!nzchar(cmd_err)) {
      out_tail <- tail(cmd_out, 20)
      cmd_err <- paste(c(sprintf("Rscript exit status: %d", exit_status), out_tail), collapse = "\n")
    }
    error_message <- cmd_err
  }

  cat(sprintf("[%s] setting=%d | status=%s | result=%s\n", format(Sys.time(), "%H:%M:%S"), setting_id, run_status, ifelse(result_exists, "FOUND", "MISSING")))

  data.frame(
    setting = setting_id,
    run_status = run_status,
    result_exists = result_exists,
    elapsed_sec = elapsed,
    checked_at = format(start_time, "%Y-%m-%d %H:%M:%S"),
    error_message = error_message,
    stringsAsFactors = FALSE
  )
}

status_rows <- lapply(target_settings, run_one_setting)
status_df <- do.call(rbind, status_rows)
status_df <- merge(status_df, role_table, by = "setting", all.x = TRUE, sort = FALSE)

write.csv(status_df, "runs-all-results-a-status.csv", row.names = FALSE)

cat("\n=== runs-all-results-a summary ===\n")
print(table(status_df$run_status, useNA = "ifany"))

failed <- status_df$setting[status_df$run_status %in% c("error", "error_result_exists", "dry_missing")]
if (length(failed) > 0) {
  cat("Potentially unfinished settings:", paste(failed, collapse = ", "), "\n")
}

cat("Outputs written:\n")
cat("- runs-all-results-a-target-pairs.csv\n")
cat("- runs-all-results-a-target-settings.csv\n")
cat("- runs-all-results-a-status.csv\n")
