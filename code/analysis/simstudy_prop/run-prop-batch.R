#########################################################################
# Execute one or more rows from settings_prop.csv.
# Row numbering is 1-based over data rows (header excluded).
#
# Usage:
#   Rscript run-prop-batch.R [rows]
#   Rscript run-prop-batch.R --rows=1:3
#   Rscript run-prop-batch.R --run_id=phase3-02
#
# If no selector is provided, runs all enabled rows ordered by priority.
#########################################################################

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = FALSE)
  script_dir <- dirname(script_path)
  if (dir.exists(script_dir)) {
    setwd(script_dir)
  }
}

source("./prop_simstudy_helpers.R", chdir = TRUE)

manifest <- load_settings_prop_manifest("settings_prop.csv")
args <- commandArgs(trailingOnly = TRUE)

rows_spec <- NULL
run_id_spec <- NULL
positional <- character(0)

for (arg in args) {
  if (grepl("^--rows=", arg)) {
    rows_spec <- sub("^--rows=", "", arg)
  } else if (grepl("^--run_id=", arg)) {
    run_id_spec <- sub("^--run_id=", "", arg)
  } else {
    positional <- c(positional, arg)
  }
}

if (length(positional) > 1L) {
  stop("Provide at most one positional rows specification.", call. = FALSE)
}
if (length(positional) == 1L) {
  if (!is.null(rows_spec)) {
    stop("Use either positional rows or --rows=..., not both.", call. = FALSE)
  }
  rows_spec <- positional[1]
}
if (!is.null(rows_spec) && !is.null(run_id_spec)) {
  stop("Use either row selection or run_id selection, not both.", call. = FALSE)
}

selected <- manifest[0, , drop = FALSE]
selection_mode <- "enabled"

if (!is.null(rows_spec)) {
  row_ids <- parse_manifest_rows_spec(rows_spec, nrow(manifest))
  selected <- manifest[manifest$row_id %in% row_ids, , drop = FALSE]
  selection_mode <- sprintf("rows=%s", paste(row_ids, collapse = ","))
} else if (!is.null(run_id_spec)) {
  run_ids <- parse_run_id_spec(run_id_spec)
  selected <- manifest[manifest$run_id %in% run_ids, , drop = FALSE]
  missing_ids <- setdiff(run_ids, selected$run_id)
  if (length(missing_ids) > 0L) {
    stop(sprintf("Unknown run_id values: %s", paste(missing_ids, collapse = ", ")), call. = FALSE)
  }
  selection_mode <- sprintf("run_id=%s", paste(run_ids, collapse = ","))
} else {
  selected <- manifest[manifest$enabled_flag, , drop = FALSE]
}

if (nrow(selected) == 0L) {
  stop("No rows selected from settings_prop.csv.", call. = FALSE)
}

selected <- selected[order(selected$priority, selected$row_id), , drop = FALSE]

output_results_dir <- file.path("output", "results")
if (!dir.exists(output_results_dir)) {
  dir.create(output_results_dir, recursive = TRUE, showWarnings = FALSE)
}

write.csv(
  selected[, c(
    "row_id", "run_id", "status", "enabled", "priority", "runner_script",
    "setting", "datasets", "methods", "prop_k", "workers", "ms_threads",
    "iters", "burn", "update", "thin", "prop_cov_update_every",
    "compare_after", "note"
  )],
  file.path(output_results_dir, "settings_prop_batch_selection.csv"),
  row.names = FALSE
)

rscript_bin <- resolve_rscript_binary()
status_rows <- vector("list", nrow(selected))

for (ii in seq_len(nrow(selected))) {
  spec <- selected[ii, , drop = FALSE]
  runner_script <- trimws(spec$runner_script[1])
  if (!nzchar(runner_script)) {
    runner_script <- "run-prop.R"
  }
  if (!file.exists(runner_script)) {
    stop(sprintf("Runner script not found for row %d: %s", spec$row_id[1], runner_script), call. = FALSE)
  }

  env_fit <- c(
    SIMSTUDY_PROP_ITERS = as.character(spec$iters[1]),
    SIMSTUDY_PROP_BURN = as.character(spec$burn[1]),
    SIMSTUDY_PROP_UPDATE = as.character(spec$update[1]),
    SIMSTUDY_PROP_THIN = as.character(spec$thin[1]),
    SIMSTUDY_PROP_COV_UPDATE_EVERY = as.character(spec$prop_cov_update_every[1]),
    SIMSTUDY_PROP_FITS_DIR = "fits"
  )

  fit_args <- c(
    sprintf("--setting=%s", spec$setting[1]),
    spec$datasets[1],
    as.character(spec$workers[1]),
    as.character(spec$ms_threads[1]),
    spec$methods[1],
    spec$prop_k[1]
  )

  cat(sprintf(
    "[batch] start row=%d run_id=%s setting=%s datasets=%s methods=%s prop_k=%s\n",
    spec$row_id[1], spec$run_id[1], spec$setting[1], spec$datasets[1], spec$methods[1], spec$prop_k[1]
  ))
  fit_status <- run_rscript_source_with_env(
    rscript_bin = rscript_bin,
    script = runner_script,
    args = fit_args,
    env = env_fit
  )
  if (as.integer(fit_status) != 0L) {
    stop(sprintf("Batch fit failed for row %d (run_id=%s).", spec$row_id[1], spec$run_id[1]), call. = FALSE)
  }

  compare_status <- NA_integer_
  if (isTRUE(spec$compare_after_flag[1])) {
    env_compare <- c(
      SIMSTUDY_PROP_FITS_DIR = "fits",
      SIMSTUDY_PROP_RESULTS_SETTINGS = as.character(spec$setting[1]),
      SIMSTUDY_PROP_RESULTS_DATASETS = spec$datasets[1],
      SIMSTUDY_PROP_RESULTS_METHODS = spec$methods[1],
      SIMSTUDY_PROP_RESULTS_K = spec$prop_k[1]
    )
    compare_status <- run_rscript_source_with_env(
      rscript_bin = rscript_bin,
      script = "results-prop.R",
      env = env_compare
    )
    if (as.integer(compare_status) != 0L) {
      stop(sprintf("Comparison step failed for row %d (run_id=%s).", spec$row_id[1], spec$run_id[1]), call. = FALSE)
    }
  }

  status_rows[[ii]] <- data.frame(
    row_id = spec$row_id[1],
    run_id = spec$run_id[1],
    selection_mode = selection_mode,
    runner_script = runner_script,
    setting = spec$setting[1],
    datasets = spec$datasets[1],
    methods = spec$methods[1],
    prop_k = spec$prop_k[1],
    compare_after = spec$compare_after[1],
    fit_status = fit_status,
    compare_status = compare_status,
    stringsAsFactors = FALSE
  )
}

batch_status <- do.call(rbind, status_rows)
write.csv(
  batch_status,
  file.path(output_results_dir, "settings_prop_batch_status.csv"),
  row.names = FALSE
)

cat(sprintf(
  "[batch] complete rows=%s output=%s\n",
  paste(selected$row_id, collapse = ","),
  file.path(output_results_dir, "settings_prop_batch_status.csv")
))
