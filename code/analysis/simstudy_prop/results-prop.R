rm(list = ls())

script_dir_override <- trimws(Sys.getenv("SIMSTUDY_PROP_SCRIPT_DIR", unset = ""))
if (nzchar(script_dir_override) && dir.exists(script_dir_override)) {
  setwd(normalizePath(script_dir_override, winslash = "/", mustWork = TRUE))
} else {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_arg) > 0) {
    script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = FALSE)
    script_dir <- dirname(script_path)
    if (dir.exists(script_dir)) {
      setwd(script_dir)
    }
  }
}

source("./prop_simstudy_helpers.R")

resolve_simstudy_prop_data <- function() {
  candidates <- c(
    "./simdata.RData",
    "../simstudy/simdata.RData"
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0L) {
    stop(
      "Unable to locate simdata.RData. Checked ./simdata.RData and ../simstudy/simdata.RData",
      call. = FALSE
    )
  }
  normalizePath(existing[1], winslash = "/", mustWork = TRUE)
}

load(resolve_simstudy_prop_data())
source("../../R/prop/auxfunctions.R")

build_baseline_result_file <- function(results_dir, setting_id, method_id, dataset_id) {
  file.path(
    results_dir,
    sprintf("%d-%d-%d.RData", as.integer(setting_id), as.integer(method_id), as.integer(dataset_id))
  )
}

parse_env_indices <- function(env_name, default_values, max_value) {
  raw <- trimws(Sys.getenv(env_name, unset = ""))
  if (!nzchar(raw)) {
    return(as.integer(default_values))
  }

  values <- parse_index_expr(raw, env_name)
  if (any(values < 1L | values > max_value)) {
    stop(sprintf("%s must be in 1..%d", env_name, max_value), call. = FALSE)
  }

  as.integer(values)
}

load_fit_object <- function(result_file) {
  if (!file.exists(result_file)) {
    return(NULL)
  }

  env <- new.env(parent = emptyenv())
  load(result_file, envir = env)

  if (exists("fit.1", envir = env, inherits = FALSE)) {
    return(get("fit.1", envir = env, inherits = FALSE))
  }
  if (exists("fit", envir = env, inherits = FALSE)) {
    return(get("fit", envir = env, inherits = FALSE))
  }

  NULL
}

resolve_prop_result_file <- function(results_dir, setting_id, method_id, dataset_id, prop_k) {
  candidates <- build_prop_result_candidates(
    results_dir = results_dir,
    setting_id = setting_id,
    method_id = method_id,
    dataset_id = dataset_id,
    prop_k = prop_k
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) > 0L) {
    return(existing[1])
  }
  candidates[1]
}

score_fit <- function(fit_obj, validate, thresholds, probs) {
  if (is.null(fit_obj) || is.null(fit_obj$yp)) {
    return(NULL)
  }

  list(
    quant = QuantScore(fit_obj$yp, probs, validate),
    brier = BrierScore(fit_obj$yp, thresholds, validate)
  )
}

classify_quantile_band <- function(quantiles) {
  bands <- rep("other", length(quantiles))
  bands[quantiles >= 0.90 & quantiles <= 0.95] <- "bulk"
  bands[quantiles >= 0.98] <- "tail"
  bands
}

safe_mean <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

empty_availability_df <- function() {
  data.frame(
    setting = integer(),
    dataset = integer(),
    method_id = integer(),
    method_key = character(),
    family = character(),
    label = character(),
    prop_k = integer(),
    baseline_file = character(),
    prop_file = character(),
    baseline_available = logical(),
    prop_available = logical(),
    stringsAsFactors = FALSE
  )
}

empty_paired_df <- function() {
  data.frame(
    setting = integer(),
    dataset = integer(),
    method_id = integer(),
    method_key = character(),
    family = character(),
    baseline_label = character(),
    proposal_label = character(),
    prop_k = integer(),
    quantile = numeric(),
    band = character(),
    baseline_brier = numeric(),
    prop_brier = numeric(),
    delta_brier = numeric(),
    improve_brier_pct = numeric(),
    baseline_quant = numeric(),
    prop_quant = numeric(),
    delta_quant = numeric(),
    improve_quant_pct = numeric(),
    stringsAsFactors = FALSE
  )
}

empty_summary_df <- function(include_quantile = TRUE) {
  base <- data.frame(
    setting = integer(),
    method_id = integer(),
    method_key = character(),
    family = character(),
    baseline_label = character(),
    proposal_label = character(),
    prop_k = integer(),
    band = character(),
    datasets_compared = integer(),
    baseline_brier = numeric(),
    prop_brier = numeric(),
    delta_brier = numeric(),
    improve_brier_pct = numeric(),
    baseline_quant = numeric(),
    prop_quant = numeric(),
    delta_quant = numeric(),
    improve_quant_pct = numeric(),
    stringsAsFactors = FALSE
  )

  if (include_quantile) {
    base$quantile <- numeric()
    base <- base[, c(
      "setting", "method_id", "method_key", "family",
      "baseline_label", "proposal_label", "prop_k",
      "quantile", "band", "datasets_compared",
      "baseline_brier", "prop_brier", "delta_brier", "improve_brier_pct",
      "baseline_quant", "prop_quant", "delta_quant", "improve_quant_pct"
    )]
  }

  base
}

summarise_metric <- function(df, group_cols) {
  split_idx <- interaction(df[, group_cols], drop = TRUE, lex.order = TRUE)
  pieces <- split(df, split_idx)

  out <- lapply(pieces, function(piece) {
    base <- piece[1, group_cols, drop = FALSE]
    base$datasets_compared <- sum(complete.cases(piece[, c("delta_brier", "delta_quant")]))
    base$baseline_brier <- safe_mean(piece$baseline_brier)
    base$prop_brier <- safe_mean(piece$prop_brier)
    base$delta_brier <- safe_mean(piece$delta_brier)
    base$improve_brier_pct <- safe_mean(piece$improve_brier_pct)
    base$baseline_quant <- safe_mean(piece$baseline_quant)
    base$prop_quant <- safe_mean(piece$prop_quant)
    base$delta_quant <- safe_mean(piece$delta_quant)
    base$improve_quant_pct <- safe_mean(piece$improve_quant_pct)
    base
  })

  do.call(rbind, out)
}

output_root <- "output"
baseline_results_dir <- trimws(Sys.getenv("SIMSTUDY_PROP_BASELINE_RESULTS_DIR", unset = "../simstudy/results"))
prop_results_dir <- trimws(Sys.getenv("SIMSTUDY_PROP_FITS_DIR", unset = ""))
if (!nzchar(prop_results_dir)) {
  prop_results_dir <- trimws(Sys.getenv("SIMSTUDY_PROP_RESULTS_DIR", unset = "fits"))
}
tables_dir <- trimws(Sys.getenv("SIMSTUDY_PROP_COMPARE_OUTPUT_DIR", unset = file.path(output_root, "tables")))
results_dir <- trimws(Sys.getenv("SIMSTUDY_PROP_COMPARE_RESULTS_DIR", unset = file.path(output_root, "results")))
if (!nzchar(prop_results_dir)) {
  prop_results_dir <- "fits"
}
if (!nzchar(tables_dir)) {
  tables_dir <- file.path(output_root, "tables")
}
if (!nzchar(results_dir)) {
  results_dir <- file.path(output_root, "results")
}
for (dir_path in c(tables_dir, results_dir)) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(dir_path)) {
    stop(sprintf("Unable to create output directory '%s'", dir_path), call. = FALSE)
  }
}

setting_ids <- parse_env_indices("SIMSTUDY_PROP_RESULTS_SETTINGS", seq_len(dim(y)[4]), dim(y)[4])
dataset_ids <- parse_env_indices("SIMSTUDY_PROP_RESULTS_DATASETS", seq_len(dim(y)[3]), dim(y)[3])
method_ids <- parse_prop_methods_spec(Sys.getenv("SIMSTUDY_PROP_RESULTS_METHODS", unset = "1:5"))
prop_k_values <- parse_prop_k_spec(Sys.getenv("SIMSTUDY_PROP_RESULTS_K", unset = "10"))

method_catalog <- get_prop_method_catalog()
method_catalog <- method_catalog[method_catalog$method_id %in% method_ids, , drop = FALSE]
method_catalog <- method_catalog[order(method_catalog$method_id), , drop = FALSE]

obs <- c(rep(TRUE, nrow(y) - ntest), rep(FALSE, ntest))
probs <- c(0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995)
bands <- classify_quantile_band(probs)

paired_rows <- list()
available_rows <- list()
missing_rows <- list()
paired_id <- 1L
available_id <- 1L
missing_id <- 1L

for (setting_id in setting_ids) {
  cat("start setting", setting_id, "\n")

  for (dataset_id in dataset_ids) {
    validate <- y[!obs, , dataset_id, setting_id]
    thresholds <- quantile(y[, , dataset_id, setting_id], probs = probs, na.rm = TRUE)

    for (ii in seq_len(nrow(method_catalog))) {
      spec <- method_catalog[ii, , drop = FALSE]
      baseline_file <- build_baseline_result_file(
        results_dir = baseline_results_dir,
        setting_id = setting_id,
        method_id = spec$method_id[1],
        dataset_id = dataset_id
      )
      baseline_fit <- load_fit_object(baseline_file)
      baseline_scores <- score_fit(baseline_fit, validate, thresholds, probs)

      for (prop_k in prop_k_values) {
        prop_file <- resolve_prop_result_file(
          results_dir = prop_results_dir,
          setting_id = setting_id,
          method_id = spec$method_id[1],
          dataset_id = dataset_id,
          prop_k = prop_k
        )
        prop_fit <- load_fit_object(prop_file)
        prop_scores <- score_fit(prop_fit, validate, thresholds, probs)

        availability <- data.frame(
          setting = setting_id,
          dataset = dataset_id,
          method_id = spec$method_id[1],
          method_key = spec$method_key[1],
          family = spec$family[1],
          label = spec$label[1],
          prop_k = as.integer(prop_k),
          baseline_file = baseline_file,
          prop_file = prop_file,
          baseline_available = !is.null(baseline_scores),
          prop_available = !is.null(prop_scores),
          stringsAsFactors = FALSE
        )

        available_rows[[available_id]] <- availability
        available_id <- available_id + 1L

        if (is.null(baseline_scores) || is.null(prop_scores)) {
          missing_rows[[missing_id]] <- availability
          missing_id <- missing_id + 1L
          next
        }

        for (qq in seq_along(probs)) {
          paired_rows[[paired_id]] <- data.frame(
            setting = setting_id,
            dataset = dataset_id,
            method_id = spec$method_id[1],
            method_key = spec$method_key[1],
            family = spec$family[1],
            baseline_label = spec$label[1],
            proposal_label = sprintf("%s + prop(k=%d)", spec$label[1], prop_k),
            prop_k = as.integer(prop_k),
            quantile = probs[qq],
            band = bands[qq],
            baseline_brier = baseline_scores$brier[qq],
            prop_brier = prop_scores$brier[qq],
            delta_brier = prop_scores$brier[qq] - baseline_scores$brier[qq],
            improve_brier_pct = if (is.finite(baseline_scores$brier[qq]) && baseline_scores$brier[qq] != 0) {
              100 * (baseline_scores$brier[qq] - prop_scores$brier[qq]) / baseline_scores$brier[qq]
            } else {
              NA_real_
            },
            baseline_quant = baseline_scores$quant[qq],
            prop_quant = prop_scores$quant[qq],
            delta_quant = prop_scores$quant[qq] - baseline_scores$quant[qq],
            improve_quant_pct = if (is.finite(baseline_scores$quant[qq]) && baseline_scores$quant[qq] != 0) {
              100 * (baseline_scores$quant[qq] - prop_scores$quant[qq]) / baseline_scores$quant[qq]
            } else {
              NA_real_
            },
            stringsAsFactors = FALSE
          )
          paired_id <- paired_id + 1L
        }
      }
    }
  }
}

comparison_prop_paired <- if (length(paired_rows) > 0L) do.call(rbind, paired_rows) else empty_paired_df()
comparison_prop_available <- if (length(available_rows) > 0L) do.call(rbind, available_rows) else empty_availability_df()
comparison_prop_missing <- if (length(missing_rows) > 0L) do.call(rbind, missing_rows) else empty_availability_df()

if (nrow(comparison_prop_paired) > 0L) {
  comparison_prop_summary <- summarise_metric(
    comparison_prop_paired,
    group_cols = c(
      "setting", "method_id", "method_key", "family",
      "baseline_label", "proposal_label", "prop_k", "quantile", "band"
    )
  )

  band_only <- comparison_prop_paired[comparison_prop_paired$band %in% c("bulk", "tail"), , drop = FALSE]
  comparison_prop_band_summary <- if (nrow(band_only) > 0L) {
    summarise_metric(
      band_only,
      group_cols = c(
        "setting", "method_id", "method_key", "family",
        "baseline_label", "proposal_label", "prop_k", "band"
      )
    )
  } else {
    empty_summary_df(include_quantile = FALSE)
  }
} else {
  comparison_prop_summary <- empty_summary_df(include_quantile = TRUE)
  comparison_prop_band_summary <- empty_summary_df(include_quantile = FALSE)
}

write.csv(
  comparison_prop_paired,
  file.path(tables_dir, "comparison_prop_paired.csv"),
  row.names = FALSE
)
write.csv(
  comparison_prop_summary,
  file.path(tables_dir, "comparison_prop_summary.csv"),
  row.names = FALSE
)
write.csv(
  comparison_prop_band_summary,
  file.path(tables_dir, "comparison_prop_band_summary.csv"),
  row.names = FALSE
)
write.csv(
  comparison_prop_available,
  file.path(tables_dir, "comparison_prop_available.csv"),
  row.names = FALSE
)
write.csv(
  comparison_prop_missing,
  file.path(tables_dir, "comparison_prop_missing.csv"),
  row.names = FALSE
)

save(
  comparison_prop_paired,
  comparison_prop_summary,
  comparison_prop_band_summary,
  comparison_prop_available,
  comparison_prop_missing,
  file = file.path(results_dir, "results-prop.RData")
)

cat(sprintf(
  "results-prop complete: settings=%s datasets=%s methods=%s prop_k=%s tables_dir=%s results_dir=%s\n",
  paste(setting_ids, collapse = ","),
  paste(dataset_ids, collapse = ","),
  paste(method_ids, collapse = ","),
  paste(prop_k_values, collapse = ","),
  tables_dir,
  results_dir
))
