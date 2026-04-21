rm(list = ls())

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = FALSE))
} else {
  getwd()
}

source(file.path(script_dir, "00_bootstrap.R"))
source(file.path(script_dir, "01_score_engine.R"))
source(file.path(script_dir, "02_comparison_tables.R"))
source(file.path(script_dir, "03_plot_results.R"))

set_working_dir_to_script()
ctx <- load_us_all_context()
ensure_us_all_output_dirs()
ensure_plot_packages()

results_bundle_path <- us_all_output_path("us-all-results-proposed.RData", subdir = "results")
if (!file.exists(results_bundle_path)) {
  stop(
    "Missing proposed results bundle: ",
    results_bundle_path,
    ". Run refactor_results/12_run_proposed_results.R first."
  )
}

bundle_env <- new.env(parent = emptyenv())
load(results_bundle_path, envir = bundle_env)

required_objects <- c(
  "available_settings",
  "comparison_brier_split",
  "comparison_paired_same_basis",
  "summary_obj"
)
missing_objects <- required_objects[!vapply(
  required_objects,
  function(nm) exists(nm, envir = bundle_env, inherits = FALSE),
  logical(1)
)]
if (length(missing_objects) > 0) {
  stop(
    "Proposed results bundle is missing required objects: ",
    paste(missing_objects, collapse = ", ")
  )
}

focus_settings <- build_mrts_focus_settings(
  settings = ctx$settings,
  comparison_paired_same_basis = get("comparison_paired_same_basis", envir = bundle_env),
  available_settings = get("available_settings", envir = bundle_env)
)
if (nrow(focus_settings) == 0) {
  stop("No MRTS comparison rows were found in us-all-results-proposed.RData.")
}

result_dirs <- resolve_plot_result_dirs()
plot_summary_draws <- parse_plot_summary_draws()
extreme_probs <- default_extreme_plot_probs()

mean_metrics <- compute_focus_mean_prediction_metrics(
  focus_settings = focus_settings,
  result_dirs = result_dirs,
  Y = ctx$Y,
  cv_lst = ctx$cv_lst,
  trans_setting_ids = 2L,
  summary_draws = plot_summary_draws,
  enforce_contract = TRUE
)
if (nrow(mean_metrics) == 0) {
  stop("Unable to compute mean-level diagnostics for the MRTS focus settings.")
}

mean_metrics_path <- us_all_output_path("mrts_mean_diagnostics.csv", subdir = "tables")
write.csv(mean_metrics, mean_metrics_path, row.names = FALSE)

plot_specs <- list(
  list(
    name = "mrts_extreme_delta_profiles.png",
    fn = function(path) {
      plot_mrts_extreme_delta_profiles(
        comparison_paired_same_basis = get("comparison_paired_same_basis", envir = bundle_env),
        focus_settings = focus_settings,
        output_path = path,
        target_probs = extreme_probs
      )
    }
  ),
  list(
    name = "mrts_extreme_split_brier.png",
    fn = function(path) {
      plot_mrts_extreme_split_brier(
        comparison_brier_split = get("comparison_brier_split", envir = bundle_env),
        focus_settings = mean_metrics,
        output_path = path,
        target_probs = extreme_probs[extreme_probs >= 0.95]
      )
    }
  ),
  list(
    name = "mrts_mean_error_gap.png",
    fn = function(path) {
      plot_mrts_mean_error_gap(
        mean_metrics = mean_metrics,
        output_path = path
      )
    }
  ),
  list(
    name = "mrts_mean_bias.png",
    fn = function(path) {
      plot_mrts_mean_bias(
        mean_metrics = mean_metrics,
        output_path = path
      )
    }
  ),
  list(
    name = "mrts_tail_mean_tradeoff.png",
    fn = function(path) {
      plot_mrts_tail_mean_tradeoff(
        comparison_paired_same_basis = get("comparison_paired_same_basis", envir = bundle_env),
        mean_metrics = mean_metrics,
        output_path = path,
        target_prob = 0.99,
        metric = "brier"
      )
    }
  )
)

written_plots <- character(0)
for (spec in plot_specs) {
  plot_path <- us_all_output_path(spec$name, subdir = "plots")
  ok <- spec$fn(plot_path)
  if (isTRUE(ok)) {
    written_plots <- c(written_plots, plot_path)
  }
}

cat("Result directories searched (priority order):\n")
for (d in result_dirs) {
  cat("-", d, "\n")
}

cat("Plot summary draws used:", plot_summary_draws, "\n")
cat("Focus settings:\n")
for (i in seq_len(nrow(mean_metrics))) {
  cat(
    "-",
    mean_metrics$model_label[i],
    "| basis:", mean_metrics$basis_label[i],
    "\n"
  )
}

cat("Outputs written:\n")
cat("- ", mean_metrics_path, "\n", sep = "")
for (path in written_plots) {
  cat("- ", path, "\n", sep = "")
}
