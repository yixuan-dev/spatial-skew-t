#########################################################################
# tables.R - Stage 2 of the time-block forecast post-fit pipeline.
#
# Reads output/results/scores<setting><suffix>.RData (from scores.R) and
# emits the lead-time evaluation *tables* of Section 4.3 of
# ar2_rethink.tex:
#
#   - lead-time curve  S_bar(h)  per method, with a standard error at
#     each lead estimated from the dispersion of the B per-block means
#     (Definition 9);
#   - the relative curve AR(2) / i.i.d. and the crossing lead -- the
#     headline quantity: the AR(2) curve should sit below the baseline
#     at short leads and converge near the memory horizon h*;
#   - the joint-structure summary (energy + variogram score per method).
#
# Plotting moved to Stage 3 (plots.R), which consumes the simresults
# artifact written here. This script produces tables only and has no
# graphics dependency.
#
# Outputs (suffix "" for simdata.RData):
#   output/tables/lead_curve<setting><suffix>.csv
#   output/tables/lead_curve_rel<setting><suffix>.csv
#   output/tables/joint_summary<setting><suffix>.csv
#   output/results/simresults<setting><suffix>.RData
#
# Usage:
#   Rscript tables.R --setting=<id> [--data=<path>]
#########################################################################

rm(list = ls())

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0L) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]),
    winslash = "/", mustWork = FALSE)
  if (dir.exists(dirname(script_path))) setwd(dirname(script_path))
}

source("./time_block_helpers.R")

cli_args <- commandArgs(trailingOnly = TRUE)
prior <- tbf_take_prior_flag(cli_args, c("data", "setting"))
parsed <- prior$parsed
flags <- parsed$values
prior_tag <- prior$prior_tag

if (is.null(flags$setting) || !nzchar(flags$setting)) {
  stop("tables.R: --setting=<id> is required.", call. = FALSE)
}
setting_id <- as.integer(parse_index_expr(flags$setting, "setting"))

data_suffix <- if (!is.null(flags$data) && nzchar(flags$data)) {
  derive_data_suffix(flags$data)
} else {
  ""
}

results_dir <- "output/results"
tables_dir <- "output/tables"
for (d in c(results_dir, tables_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

scores_file <- tbf_score_cache_file(setting_id, data_suffix, prior_tag,
  dir = results_dir)
if (!file.exists(scores_file)) {
  stop(sprintf("Score cache not found: %s\n  Run scores.R --setting=%d first.",
    scores_file, setting_id), call. = FALSE)
}
load(scores_file)
# Provides: crps.lead, brier.lead (full-series thresholds since
#           2026-07-31; brier_threshold_basis records the rule),
#           brier.lead.blockq, pexceed.mean, brier.thresholds,
#           exceed.rate.lead, energy.score, vario.score, elapsed_sec,
#           probs, datasets, methods, setting, block_H, block_seams, ...

catalog <- get_tbf_method_catalog()
H <- block_H

cat(sprintf("tables: setting=%d cache=%s\n", setting_id, scores_file))

# ---- lead-time curve S_bar(h) ----------------------------------------
# crps.lead is [lead, dataset, method, block]. Per Definition 9, average
# over datasets to get a per-block mean, then take the mean and the
# standard error across the B blocks at each lead.
lead_curve <- function(arr, score_name) {
  rows <- list()
  for (mi in seq_along(methods)) {
    for (h in seq_len(H)) {
      block_means <- apply(arr[h, , mi, , drop = FALSE], 4, mean, na.rm = TRUE)
      block_means <- block_means[is.finite(block_means)]
      mu <- mean(block_means)
      se <- if (length(block_means) > 1L) {
        sd(block_means) / sqrt(length(block_means))
      } else {
        NA_real_
      }
      rows[[length(rows) + 1L]] <- data.frame(
        score = score_name, method = methods[mi], lead = h,
        mean = mu, se = se
      )
    }
  }
  do.call(rbind, rows)
}

crps_curve <- lead_curve(crps.lead, "crps")
# average the Brier lead curve over the threshold-quantile grid
brier_mean_lead <- apply(brier.lead, c(1, 3, 4, 5), mean, na.rm = TRUE)
brier_curve <- lead_curve(brier_mean_lead, "brier")
curve_table <- rbind(crps_curve, brier_curve)
write.csv(curve_table,
  file.path(tables_dir, sprintf("lead_curve%d%s.csv", setting_id, data_suffix)),
  row.names = FALSE)

# ---- relative curve + crossing lead ----------------------------------
baseline <- min(methods)            # method 1 = i.i.d. baseline
rel_rows <- list()
crossing_rows <- list()
for (score_name in c("crps", "brier")) {
  ct <- curve_table[curve_table$score == score_name, ]
  base_mean <- ct$mean[ct$method == baseline][order(ct$lead[ct$method == baseline])]
  for (m in setdiff(methods, baseline)) {
    m_mean <- ct$mean[ct$method == m][order(ct$lead[ct$method == m])]
    rel <- m_mean / base_mean
    rel_rows[[length(rel_rows) + 1L]] <- data.frame(
      score = score_name, method = m, lead = seq_len(H),
      rel_mean = rel
    )
    # crossing lead: first lead at which AR(2) stops beating the baseline
    worse <- which(rel >= 1)
    crossing <- if (length(worse) == 0L) NA_integer_ else min(worse)
    crossing_rows[[length(crossing_rows) + 1L]] <- data.frame(
      score = score_name, method = m,
      crossing_lead = crossing,
      mean_rel_short = mean(rel[seq_len(min(5L, H))], na.rm = TRUE)
    )
  }
}
rel_table <- do.call(rbind, rel_rows)
crossing_table <- do.call(rbind, crossing_rows)
write.csv(rel_table,
  file.path(tables_dir, sprintf("lead_curve_rel%d%s.csv", setting_id, data_suffix)),
  row.names = FALSE)

# ---- joint-structure summary -----------------------------------------
setting_cat <- get_tbf_setting_catalog()
phi_row <- setting_cat[setting_cat$setting_id == setting_id, , drop = FALSE]
# ARFIMA settings have no AR(2) memory horizon: the ACF decays
# hyperbolically, so h_star is NA and a Yule-Walker AR(2) projection stands
# in as a comparable proxy. tbf_memory_horizon() dispatches on the family.
hs <- tbf_memory_horizon(setting_id)
h_star <- hs$h_star
if (!is.na(hs$note)) cat("  ", hs$note, "\n", sep = "")

joint_rows <- list()
for (mi in seq_along(methods)) {
  joint_rows[[mi]] <- data.frame(
    method = methods[mi],
    label = catalog$label[catalog$method_id == methods[mi]],
    energy_mean = mean(energy.score[, mi, ], na.rm = TRUE),
    vario_mean = mean(vario.score[, mi, ], na.rm = TRUE),
    elapsed_sec_mean = mean(elapsed_sec[, mi], na.rm = TRUE)
  )
}
joint_table <- do.call(rbind, joint_rows)
joint_table$h_star <- h_star
joint_table$h_star_basis <- hs$basis
joint_table$h_star_proxy <- hs$proxy
joint_table$family <- phi_row$family[1]
joint_table$d <- phi_row$d[1]
joint_table$hurst <- phi_row$hurst[1]
joint_table <- merge(joint_table, crossing_table[crossing_table$score == "crps",
  c("method", "crossing_lead")], by = "method", all.x = TRUE)
write.csv(joint_table,
  file.path(tables_dir, sprintf("joint_summary%d%s.csv", setting_id, data_suffix)),
  row.names = FALSE)

# ---- aggregated artifact (consumed by plots.R, Stage 3) --------------
simresults_file <- tbf_score_cache_file(setting_id, data_suffix, prior_tag,
  dir = results_dir, prefix = "simresults")
save(curve_table, rel_table, crossing_table, joint_table,
  h_star, hs, probs, methods, datasets, setting, block_H, block_seams,
  data_suffix,
  file = simresults_file)

cat("\nWrote tables to ", tables_dir, "\n", sep = "")
cat("Saved analysis objects to ", simresults_file, "\n", sep = "")
cat("Next: Rscript plots.R --setting=", setting_id,
  " for the lead-time figures\n", sep = "")
