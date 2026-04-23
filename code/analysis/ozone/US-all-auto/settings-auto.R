# settings-auto.R
#
# Reads the existing settings.csv and performance tables from the sibling US-all
# experiment, then produces settings-auto.csv — an evidence-based, focused
# settings grid for the US-all-auto exceedance auto-selection experiment.
#
# Design principles:
#   - No 2-fold CV is needed for the auto-select experiment.  The grid here only
#     lists settings whose result files (us-all-<setting>.RData) already exist in
#     US-all.  The auto-selection script scores them on a held-out validation split
#     and picks the best one — no re-training, no CV averaging.
#
#   - Evidence from comparison_full_table.csv drives tier assignment.
#     Composite tail Brier = mean of rel_to_gaussian at q ∈ {0.95, 0.98, 0.99, 0.995}.
#     Lower is better.  Tier boundary is set so that Tier 1 is approximately the
#     top third of scored settings and Tier 2 covers the rest.
#
#   - Settings with MRTS k=800 are excluded (catastrophically bad at all quantiles).
#
#   - Non-TS, non-AR2, non-MRTS baseline settings (1–50) are assigned Tier 3 and
#     excluded from the default search range, but retained in the file for reference.
#     Exception: Setting 1 (Gaussian) is always Tier 0 (reference).
#
# Outputs:
#   settings-auto.csv   — columns: setting, method, knots, thresh, CMAQ, TS, ar2,
#                         mrts, lane, tier, composite_tail_brier, best_tail_q,
#                         n_top5_appearances
#
# How to run:
#   cd D:\Github\spatial-skew-t\code\analysis\ozone\US-all-auto
#   Rscript settings-auto.R

rm(list = ls())

# ---- Locate paths ----

.this_dir     <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg) > 0)
    dirname(normalizePath(sub("^--file=", "", arg[1]), winslash = "/", mustWork = FALSE))
  else
    normalizePath(".", winslash = "/", mustWork = FALSE)
})
setwd(.this_dir)

.us_all_root  <- normalizePath(file.path(.this_dir, "../US-all"), winslash = "/", mustWork = FALSE)
.tables_dir   <- file.path(.us_all_root, "output", "us-all", "tables")
.results_dir  <- file.path(.us_all_root, "results")

cat("=== settings-auto.R ===\n")
cat("US-all root :", .us_all_root, "\n")
cat("Tables dir  :", .tables_dir, "\n\n")

# ---- Load inputs ----

settings_path <- file.path(.us_all_root, "settings.csv")
if (!file.exists(settings_path)) stop("Not found: ", settings_path)
settings_raw <- read.csv(settings_path, stringsAsFactors = FALSE)
cat("settings.csv: ", nrow(settings_raw), "rows\n")

full_table_path <- file.path(.tables_dir, "comparison_full_table.csv")
if (!file.exists(full_table_path)) stop("Not found: ", full_table_path)
full_table <- read.csv(full_table_path, stringsAsFactors = FALSE)
cat("comparison_full_table.csv: ", nrow(full_table), "rows\n\n")

# ---- Parse settings ----

settings <- settings_raw
settings$setting_num <- suppressWarnings(as.integer(settings$setting))
settings <- settings[!is.na(settings$setting_num), ]   # drop rows like "5a", "6a"

normalize_yn <- function(x) tolower(trimws(as.character(x)))
parse_int_or_na <- function(x) {
  v <- suppressWarnings(as.integer(trimws(as.character(x))))
  v[!is.finite(v)] <- NA_integer_
  v
}

settings$TS    <- normalize_yn(settings$TS)
settings$ar2   <- normalize_yn(settings$ar2)
settings$CMAQ  <- normalize_yn(settings$CMAQ)
settings$mrts_k <- parse_int_or_na(settings$mrts)
settings$thresh <- suppressWarnings(as.numeric(settings$thresh))
settings$knots  <- suppressWarnings(as.integer(settings$knots))

# ---- Identify available result files ----

settings$has_result <- file.exists(
  file.path(.results_dir, sprintf("us-all-%d.RData", settings$setting_num))
)
cat(sprintf("Result files found: %d / %d settings\n\n",
    sum(settings$has_result), nrow(settings)))

# ---- Compute composite tail Brier score from the full table ----
# composite_tail_brier = mean rel_to_gaussian at Brier q ∈ {0.95, 0.98, 0.99, 0.995}

tail_qs     <- c(0.95, 0.98, 0.99, 0.995)
brier_tail  <- full_table[
  full_table$metric == "brier" &
  full_table$quantile %in% tail_qs &
  is.finite(full_table$rel_to_gaussian), ]

# Per-setting summary
score_summary <- do.call(rbind, lapply(split(brier_tail, brier_tail$setting), function(sub) {
  rels <- sub$rel_to_gaussian
  qs   <- sub$quantile
  data.frame(
    setting              = sub$setting[1],
    composite_tail_brier = mean(rels, na.rm = TRUE),
    best_tail_q          = qs[which.min(rels)],
    n_top5_appearances   = 0L,   # filled in below
    stringsAsFactors     = FALSE
  )
}))
score_summary$setting <- as.integer(score_summary$setting)

# Count top-5 appearances at each tail quantile
for (q in tail_qs) {
  sub_q   <- brier_tail[brier_tail$quantile == q, ]
  sub_q   <- sub_q[order(sub_q$rel_to_gaussian), ]
  top5_ids <- as.integer(head(sub_q$setting, 5))
  hit      <- score_summary$setting %in% top5_ids
  score_summary$n_top5_appearances[hit] <- score_summary$n_top5_appearances[hit] + 1L
}

# ---- Assign lanes ----
# Each setting belongs to exactly one lane.

assign_lane <- function(s) {
  dplyr_free_case <- function(row) {
    if (row$setting_num == 1)          return("gaussian_reference")
    if (!is.na(row$mrts_k))            return(sprintf("mrts_k%d", row$mrts_k))
    if (row$ar2 == "yes")              return("ar2")
    if (row$TS  == "yes")              return("ts_baseline")
    return("baseline")
  }
  vapply(seq_len(nrow(s)), function(i) dplyr_free_case(s[i, ]), character(1))
}

settings$lane <- assign_lane(settings)

# ---- Assign tiers ----
#
# Tier 0  Gaussian reference (setting 1)
# Tier 1  Core candidates:
#           top-20 by composite_tail_brier among scored settings AND
#           lane ∈ {ts_baseline, ar2, mrts_k5, mrts_k10}
# Tier 2  Extended candidates:
#           remaining ts_baseline / ar2 / mrts_k5 settings not in Tier 1, plus
#           mrts_k15 (mixed evidence, include for reference)
# Tier 3  Baseline (non-TS) settings: informative but not primary candidates
# excl    Excluded: mrts_k800 (catastrophically bad)

settings <- merge(settings, score_summary[, c("setting", "composite_tail_brier",
                                               "best_tail_q", "n_top5_appearances")],
                  by.x = "setting_num", by.y = "setting", all.x = TRUE)

settings$tier <- NA_integer_

# Tier 0: Gaussian reference
settings$tier[settings$setting_num == 1] <- 0L

# Exclude: MRTS k=800
settings$tier[!is.na(settings$mrts_k) & settings$mrts_k == 800] <- NA_integer_

# Identify top-20 scored settings (excluding Gaussian and mrts_k800)
candidates <- settings[
  is.na(settings$tier) &
  settings$has_result &
  settings$lane %in% c("ts_baseline", "ar2", "mrts_k5", "mrts_k10", "mrts_k15", "baseline") &
  is.finite(settings$composite_tail_brier), ]

candidates <- candidates[order(candidates$composite_tail_brier), ]
tier1_top_n <- 20L

tier1_ids <- candidates$setting_num[seq_len(min(tier1_top_n, nrow(candidates)))]

# Always include MRTS k=10 in Tier 1 (best MRTS variant)
mrts10_ids <- settings$setting_num[
  !is.na(settings$mrts_k) & settings$mrts_k == 10 & settings$has_result]
tier1_ids  <- sort(unique(c(tier1_ids, mrts10_ids)))

# Always include t + thresh=75 at key knots in Tier 1
# (dominant pattern at q=0.99-0.995; knots 7, 9, 15 from evidence)
key_t75 <- settings$setting_num[
  settings$method == "t" &
  settings$thresh == 75 &
  settings$knots %in% c(7, 9, 15) &
  settings$lane %in% c("ts_baseline", "ar2") &
  settings$has_result]
tier1_ids <- sort(unique(c(tier1_ids, key_t75)))

# Always include skew-t + thresh=50 at key knots in Tier 1
# (dominant pattern at q=0.95-0.98; knots 5, 7, 9, 10 from evidence)
key_st50 <- settings$setting_num[
  settings$method == "skew-t" &
  settings$thresh == 50 &
  settings$knots %in% c(5, 7, 9, 10) &
  settings$lane %in% c("ts_baseline", "ar2", "baseline") &
  settings$has_result]
tier1_ids <- sort(unique(c(tier1_ids, key_st50)))

settings$tier[settings$setting_num %in% tier1_ids & is.na(settings$tier)] <- 1L

# Tier 2: remaining TS / AR2 / MRTS k=5 / MRTS k=15 settings with result files
tier2_lanes <- c("ts_baseline", "ar2", "mrts_k5", "mrts_k15")
settings$tier[
  is.na(settings$tier) &
  settings$lane %in% tier2_lanes &
  settings$has_result] <- 2L

# Tier 3: baseline (non-TS) settings with result files
settings$tier[
  is.na(settings$tier) &
  settings$lane == "baseline" &
  settings$has_result] <- 3L

# Remaining with result files but no score (e.g. no result in full_table)
settings$tier[is.na(settings$tier) & settings$has_result] <- 3L

# ---- Build settings-auto output ----

out_cols <- c("setting_num", "method", "knots", "thresh", "CMAQ", "TS", "ar2",
              "mrts_k", "lane", "tier",
              "composite_tail_brier", "best_tail_q", "n_top5_appearances",
              "has_result")

settings_auto <- settings[!is.na(settings$tier), out_cols]
names(settings_auto)[names(settings_auto) == "setting_num"] <- "setting"
names(settings_auto)[names(settings_auto) == "mrts_k"]      <- "mrts"

settings_auto <- settings_auto[order(settings_auto$tier, settings_auto$setting), ]
rownames(settings_auto) <- NULL

# ---- Print summary ----

cat("=== settings-auto grid summary ===\n\n")
cat(sprintf("%-6s  %-20s  %s\n", "Tier", "Lane", "N settings"))
cat(strrep("-", 50), "\n")
for (tier_i in sort(unique(settings_auto$tier))) {
  sub <- settings_auto[settings_auto$tier == tier_i, ]
  lanes_str <- paste(sort(unique(sub$lane)), collapse = ", ")
  cat(sprintf("%-6s  %-20s  %d\n", tier_i, lanes_str, nrow(sub)))
}
cat("\n")

cat("--- Tier 1 settings (core candidates) ---\n")
t1 <- settings_auto[settings_auto$tier == 1, ]
t1_display <- t1[, c("setting","method","knots","thresh","lane",
                      "composite_tail_brier","best_tail_q","n_top5_appearances")]
t1_display <- t1_display[order(t1_display$composite_tail_brier, na.last = TRUE), ]
print(t1_display, row.names = FALSE)
cat("\n")

cat("Evidence interpretation (from comparison_full_table.csv):\n")
cat("  thresh=50 + skew-t  => best at q=0.95-0.98 Brier\n")
cat("  thresh=75 + t model => best at q=0.99-0.995 Brier\n")
cat("  knots 5-10          => composite sweet spot\n")
cat("  MRTS k=10           => best MRTS variant; k=800 excluded (harmful)\n")
cat("  AR2 vs TS baseline  => AR2 often +1-2% at extremes; include both\n\n")

# ---- Write output ----

out_path <- file.path(.this_dir, "settings-auto.csv")
write.csv(settings_auto, file = out_path, row.names = FALSE)
cat("Written:", out_path, "\n")
cat("Total settings:", nrow(settings_auto), "\n")
cat("  Tier 0:", sum(settings_auto$tier == 0), " (Gaussian reference)\n")
cat("  Tier 1:", sum(settings_auto$tier == 1), " (core exceedance candidates)\n")
cat("  Tier 2:", sum(settings_auto$tier == 2), " (extended coverage)\n")
cat("  Tier 3:", sum(settings_auto$tier == 3), " (baseline, not recommended for default run)\n")
cat("\nSuggested US_ALL_AUTOSELECT_SETTINGS for Tier 1:\n")
t1_ids <- sort(settings_auto$setting[settings_auto$tier == 1])
cat(" ", paste(
  vapply(split(t1_ids, cumsum(c(1, diff(t1_ids)) != 1)), function(g) {
    if (length(g) == 1) as.character(g)
    else sprintf("%d:%d", g[1], g[length(g)])
  }, character(1)),
  collapse = ","
), "\n")
cat("\nSuggested US_ALL_AUTOSELECT_SETTINGS for Tier 1+2:\n")
t12_ids <- sort(settings_auto$setting[settings_auto$tier %in% c(1, 2)])
cat(" ", paste(
  vapply(split(t12_ids, cumsum(c(1, diff(t12_ids)) != 1)), function(g) {
    if (length(g) == 1) as.character(g)
    else sprintf("%d:%d", g[1], g[length(g)])
  }, character(1)),
  collapse = ","
), "\n")
