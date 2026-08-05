#########################################################################
# expA_threeway.R -- Experiment A (5 blocks) three-way comparison tables
# and paired one-sided tests, for BOTH CRPS and Brier q95.
#
# Reads the merged score cache written by scores.R / merge_score_caches.R
# and emits the artifacts the time-block score report consumes:
#
#   output/tables/expA_lead_<setting>.csv   per-lead means for each method
#                                           plus the three gaps
#   output/tables/expA_tests_<setting>.csv  paired one-sided t and Wilcoxon
#                                           for 3 contrasts x 2 scores x
#                                           2 windows x 2 pairing units
#   output/tables/expA_diag_<setting>.csv   per (dataset, block, method)
#                                           diagnostics, from the cache
#
# Pairing units. Primary is the DATASET (n = 10): the five blocks of one
# dataset are averaged first, which absorbs their dependence and matches
# the unit used by the block-1 Experiment B tables. The cell unit
# (dataset x block, n = 50) is reported alongside as a power sensitivity;
# its replicates are not independent, so its p-values are optimistic.
#
# Usage:
#   Rscript expA_threeway.R --setting=<id> [--data=<path>]
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
prior <- tbf_take_prior_flag(cli_args, c("setting", "data"))
parsed <- prior$parsed
flags <- parsed$values
prior_tag <- prior$prior_tag
if (is.null(flags$setting) || !nzchar(flags$setting)) {
  stop("expA_threeway.R: --setting=<id> is required.", call. = FALSE)
}
setting_id <- as.integer(parse_index_expr(flags$setting, "setting"))
data_suffix <- if (!is.null(flags$data) && nzchar(flags$data)) {
  derive_data_suffix(flags$data)
} else {
  ""
}

tables_dir <- "output/tables"
if (!dir.exists(tables_dir)) dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

scores_file <- tbf_score_cache_file(setting_id, data_suffix, prior_tag)
if (!file.exists(scores_file)) {
  stop(sprintf("score cache not found: %s (run scores.R --prior=%s)",
    scores_file, prior_tag), call. = FALSE)
}
load(scores_file)
# Provides: crps.lead, brier.lead (full-series thresholds since
#           2026-07-31; brier_threshold_basis records the rule),
#           brier.lead.blockq, pexceed.mean, brier.thresholds,
#           exceed.rate.lead, fit.diag, post.summary, hn_prior,
#           es_max_draws, probs, datasets, methods, setting, block_H, ...

catalog <- get_tbf_method_catalog()
mlab <- c("1" = "iid", "2" = "AR2", "4" = "AR1")
H <- block_H
q95 <- match(0.95, probs)
if (is.na(q95)) stop("the probs grid has no 0.95 threshold", call. = FALSE)

cat(sprintf("expA_threeway: setting %d, %d datasets, methods [%s], %d blocks\n",
  setting_id, length(datasets), paste(methods, collapse = ","),
  dim(crps.lead)[4]))
cat(sprintf("  hn_prior = %s, es_max_draws = %s\n",
  format(hn_prior), format(es_max_draws)))

# score arrays as [lead, dataset, method, block]
SCORES <- list(
  CRPS = crps.lead,
  Brier95 = brier.lead[, q95, , , , drop = TRUE]
)

mi <- function(m) match(as.character(m), dimnames(crps.lead)[["method"]])

# ---- per-lead table ----------------------------------------------------
lead_rows <- list()
for (sc in names(SCORES)) {
  arr <- SCORES[[sc]]
  wide <- sapply(methods, function(m) apply(arr[, , mi(m), , drop = FALSE],
    1, mean, na.rm = TRUE))
  colnames(wide) <- mlab[as.character(methods)]
  df <- data.frame(setting = setting_id, score = sc, lead = seq_len(H),
    wide, check.names = FALSE)
  if (all(c("AR1", "iid") %in% colnames(wide))) df$`AR1-iid` <- wide[, "AR1"] - wide[, "iid"]
  if (all(c("AR2", "iid") %in% colnames(wide))) df$`AR2-iid` <- wide[, "AR2"] - wide[, "iid"]
  if (all(c("AR2", "AR1") %in% colnames(wide))) df$`AR2-AR1` <- wide[, "AR2"] - wide[, "AR1"]
  lead_rows[[sc]] <- df
}
lead_tab <- do.call(rbind, lead_rows)
lead_file <- file.path(tables_dir, sprintf("expA_lead_%d%s.csv", setting_id, data_suffix))
write.csv(lead_tab, lead_file, row.names = FALSE)

# ---- paired one-sided tests -------------------------------------------
# H1: the first method of the contrast is BETTER, i.e. its score is lower,
# i.e. the paired difference is negative -- hence alternative = "less".
one_sided <- function(d) {
  d <- d[is.finite(d)]
  if (length(d) < 3L || isTRUE(all.equal(sd(d), 0))) {
    return(c(n = length(d), gap = if (length(d)) mean(d) else NA_real_,
             t_p = NA_real_, w_p = NA_real_))
  }
  tp <- tryCatch(t.test(d, alternative = "less")$p.value, error = function(e) NA_real_)
  wp <- tryCatch(suppressWarnings(
    wilcox.test(d, alternative = "less", exact = FALSE)$p.value),
    error = function(e) NA_real_)
  c(n = length(d), gap = mean(d), t_p = tp, w_p = wp)
}

# window x unit -> the vector of paired differences
paired_vec <- function(arr, mA, mB, leads, unit) {
  a <- arr[leads, , mi(mA), , drop = FALSE]
  b <- arr[leads, , mi(mB), , drop = FALSE]
  # average over the lead window first: the window is the endpoint, not a
  # replicate axis
  a <- apply(a, c(2, 4), mean, na.rm = TRUE)   # [dataset, block]
  b <- apply(b, c(2, 4), mean, na.rm = TRUE)
  d <- a - b
  if (identical(unit, "dataset")) rowMeans(d, na.rm = TRUE) else as.numeric(d)
}

CONTRASTS <- list(
  c("4", "1", "AR1-iid"),
  c("2", "1", "AR2-iid"),
  c("2", "4", "AR2-AR1")
)
WINDOWS <- list("1-5" = 1:5, "1-15" = seq_len(H))

test_rows <- list()
for (ct in CONTRASTS) {
  mA <- as.integer(ct[1]); mB <- as.integer(ct[2])
  if (!all(c(mA, mB) %in% methods)) next
  for (sc in names(SCORES)) {
    for (wn in names(WINDOWS)) {
      for (unit in c("dataset", "cell")) {
        v <- paired_vec(SCORES[[sc]], mA, mB, WINDOWS[[wn]], unit)
        st <- one_sided(v)
        test_rows[[length(test_rows) + 1L]] <- data.frame(
          setting = setting_id, score = sc, contrast = ct[3], window = wn,
          unit = unit, primary = unit == "dataset",
          n = st[["n"]], gap = st[["gap"]],
          t_p = st[["t_p"]], w_p = st[["w_p"]]
        )
      }
    }
  }
}
test_tab <- do.call(rbind, test_rows)
test_file <- file.path(tables_dir, sprintf("expA_tests_%d%s.csv", setting_id, data_suffix))
write.csv(test_tab, test_file, row.names = FALSE)

# ---- diagnostics, flattened from the cache ----------------------------
diag_file <- NA_character_
if (exists("fit.diag") && !all(is.na(fit.diag))) {
  dn <- dimnames(fit.diag)
  grid <- expand.grid(dataset = dn$dataset, method = dn$method, block = dn$block,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  vals <- t(vapply(seq_len(nrow(grid)), function(i) {
    fit.diag[, grid$dataset[i], grid$method[i], grid$block[i]]
  }, numeric(length(dn$stat))))
  colnames(vals) <- dn$stat
  diag_tab <- cbind(setting = setting_id, grid, as.data.frame(vals))
  diag_file <- file.path(tables_dir, sprintf("expA_diag_%d%s.csv", setting_id, data_suffix))
  write.csv(diag_tab, diag_file, row.names = FALSE)
}

# ---- console summary ---------------------------------------------------
cat("\n==== per-lead means ====\n")
print(lead_tab[lead_tab$score == "CRPS", ], row.names = FALSE, digits = 4)
cat("\n==== paired one-sided tests (H1: first method better) ====\n")
show <- test_tab[test_tab$unit == "dataset", ]
cat(sprintf("%-8s %-9s %-6s %4s %10s %9s %9s\n",
  "score", "contrast", "window", "n", "gap", "t p", "Wilcoxon p"))
for (i in seq_len(nrow(show))) {
  r <- show[i, ]
  cat(sprintf("%-8s %-9s %-6s %4d %10.4f %9.4f %9.4f\n",
    r$score, r$contrast, r$window, r$n, r$gap, r$t_p, r$w_p))
}
cat("\n(cell-unit rows, n = 50, are in the CSV as the power sensitivity)\n")

cat("\nwrote ", lead_file, "\n", sep = "")
cat("wrote ", test_file, "\n", sep = "")
if (!is.na(diag_file)) cat("wrote ", diag_file, "\n", sep = "")
