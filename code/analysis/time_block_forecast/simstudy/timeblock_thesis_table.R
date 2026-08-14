#########################################################################
# timeblock_thesis_table.R -- the thesis table for the time-block arm.
#
# Emits the body of tab:timeblock-brier from the Experiment A caches, so
# that no number in the thesis is transcribed by hand. Shape mirrors the
# spatial hold-out table: mean level per method, then the two contrasts
# that carry the argument, each with a two-sided 95% paired-t interval.
#
#   rows    3 generating designs x 2 lead windows, twice: once per
#           co-primary threshold (q90 then q95)
#   columns i.i.d. | AR(1) | AR(2) | AR(2)-i.i.d. [CI] | AR(2)-AR(1) [CI]
#
# Pairing unit is the DATA SET (n = 10): the lead window is averaged
# first, then the five blocks, leaving one number per data set. Blocks
# share training data and are not independent replicates, so they are
# never used as the unit. Tests are TWO-sided here, matching the thesis
# convention of sec:sim-uq; expA_breakdown.R reports the one-sided form
# for the standalone score reports.
#
# Differences are printed in units of 1e-3, as in the spatial table and
# in the lead-difference figure.
#
# Usage:
#   Rscript timeblock_thesis_table.R [--settings=(4,5,7)] [--prior=hn]
#
# Outputs (output/tables):
#   timeblock_thesis_table.csv        every number, machine readable
#   timeblock_thesis_table_body.tex   the LaTeX body for the thesis
#########################################################################

rm(list = ls())

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0L) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]),
    winslash = "/", mustWork = FALSE
  )
  if (dir.exists(dirname(script_path))) setwd(dirname(script_path))
}
source("./time_block_helpers.R")

cli_args <- commandArgs(trailingOnly = TRUE)
# The thesis reports the HN arm, so that is the default here; the shared
# helper defaults to the backend prior, which this table never uses.
if (!any(grepl("^--(hn|prior=)", cli_args))) cli_args <- c(cli_args, "--hn")
prior <- tbf_take_prior_flag(cli_args, c("settings", "data"))
flags <- prior$parsed$values
prior_tag <- prior$prior_tag
setting_ids <- if (!is.null(flags$settings) && nzchar(flags$settings)) {
  as.integer(parse_index_expr(flags$settings, "settings"))
} else {
  c(4L, 5L, 7L)
}

tables_dir <- "output/tables"
if (!dir.exists(tables_dir)) {
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
}

PRIMARY <- c(0.90, 0.95)
MLAB <- c("1" = "iid", "2" = "AR2", "4" = "AR1")
WINDOWS <- list("1-5" = 1:5, "1-15" = 1:15)
# display form: math-mode dashes, as in the spatial table
WLAB <- c("1-5" = "$1$--$5$", "1-15" = "$1$--$15$")
CONTRASTS <- list(
  list(a = "AR2", b = "iid", key = "AR2_iid"),
  list(a = "AR2", b = "AR1", key = "AR2_AR1")
)
DESIGN <- c(
  "4" = "AR(2) $(0.80, -0.35)$",
  "5" = "AR(2) $(0.15, 0.80)$",
  "7" = "ARFIMA, $d = 0.45$"
)

# ---- load, with the provenance gates the campaign writes --------------
CACHE <- list()
for (S in setting_ids) {
  f <- tbf_score_cache_file(S, "", prior_tag)
  if (!file.exists(f)) {
    stop(sprintf("score cache not found: %s", f), call. = FALSE)
  }
  e <- new.env(parent = emptyenv())
  load(f, envir = e)
  if (!identical(e$brier_threshold_basis, "full_series")) {
    stop(sprintf("%s: brier_threshold_basis is '%s', expected 'full_series'",
      f, format(e$brier_threshold_basis)), call. = FALSE)
  }
  if (!isTRUE(e$hn_prior)) {
    stop(sprintf("%s: hn_prior is not TRUE", f), call. = FALSE)
  }
  if (anyNA(e$brier.lead)) {
    stop(sprintf("%s: brier.lead has NA cells", f), call. = FALSE)
  }
  CACHE[[as.character(S)]] <- e
  cat(sprintf("loaded %s: %d data sets, methods [%s], %d blocks\n",
    f, length(e$datasets), paste(e$methods, collapse = ","),
    dim(e$brier.lead)[5]))
}

# per data set: average the lead window, then the blocks
per_dataset <- function(e, prob, leads) {
  qi <- which.min(abs(e$probs - prob))
  a <- apply(e$brier.lead[leads, qi, , , , drop = FALSE], c(3, 4), mean)
  colnames(a) <- MLAB[as.character(e$methods)]
  a
}

paired <- function(d) {
  tt <- t.test(d)                                   # two-sided
  wt <- suppressWarnings(wilcox.test(d, exact = FALSE))
  list(delta = mean(d), lo = tt$conf.int[1], hi = tt$conf.int[2],
       t_p = tt$p.value, w_p = wt$p.value, n = length(d))
}

rows <- list()
for (p in PRIMARY) {
  for (S in setting_ids) {
    e <- CACHE[[as.character(S)]]
    for (wn in names(WINDOWS)) {
      a <- per_dataset(e, p, WINDOWS[[wn]])
      lv <- colMeans(a)
      r <- data.frame(prob = p, setting = S, leads = wn, n = nrow(a),
        iid = lv[["iid"]], AR1 = lv[["AR1"]], AR2 = lv[["AR2"]],
        stringsAsFactors = FALSE)
      for (ct in CONTRASTS) {
        st <- paired(a[, ct$a] - a[, ct$b])
        r[[paste0("delta_", ct$key)]] <- st$delta
        r[[paste0("lo_", ct$key)]] <- st$lo
        r[[paste0("hi_", ct$key)]] <- st$hi
        r[[paste0("tp_", ct$key)]] <- st$t_p
        r[[paste0("wp_", ct$key)]] <- st$w_p
      }
      rows[[length(rows) + 1L]] <- r
    }
  }
}
tab <- do.call(rbind, rows)
csv_file <- file.path(tables_dir, "timeblock_thesis_table.csv")
write.csv(tab, csv_file, row.names = FALSE)
cat("wrote ", csv_file, "\n", sep = "")

# ---- LaTeX body -------------------------------------------------------
fmt_lv <- function(x) formatC(x, format = "f", digits = 4)
cell <- function(delta, lo, hi) {
  star <- if (lo > 0 || hi < 0) "^{*}" else ""
  sprintf("$%+.2f%s$ {\\tiny $[%+.2f,\\,%+.2f]$}",
    delta * 1000, star, lo * 1000, hi * 1000)
}

lines <- c(
  "% Auto-generated by timeblock_thesis_table.R -- do not edit by hand.",
  sprintf("%% lambda-prior arm: %s; thresholds from the full series.", prior_tag),
  "\\begin{tabular}{@{}llccccc@{}}",
  "\\toprule",
  "& & \\multicolumn{3}{c}{Mean Brier score} & \\multicolumn{2}{c}{$\\Delta \\times 10^{3}$ [95\\% CI]} \\\\",
  "\\cmidrule(lr){3-5}\\cmidrule(lr){6-7}",
  "Generating design & Leads & i.i.d. & AR(1) & AR(2) & AR(2)$-$i.i.d. & AR(2)$-$AR(1) \\\\"
)
for (pi in seq_along(PRIMARY)) {
  p <- PRIMARY[pi]
  lines <- c(lines, "\\midrule",
    sprintf("\\multicolumn{7}{@{}l}{\\emph{Threshold quantile $%.2f$}} \\\\", p))
  for (S in setting_ids) {
    sub <- tab[tab$prob == p & tab$setting == S, ]
    for (i in seq_len(nrow(sub))) {
      r <- sub[i, ]
      lab <- if (i == 1L) DESIGN[[as.character(S)]] else ""
      lines <- c(lines, sprintf("%s & %s & %s & %s & %s & %s & %s \\\\",
        lab, WLAB[[r$leads]], fmt_lv(r$iid), fmt_lv(r$AR1), fmt_lv(r$AR2),
        cell(r$delta_AR2_iid, r$lo_AR2_iid, r$hi_AR2_iid),
        cell(r$delta_AR2_AR1, r$lo_AR2_AR1, r$hi_AR2_AR1)))
    }
    if (S != setting_ids[length(setting_ids)]) lines <- c(lines, "\\addlinespace")
  }
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}")

body_file <- file.path(tables_dir, "timeblock_thesis_table_body.tex")
writeLines(lines, body_file)
cat("wrote ", body_file, "\n", sep = "")

# ---- what the prose needs, printed for the writer ---------------------
wp <- c(tab$wp_AR2_iid, tab$wp_AR2_AR1)
tp <- c(tab$tp_AR2_iid, tab$tp_AR2_AR1)
cat(sprintf("\ncells in the table: %d\n", 2 * nrow(tab)))
cat(sprintf("intervals excluding zero: %d\n",
  sum(tab$lo_AR2_iid > 0 | tab$hi_AR2_iid < 0) +
    sum(tab$lo_AR2_AR1 > 0 | tab$hi_AR2_AR1 < 0)))
cat(sprintf("two-sided Wilcoxon p range: %.3f to %.3f\n", min(wp), max(wp)))
cat(sprintf("two-sided paired-t p range: %.4f to %.3f\n", min(tp), max(tp)))
cat(sprintf("Bonferroni threshold at %d tests: %.5f; smallest t p = %.5f\n",
  2 * nrow(tab), 0.05 / (2 * nrow(tab)), min(tp)))
