#########################################################################
# expA_breakdown.R -- per-block and per-dataset Brier breakdowns for the
# Experiment A campaign, under BOTH threshold rules, from ONE cache.
#
# scores.R scores every cell twice from the same predictive draws:
#   brier.lead        thresholds = quantile(y[, , set, setting], p)
#                     -- the dataset's FULL SERIES. Primary rule.
#   brier.lead.blockq thresholds = quantile(y_val, p)
#                     -- the validation block's own window. Demoted:
#                     look-ahead, and it pins the exceedance base rate.
# Both arrays therefore describe the SAME fits and differ ONLY in the
# threshold rule, which is what makes the contrast in section 6 of the
# blockq report a controlled comparison rather than two campaigns.
#
# Emits, for each rule R in {full, blockq}:
#   output/tables/expA_block_<R>.csv    + _body.tex   block x method
#   output/tables/expA_dataset_<R>.csv  + _body.tex   dataset x method
#   output/tables/expA_tests_<R>.csv    + _body.tex   paired one-sided
#   output/tables/expA_qsweep_<R>.csv   + _body.tex   all 11 quantiles
# and once:
#   output/tables/expA_coverage.csv     + _body.tex   degeneracy ledger
#   output/tables/expA_rulecontrast.csv + _body.tex   the two rules paired
#   output/tables/expA_pinning.csv      + _body.tex   base-rate pinning
#
# Pairing unit is the DATASET (n = 10): the lead window is averaged
# first (it is the endpoint, not a replicate axis), then the blocks,
# leaving one number per dataset. Same convention as expA_threeway.R.
#
# No CRPS. These reports are Brier-only by design; CRPS is invariant to
# the threshold rule and is recorded in expA_lead_<setting>.csv and in
# tex/timeblock_expABC_legacy.
#
# Usage:
#   Rscript expA_breakdown.R [--settings=5,7] [--data=<path>]
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
parsed <- extract_leading_flags(cli_args, c("settings", "data"))
flags <- parsed$values

# Plural --settings, not the usual --setting=<id>: every table below is a
# cross-setting panel, so a per-setting invocation would force the report
# to hand-assemble panels and reintroduce the transcription this script
# exists to remove.
setting_ids <- if (!is.null(flags$settings) && nzchar(flags$settings)) {
  as.integer(parse_index_expr(flags$settings, "settings"))
} else {
  c(5L, 7L)
}
data_suffix <- if (!is.null(flags$data) && nzchar(flags$data)) {
  derive_data_suffix(flags$data)
} else {
  ""
}

tables_dir <- "output/tables"
if (!dir.exists(tables_dir)) dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# Co-primary thresholds. q0.95 for continuity with every earlier table;
# q0.90 because its coverage is materially better (see the ledger) and
# the two disagree substantively at setting 5.
PRIMARY <- c(0.90, 0.95)
RULES <- c(full = "brier.lead", blockq = "brier.lead.blockq")
RULE_LAB <- c(full = "full-series", blockq = "per-block")
MLAB <- c("1" = "iid", "2" = "AR2", "4" = "AR1")
DISPLAY <- c("1", "4", "2")                     # iid, AR(1), AR(2)
CONTRASTS <- list(c("4", "1", "AR1-iid"),
                  c("2", "1", "AR2-iid"),
                  c("2", "4", "AR2-AR1"))
WINDOWS <- list("1-5" = 1:5, "1-15" = NULL)     # 1-15 filled from block_H

# ---- load every setting's cache, with the provenance gates ------------
CACHE <- list()
for (S in setting_ids) {
  f <- file.path("output/results", sprintf("scores%d%s.RData", S, data_suffix))
  if (!file.exists(f)) stop(sprintf("score cache not found: %s", f), call. = FALSE)
  e <- new.env(parent = emptyenv())
  load(f, envir = e)

  # The tripwire scores.R writes. If it ever reads anything else, the
  # brier.lead array is not the full-series rule and every "full" table
  # below would be mislabelled -- fail loudly rather than emit it.
  if (!identical(e$brier_threshold_basis, "full_series")) {
    stop(sprintf("%s: brier_threshold_basis is '%s', expected 'full_series'",
      f, format(e$brier_threshold_basis)), call. = FALSE)
  }
  for (nm in unname(RULES)) {
    if (is.null(e[[nm]])) stop(sprintf("%s: %s is missing", f, nm), call. = FALSE)
    if (anyNA(e[[nm]])) {
      stop(sprintf("%s: %s has %d NA cells -- the campaign is incomplete",
        f, nm, sum(is.na(e[[nm]]))), call. = FALSE)
    }
  }
  if (!all(PRIMARY %in% e$probs)) {
    stop(sprintf("%s: probs grid lacks %s", f,
      paste(setdiff(PRIMARY, e$probs), collapse = ", ")), call. = FALSE)
  }
  CACHE[[as.character(S)]] <- e
  cat(sprintf("loaded %s: %d datasets, methods [%s], %d blocks, hn_prior=%s\n",
    f, length(e$datasets), paste(e$methods, collapse = ","),
    dim(e$crps.lead)[4], format(e$hn_prior)))
}

H <- CACHE[[1]]$block_H
WINDOWS[["1-15"]] <- seq_len(H)

# arr(e, rule, prob) -> [lead, dataset, method, block] at one threshold
arr_at <- function(e, rule, p) {
  qi <- match(p, e$probs)
  e[[RULES[[rule]]]][, qi, , , ]
}

one_sided <- function(d) {
  d <- d[is.finite(d)]
  if (length(d) < 3L || isTRUE(all.equal(sd(d), 0))) {
    return(c(n = length(d), gap = if (length(d)) mean(d) else NA_real_,
             t_p = NA_real_, w_p = NA_real_))
  }
  c(n = length(d), gap = mean(d),
    t_p = tryCatch(t.test(d, alternative = "less")$p.value,
                   error = function(err) NA_real_),
    w_p = tryCatch(suppressWarnings(
      wilcox.test(d, alternative = "less", exact = FALSE)$p.value),
      error = function(err) NA_real_))
}

# Paired differences at the dataset unit: average the lead window, then
# the blocks, leaving one value per dataset. H1 is "mA better", i.e. the
# difference is negative, hence alternative = "less" above.
paired_vec <- function(arr, mA, mB, leads) {
  a <- apply(arr[leads, , mA, , drop = FALSE], c(2, 4), mean, na.rm = TRUE)
  b <- apply(arr[leads, , mB, , drop = FALSE], c(2, 4), mean, na.rm = TRUE)
  rowMeans(a - b, na.rm = TRUE)
}

# ---- LaTeX helpers ----------------------------------------------------
# Gaps are printed in units of 1e-3 throughout: the levels are ~0.05-0.13
# and the gaps ~0.0002-0.03, so a shared format either loses the gaps or
# doubles the column width.
fmt <- function(x, d = 4) ifelse(is.na(x), "---", formatC(x, format = "f", digits = d))
fmt_gap <- function(x, d = 1) {
  ifelse(is.na(x), "---",
    sprintf("$%s%s$", ifelse(x < 0, "-", "+"),
      formatC(abs(x) * 1000, format = "f", digits = d)))
}
fmt_p <- function(x, d = 3) {
  ifelse(is.na(x), "---",
    ifelse(x < 0.001, "$<\\!0.001$",
      ifelse(x < 0.05,
        sprintf("$\\mathbf{%s}$", formatC(x, format = "f", digits = d)),
        formatC(x, format = "f", digits = d))))
}
write_body <- function(lines, file) {
  writeLines(c("% Auto-generated by expA_breakdown.R -- do not edit by hand.",
    lines), file)
  cat("wrote ", file, "\n", sep = "")
}
setting_desc <- function(S) {
  switch(as.character(S),
    "5" = "near-unit-root, $\\phi=(0.15,0.80)$",
    "7" = "long memory, ARFIMA$(0,0.45,0)$",
    sprintf("setting %d", S))
}

# =======================================================================
# 1. per-block and per-dataset breakdowns
# =======================================================================
margin_table <- function(rule, margin) {
  # margin = "block" (collapse leads+datasets) or "dataset" (leads+blocks)
  mdim <- if (identical(margin, "block")) 4L else 2L
  rows <- list()
  for (S in setting_ids) {
    e <- CACHE[[as.character(S)]]
    for (p in PRIMARY) {
      a <- arr_at(e, rule, p)
      m <- apply(a, c(3L, mdim), mean, na.rm = TRUE)   # [method, margin]
      for (k in colnames(m)) {
        rows[[length(rows) + 1L]] <- data.frame(
          setting = S, rule = RULE_LAB[[rule]], prob = p, margin = margin,
          level = k,
          iid = m["1", k], AR1 = m["4", k], AR2 = m["2", k],
          `AR1-iid` = m["4", k] - m["1", k],
          `AR2-iid` = m["2", k] - m["1", k],
          `AR2-AR1` = m["2", k] - m["4", k],
          check.names = FALSE)
      }
    }
  }
  do.call(rbind, rows)
}

margin_body <- function(tab, rule, margin) {
  lab <- if (identical(margin, "block")) "block" else "data set"
  head <- c(
    "\\begin{table}[htbp]", "\\centering",
    sprintf("\\caption{Experiment~A, %s threshold rule: mean Brier score by %s,",
      RULE_LAB[[rule]], lab),
    sprintf("leads $1$--$%d$ and the remaining axis pooled. Levels in native units;", H),
    "the AR(2)$-$i.i.d.\\ gap in units of $10^{-3}$ (negative favours AR(2)).",
    "Co-primary thresholds $q_{90}$ and $q_{95}$ side by side.}",
    sprintf("\\label{tab:%s-%s}", margin, rule),
    "\\small",
    "\\begin{tabular}{clrrrr rrrr}", "\\toprule",
    sprintf(" & & \\multicolumn{4}{c}{$q_{90}$} & \\multicolumn{4}{c}{$q_{95}$} \\\\"),
    "\\cmidrule(lr){3-6}\\cmidrule(lr){7-10}",
    sprintf("setting & %s & i.i.d. & AR(1) & AR(2) & gap & i.i.d. & AR(1) & AR(2) & gap \\\\",
      lab),
    "\\midrule")
  body <- character(0)
  for (S in setting_ids) {
    sub <- tab[tab$setting == S, ]
    lv <- unique(sub$level)
    for (i in seq_along(lv)) {
      r90 <- sub[sub$prob == 0.90 & sub$level == lv[i], ]
      r95 <- sub[sub$prob == 0.95 & sub$level == lv[i], ]
      first <- if (i == 1L) sprintf("\\multirow{%d}{*}{%d}", length(lv), S) else ""
      body <- c(body, sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
        first, lv[i],
        fmt(r90$iid), fmt(r90$AR1), fmt(r90$AR2), fmt_gap(r90$`AR2-iid`),
        fmt(r95$iid), fmt(r95$AR1), fmt(r95$AR2), fmt_gap(r95$`AR2-iid`)))
    }
    if (S != setting_ids[length(setting_ids)]) body <- c(body, "\\midrule")
  }
  c(head, body, "\\bottomrule", "\\end{tabular}", "\\end{table}")
}

for (rule in names(RULES)) {
  for (margin in c("block", "dataset")) {
    tab <- margin_table(rule, margin)
    write.csv(tab, file.path(tables_dir,
      sprintf("expA_%s_%s%s.csv", margin, rule, data_suffix)), row.names = FALSE)
    write_body(margin_body(tab, rule, margin), file.path(tables_dir,
      sprintf("expA_%s_%s_body.tex", margin, rule)))
  }
}

# =======================================================================
# 2. paired one-sided tests, dataset unit
# =======================================================================
tests_table <- function(rule) {
  rows <- list()
  for (S in setting_ids) {
    e <- CACHE[[as.character(S)]]
    for (ct in CONTRASTS) {
      if (!all(as.integer(ct[1:2]) %in% e$methods)) next
      for (wn in names(WINDOWS)) for (p in PRIMARY) {
        st <- one_sided(paired_vec(arr_at(e, rule, p), ct[1], ct[2], WINDOWS[[wn]]))
        rows[[length(rows) + 1L]] <- data.frame(
          setting = S, rule = RULE_LAB[[rule]], contrast = ct[3], window = wn,
          prob = p, n = st[["n"]], gap = st[["gap"]],
          t_p = st[["t_p"]], w_p = st[["w_p"]])
      }
    }
  }
  do.call(rbind, rows)
}

tests_body <- function(tab, rule) {
  head <- c(
    "\\begin{table}[htbp]", "\\centering",
    sprintf("\\caption{Experiment~A, %s threshold rule: paired one-sided tests",
      RULE_LAB[[rule]]),
    "($H_1$: the first method of the contrast scores lower), data-set unit",
    "$n=10$ --- the lead window is averaged first, then the five blocks.",
    "Gaps in units of $10^{-3}$; $p<0.05$ in bold.}",
    sprintf("\\label{tab:tests-%s}", rule),
    "\\small",
    "\\begin{tabular}{cllrrr rrr}", "\\toprule",
    " & & & \\multicolumn{3}{c}{$q_{90}$} & \\multicolumn{3}{c}{$q_{95}$} \\\\",
    "\\cmidrule(lr){4-6}\\cmidrule(lr){7-9}",
    "setting & contrast & leads & gap & $t$ $p$ & $W$ $p$ & gap & $t$ $p$ & $W$ $p$ \\\\",
    "\\midrule")
  body <- character(0)
  for (S in setting_ids) {
    sub <- tab[tab$setting == S, ]
    key <- unique(sub[, c("contrast", "window")])
    for (i in seq_len(nrow(key))) {
      r90 <- sub[sub$contrast == key$contrast[i] & sub$window == key$window[i] &
                   sub$prob == 0.90, ]
      r95 <- sub[sub$contrast == key$contrast[i] & sub$window == key$window[i] &
                   sub$prob == 0.95, ]
      first <- if (i == 1L) sprintf("\\multirow{%d}{*}{%d}", nrow(key), S) else ""
      cl <- gsub("-", "$-$", key$contrast[i], fixed = TRUE)
      cl <- gsub("AR1", "AR(1)", gsub("AR2", "AR(2)", cl))
      cl <- gsub("iid", "i.i.d.", cl)
      body <- c(body, sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
        first, cl, sub("-", "--", key$window[i], fixed = TRUE),
        fmt_gap(r90$gap), fmt_p(r90$t_p), fmt_p(r90$w_p),
        fmt_gap(r95$gap), fmt_p(r95$t_p), fmt_p(r95$w_p)))
    }
    if (S != setting_ids[length(setting_ids)]) body <- c(body, "\\midrule")
  }
  c(head, body, "\\bottomrule", "\\end{tabular}", "\\end{table}")
}

TESTS <- list()
for (rule in names(RULES)) {
  tab <- tests_table(rule)
  TESTS[[rule]] <- tab
  write.csv(tab, file.path(tables_dir,
    sprintf("expA_tests_%s%s.csv", rule, data_suffix)), row.names = FALSE)
  write_body(tests_body(tab, rule),
    file.path(tables_dir, sprintf("expA_tests_%s_body.tex", rule)))
}

# =======================================================================
# 3. quantile sweep over the full probs grid
# =======================================================================
sweep_table <- function(rule) {
  rows <- list()
  for (S in setting_ids) {
    e <- CACHE[[as.character(S)]]
    for (p in e$probs) {
      a <- arr_at(e, rule, p)
      lv <- apply(a, 3, mean, na.rm = TRUE)
      st <- one_sided(paired_vec(a, "2", "1", WINDOWS[["1-15"]]))
      pd <- apply(a, c(2, 3), mean, na.rm = TRUE)
      pb <- apply(a, c(3, 4), mean, na.rm = TRUE)
      rows[[length(rows) + 1L]] <- data.frame(
        setting = S, rule = RULE_LAB[[rule]], prob = p,
        iid = lv["1"], AR1 = lv["4"], AR2 = lv["2"],
        gap_AR2_iid = st[["gap"]], t_p = st[["t_p"]], w_p = st[["w_p"]],
        datasets_won = sum(pd[, "2"] - pd[, "1"] < 0),
        blocks_won = sum(pb["2", ] - pb["1", ] < 0))
    }
  }
  do.call(rbind, rows)
}

sweep_body <- function(tab, rule) {
  head <- c(
    "\\begin{table}[htbp]", "\\centering",
    sprintf("\\caption{Experiment~A, %s threshold rule: the whole threshold grid,",
      RULE_LAB[[rule]]),
    "leads $1$--$15$. Gap $=$ AR(2)$-$i.i.d.\\ in units of $10^{-3}$, data-set",
    "unit $n=10$. ``won'' counts the data sets (of $10$) and blocks (of $5$)",
    "in which AR(2) scores below i.i.d.}",
    sprintf("\\label{tab:qsweep-%s}", rule),
    "\\small",
    "\\begin{tabular}{clrrrrrrcc}", "\\toprule",
    "setting & $p$ & i.i.d. & AR(1) & AR(2) & gap & $t$ $p$ & $W$ $p$ & ds won & blk won \\\\",
    "\\midrule")
  body <- character(0)
  for (S in setting_ids) {
    sub <- tab[tab$setting == S, ]
    for (i in seq_len(nrow(sub))) {
      r <- sub[i, ]
      first <- if (i == 1L) sprintf("\\multirow{%d}{*}{%d}", nrow(sub), S) else ""
      body <- c(body, sprintf(
        "%s & %s & %s & %s & %s & %s & %s & %s & %d/10 & %d/5 \\\\",
        first, formatC(r$prob, format = "f", digits = 3),
        fmt(r$iid), fmt(r$AR1), fmt(r$AR2), fmt_gap(r$gap_AR2_iid),
        fmt_p(r$t_p), fmt_p(r$w_p), r$datasets_won, r$blocks_won))
    }
    if (S != setting_ids[length(setting_ids)]) body <- c(body, "\\midrule")
  }
  c(head, body, "\\bottomrule", "\\end{tabular}", "\\end{table}")
}

for (rule in names(RULES)) {
  tab <- sweep_table(rule)
  write.csv(tab, file.path(tables_dir,
    sprintf("expA_qsweep_%s%s.csv", rule, data_suffix)), row.names = FALSE)
  write_body(sweep_body(tab, rule),
    file.path(tables_dir, sprintf("expA_qsweep_%s_body.tex", rule)))
}

# =======================================================================
# 4. coverage / degeneracy ledger (full-series rule)
# =======================================================================
# exceed.rate.lead is a function of the DATA and the full-series
# thresholds only, so it is complete regardless of which fits exist. A
# cell whose realised rate is 0 contributes only the predicted
# probability to the Brier score: there is no event to discriminate, so
# the score there has degenerated onto climatology.
cov_rows <- list()
for (S in setting_ids) {
  e <- CACHE[[as.character(S)]]
  for (qi in seq_along(e$probs)) {
    er <- e$exceed.rate.lead[, qi, , ]
    px <- e$pexceed.mean[, qi, , , ]
    cov_rows[[length(cov_rows) + 1L]] <- data.frame(
      setting = S, prob = e$probs[qi],
      zero_share = mean(er == 0), realised_rate = mean(er),
      pred_rate_iid = mean(px[, , "1", ]), pred_rate_AR2 = mean(px[, , "2", ]),
      threshold_mean = mean(e$brier.thresholds[qi, ]))
  }
}
cov_tab <- do.call(rbind, cov_rows)
write.csv(cov_tab, file.path(tables_dir, sprintf("expA_coverage%s.csv", data_suffix)),
  row.names = FALSE)

cov_body <- c(
  "\\begin{table}[htbp]", "\\centering",
  "\\caption{Coverage and degeneracy ledger, full-series rule. ``zero'' is the",
  "share of (lead, data set, block) cells in which the validation window",
  "realises \\emph{no} exceedance of the fixed threshold, so the Brier score",
  "there degenerates onto climatology. ``realised'' is the mean realised",
  "exceedance rate and ``pred'' the mean predicted exceedance probability",
  "(\\code{pexceed.mean}); predicted far below realised is a level error.}",
  "\\label{tab:coverage}", "\\small",
  "\\begin{tabular}{lrrrr rrrr}", "\\toprule",
  " & \\multicolumn{4}{c}{setting 5} & \\multicolumn{4}{c}{setting 7} \\\\",
  "\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}",
  "$p$ & zero & realised & pred i.i.d. & pred AR(2) & zero & realised & pred i.i.d. & pred AR(2) \\\\",
  "\\midrule")
for (p in sort(unique(cov_tab$prob))) {
  r5 <- cov_tab[cov_tab$setting == 5 & cov_tab$prob == p, ]
  r7 <- cov_tab[cov_tab$setting == 7 & cov_tab$prob == p, ]
  cell <- function(r) if (nrow(r) == 0) "--- & --- & --- & ---" else
    sprintf("%.1f\\%% & %s & %s & %s", 100 * r$zero_share,
      fmt(r$realised_rate), fmt(r$pred_rate_iid), fmt(r$pred_rate_AR2))
  cov_body <- c(cov_body, sprintf("%s & %s & %s \\\\",
    formatC(p, format = "f", digits = 3), cell(r5), cell(r7)))
}
cov_body <- c(cov_body, "\\bottomrule", "\\end{tabular}", "\\end{table}")
write_body(cov_body, file.path(tables_dir, "expA_coverage_body.tex"))

# =======================================================================
# 5. the rule contrast -- same fits, two thresholds
# =======================================================================
rc <- merge(
  TESTS[["full"]][, c("setting", "contrast", "window", "prob", "n", "gap", "t_p", "w_p")],
  TESTS[["blockq"]][, c("setting", "contrast", "window", "prob", "gap", "t_p", "w_p")],
  by = c("setting", "contrast", "window", "prob"),
  suffixes = c("_full", "_blockq"))
rc <- rc[order(rc$setting, rc$contrast, rc$window, rc$prob), ]
write.csv(rc, file.path(tables_dir, sprintf("expA_rulecontrast%s.csv", data_suffix)),
  row.names = FALSE)

rc_body <- c(
  "\\begin{table}[htbp]", "\\centering",
  "\\caption{The two threshold rules on \\emph{identical fits}: paired",
  "one-sided tests at the data-set unit ($n=10$), leads $1$--$15$. The",
  "per-block rule does not shrink the gap --- at setting 5, $q_{95}$ it is",
  "\\emph{twice} the full-series gap --- yet it reaches no significance,",
  "because the threshold is a statistic of the very window being scored and",
  "so injects a per-block random component that breaks the pairing.",
  "Gaps in units of $10^{-3}$; $p<0.05$ in bold.}",
  "\\label{tab:rulecontrast}", "\\small",
  "\\begin{tabular}{clcrr rr}", "\\toprule",
  " & & & \\multicolumn{2}{c}{full-series} & \\multicolumn{2}{c}{per-block} \\\\",
  "\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}",
  "setting & contrast & $p$ & gap & $t$ $p$ & gap & $t$ $p$ \\\\",
  "\\midrule")
sub_all <- rc[rc$window == "1-15", ]
for (S in setting_ids) {
  sub <- sub_all[sub_all$setting == S, ]
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    first <- if (i == 1L) sprintf("\\multirow{%d}{*}{%d}", nrow(sub), S) else ""
    cl <- gsub("-", "$-$", r$contrast, fixed = TRUE)
    cl <- gsub("AR1", "AR(1)", gsub("AR2", "AR(2)", cl))
    cl <- gsub("iid", "i.i.d.", cl)
    rc_body <- c(rc_body, sprintf("%s & %s & %s & %s & %s & %s & %s \\\\",
      first, cl, formatC(r$prob, format = "f", digits = 2),
      fmt_gap(r$gap_full), fmt_p(r$t_p_full),
      fmt_gap(r$gap_blockq), fmt_p(r$t_p_blockq)))
  }
  if (S != setting_ids[length(setting_ids)]) rc_body <- c(rc_body, "\\midrule")
}
rc_body <- c(rc_body, "\\bottomrule", "\\end{tabular}", "\\end{table}")
write_body(rc_body, file.path(tables_dir, "expA_rulecontrast_body.tex"))

# =======================================================================
# 6. base-rate pinning -- recomputed from the data, not the fits
# =======================================================================
# The claim is that quantile(y_val, p) forces the realised exceedance rate
# of the block to 1-p regardless of what the block actually did. That is
# a statement about the DATA alone, so it is checked here directly rather
# than asserted.
data_path <- resolve_simstudy_data_path(flags$data)
denv <- new.env(parent = emptyenv())
load(data_path, envir = denv)
blocks <- tbf_blocks(denv$block_seams, denv$block_H, denv$nt)

pin_rows <- list()
for (S in setting_ids) {
  e <- CACHE[[as.character(S)]]
  for (p in e$probs) {
    blk_rate <- full_rate <- numeric(0)
    for (d in as.integer(e$datasets)) {
      thr_full <- quantile(denv$y[, , d, S], probs = p, na.rm = TRUE, names = FALSE)
      for (b in seq_along(blocks)) {
        yv <- denv$y[, blocks[[b]]$test_times, d, S]
        blk_rate <- c(blk_rate, mean(yv > quantile(yv, probs = p, na.rm = TRUE,
          names = FALSE), na.rm = TRUE))
        full_rate <- c(full_rate, mean(yv > thr_full, na.rm = TRUE))
      }
    }
    pin_rows[[length(pin_rows) + 1L]] <- data.frame(
      setting = S, prob = p, target = 1 - p,
      blockq_mean = mean(blk_rate), blockq_min = min(blk_rate),
      blockq_max = max(blk_rate), blockq_sd = sd(blk_rate),
      full_mean = mean(full_rate), full_min = min(full_rate),
      full_max = max(full_rate), full_sd = sd(full_rate))
  }
}
pin_tab <- do.call(rbind, pin_rows)
write.csv(pin_tab, file.path(tables_dir, sprintf("expA_pinning%s.csv", data_suffix)),
  row.names = FALSE)

pin_body <- c(
  "\\begin{table}[htbp]", "\\centering",
  "\\caption{Base-rate pinning, computed from the data alone (no fits are",
  "involved). For each of the $50$ (data set, block) cells the realised",
  "exceedance rate of the validation window is evaluated under each rule.",
  "The per-block rule holds it at the nominal $1-p$ with almost no spread;",
  "the full-series rule lets it vary with what the block actually did,",
  "which is the whole reason a level error is visible under one rule and",
  "invisible under the other.}",
  "\\label{tab:pinning}", "\\small",
  "\\begin{tabular}{clrrr rrr}", "\\toprule",
  " & & & \\multicolumn{2}{c}{per-block rule} & \\multicolumn{3}{c}{full-series rule} \\\\",
  "\\cmidrule(lr){4-5}\\cmidrule(lr){6-8}",
  "setting & $p$ & $1-p$ & mean & sd & mean & sd & range \\\\",
  "\\midrule")
for (S in setting_ids) {
  sub <- pin_tab[pin_tab$setting == S & pin_tab$prob %in% PRIMARY, ]
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    first <- if (i == 1L) sprintf("\\multirow{%d}{*}{%d}", nrow(sub), S) else ""
    pin_body <- c(pin_body, sprintf(
      "%s & %s & %s & %s & %s & %s & %s & %s--%s \\\\",
      first, formatC(r$prob, format = "f", digits = 2), fmt(r$target, 3),
      fmt(r$blockq_mean, 4), fmt(r$blockq_sd, 4),
      fmt(r$full_mean, 4), fmt(r$full_sd, 4),
      fmt(r$full_min, 3), fmt(r$full_max, 3)))
  }
  if (S != setting_ids[length(setting_ids)]) pin_body <- c(pin_body, "\\midrule")
}
pin_body <- c(pin_body, "\\bottomrule", "\\end{tabular}", "\\end{table}")
write_body(pin_body, file.path(tables_dir, "expA_pinning_body.tex"))

# =======================================================================
# 7. verification against expA_threeway.R
# =======================================================================
# expA_threeway.R writes the same q95 dataset-unit tests from the same
# cache. If this script's aggregation disagrees with it, one of the two
# is wrong -- so check rather than trust.
cat("\n==== verification: q95 tests vs expA_tests_<setting>.csv ====\n")
worst <- 0
for (S in setting_ids) {
  ref_file <- file.path(tables_dir, sprintf("expA_tests_%d%s.csv", S, data_suffix))
  if (!file.exists(ref_file)) {
    cat(sprintf("  setting %d: %s absent, skipped\n", S, basename(ref_file)))
    next
  }
  ref <- read.csv(ref_file, stringsAsFactors = FALSE)
  ref <- ref[ref$score == "Brier95" & ref$unit == "dataset", ]
  mine <- TESTS[["full"]]
  mine <- mine[mine$setting == S & mine$prob == 0.95, ]
  for (i in seq_len(nrow(ref))) {
    m <- mine[mine$contrast == ref$contrast[i] & mine$window == ref$window[i], ]
    if (nrow(m) != 1L) { cat(sprintf("  MISSING %s %s\n", ref$contrast[i], ref$window[i])); next }
    d <- max(abs(c(m$gap - ref$gap[i], m$t_p - ref$t_p[i], m$w_p - ref$w_p[i])))
    worst <- max(worst, d)
    cat(sprintf("  setting %d  %-9s %-5s  max|diff| = %.3e %s\n",
      S, ref$contrast[i], ref$window[i], d, if (d < 1e-12) "OK" else "MISMATCH"))
  }
}
cat(sprintf("\nworst absolute disagreement: %.3e -- %s\n", worst,
  if (worst < 1e-12) "reproduces expA_threeway.R exactly" else "INVESTIGATE"))

cat("\n==== headline (full-series, dataset unit n=10, leads 1-15) ====\n")
h <- TESTS[["full"]]
h <- h[h$window == "1-15" & h$contrast %in% c("AR2-iid", "AR2-AR1"), ]
for (i in seq_len(nrow(h))) {
  cat(sprintf("  s%d %-8s q%.2f  gap %+8.5f  t p %.4f  W p %.4f\n",
    h$setting[i], h$contrast[i], h$prob[i], h$gap[i], h$t_p[i], h$w_p[i]))
}
