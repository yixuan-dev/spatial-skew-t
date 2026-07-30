#########################################################################
# collect_diag.R -- bind the per-cell diagnostics and posterior summaries
# written by run-settings.R into two tracked ledgers.
#
# run-settings.R writes one small file per (setting, method, dataset) cell
# because its workers run concurrently and would interleave writes into a
# shared file. Those per-cell files are the copies that survive the fit
# deletion, so this script only ever reads them -- it never needs the fits.
#
# Outputs:
#   output/diag/fit_diag_inline_<setting>.csv   one row per (cell, block)
#   output/tables/recovery_expA_<setting>.csv   posterior mean/median/CI,
#                                               truth, bias and coverage
#
# Usage:
#   Rscript collect_diag.R --setting=<id>
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
parsed <- extract_leading_flags(cli_args, c("setting"))
if (is.null(parsed$values$setting) || !nzchar(parsed$values$setting)) {
  stop("collect_diag.R: --setting=<id> is required.", call. = FALSE)
}
setting_id <- as.integer(parse_index_expr(parsed$values$setting, "setting"))

diag_dir <- "output/diag"
tables_dir <- "output/tables"
for (d in c(diag_dir, tables_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

read_all <- function(pattern) {
  files <- Sys.glob(file.path(diag_dir, pattern))
  if (!length(files)) return(NULL)
  do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
}

# ---- diagnostics ledger ------------------------------------------------
diag_tab <- read_all(sprintf("diag_%d-*.csv", setting_id))
if (is.null(diag_tab)) {
  stop(sprintf("no output/diag/diag_%d-*.csv found -- run run-settings.R first",
    setting_id), call. = FALSE)
}
diag_tab <- diag_tab[order(diag_tab$method, diag_tab$dataset, diag_tab$block), ]
diag_file <- file.path(diag_dir, sprintf("fit_diag_inline_%d.csv", setting_id))
write.csv(diag_tab, diag_file, row.names = FALSE)

cat(sprintf("setting %d: %d block-cells from %d cells\n", setting_id,
  nrow(diag_tab), length(unique(paste(diag_tab$method, diag_tab$dataset)))))
cat(sprintf("  lambda: mean %+.2f, range [%+.2f, %+.2f], negative in %d\n",
  mean(diag_tab$lambda), min(diag_tab$lambda), max(diag_tab$lambda),
  sum(diag_tab$lambda < 0)))
# Descriptive only. Assertion C is not a gate in this run and nothing is
# excluded on the basis of these counts -- the HN prior is the protection.
cat(sprintf("  assertion failures (descriptive): A %d, A' %d, B %d, C %d\n",
  sum(!diag_tab$A_zconsist), sum(!diag_tab$Aprime_sdz),
  sum(!diag_tab$B_truth, na.rm = TRUE), sum(!diag_tab$C_spread)))
by_block <- aggregate(cbind(B_fail = !B_truth, C_fail = !C_spread) ~ block + seam,
  data = diag_tab, FUN = sum)
print(by_block, row.names = FALSE)

# ---- parameter recovery ------------------------------------------------
post_tab <- read_all(sprintf("post_%d-*.csv", setting_id))
if (is.null(post_tab)) {
  cat("\nno post_*.csv files -- skipping the recovery table\n")
} else {
  truth <- get_tbf_truth(setting_id)
  # phi is a pseudo-truth under ARFIMA: the Yule-Walker AR(2) projection is
  # what an AR(2) analysis model converges to, not a data-generating value,
  # so it is reported as a reference and coverage is left NA.
  phi_ref <- if (identical(truth$family, "arfima")) {
    tbf_arfima_yw_ar2(truth$d)
  } else {
    c(phi1 = truth$phi[1], phi2 = truth$phi[2])
  }
  phi_is_truth <- !identical(truth$family, "arfima")

  truth_of <- function(param) {
    if (param %in% names(truth) && is.numeric(truth[[param]])) {
      return(c(value = truth[[param]], is_truth = 1))
    }
    if (grepl("^phi\\.", param)) {
      # phi.w has no data-generating value at K = 1 (membership is fixed,
      # so w carries no signal; setup.R Remark 8)
      if (grepl("^phi\\.w", param)) return(c(value = NA_real_, is_truth = 0))
      lag <- if (grepl("1$", param)) 1L else 2L
      return(c(value = unname(phi_ref[lag]), is_truth = as.numeric(phi_is_truth)))
    }
    c(value = NA_real_, is_truth = 0)
  }
  tv <- t(vapply(post_tab$param, truth_of, c(value = 0, is_truth = 0)))
  post_tab$truth <- tv[, "value"]
  post_tab$truth_is_generating <- tv[, "is_truth"] == 1
  post_tab$covered <- ifelse(post_tab$truth_is_generating,
    post_tab$q025 <= post_tab$truth & post_tab$truth <= post_tab$q975, NA)

  post_file <- file.path(tables_dir, sprintf("recovery_expA_%d.csv", setting_id))
  write.csv(post_tab, post_file, row.names = FALSE)

  cat(sprintf("\n==== setting %d parameter recovery (%s) ====\n",
    setting_id, truth$family))
  cat(sprintf("%-10s %5s %8s %9s %9s %9s   %s\n",
    "param", "block", "truth", "post.mean", "bias", "post.sd", "95% CI coverage"))
  for (nm in tbf_post_params()) {
    for (b in sort(unique(post_tab$block))) {
      v <- post_tab[post_tab$param == nm & post_tab$block == b, ]
      v <- v[!is.na(v$mean), ]
      if (!nrow(v)) next
      cov_txt <- if (all(is.na(v$covered))) {
        "-- (no generating truth)"
      } else {
        sprintf("%d/%d", sum(v$covered, na.rm = TRUE), sum(!is.na(v$covered)))
      }
      cat(sprintf("%-10s %5d %8s %9.3f %+9.3f %9.3f   %s\n",
        nm, b,
        if (is.na(v$truth[1])) "--" else sprintf("%.2f", v$truth[1]),
        mean(v$mean),
        if (is.na(v$truth[1])) NA_real_ else mean(v$mean) - v$truth[1],
        mean(v$sd), cov_txt))
    }
  }
  cat("\nwrote ", post_file, "\n", sep = "")
}

cat("wrote ", diag_file, "\n", sep = "")
