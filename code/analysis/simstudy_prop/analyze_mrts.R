#########################################################################
# analyze_mrts.R
#
# Score fits in `fits/` (filenames `<setting>-<method>-<dataset>-p<K>.RData`)
# and analyze how Brier / Quantile scores vary with the MRTS basis rank K.
#
# Usage:
#   Rscript analyze_mrts.R                  # all fits found
#   Rscript analyze_mrts.R --setting=4 \    # restrict to one (setting,method,dataset)
#                          --method=2 \
#                          --dataset=1
#
# Outputs (under output/tables/ and output/results/):
#   - mrts_brier_long.csv          per-K x per-quantile Brier + Quantile scores
#   - mrts_brier_summary.csv       per-K summary (mean over all quantiles, plus bulk/tail bands)
#   - mrts_brier_analysis.RData    R objects: brier_long, brier_summary, brier_wide
#########################################################################

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
source("../../R/prop/auxfunctions.R")

# ---- locate simdata --------------------------------------------------
resolve_simdata <- function() {
  candidates <- c("./simdata.RData", "../simstudy/simdata.RData")
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0L) {
    stop("Unable to locate simdata.RData (looked in ./ and ../simstudy/).", call. = FALSE)
  }
  normalizePath(existing[1], winslash = "/", mustWork = TRUE)
}
load(resolve_simdata())

# ---- CLI / env -------------------------------------------------------
parse_kv_arg <- function(args, key) {
  pat <- sprintf("^--%s=", key)
  hit <- grep(pat, args, value = TRUE)
  if (length(hit) == 0L) return(NA_integer_)
  as.integer(sub(pat, "", hit[1]))
}

cli_args <- commandArgs(trailingOnly = TRUE)
filter_setting <- parse_kv_arg(cli_args, "setting")
filter_method  <- parse_kv_arg(cli_args, "method")
filter_dataset <- parse_kv_arg(cli_args, "dataset")

fits_dir <- trimws(Sys.getenv("SIMSTUDY_PROP_FITS_DIR", unset = "fits"))
if (!nzchar(fits_dir)) fits_dir <- "fits"
if (!dir.exists(fits_dir)) {
  stop(sprintf("fits directory not found: %s", fits_dir), call. = FALSE)
}

output_root <- "output"
tables_dir  <- file.path(output_root, "tables")
results_dir <- file.path(output_root, "results")
for (d in c(tables_dir, results_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ---- discover fit files ---------------------------------------------
fit_files <- list.files(fits_dir, pattern = "^\\d+-\\d+-\\d+-[pP]\\d+\\.RData$", full.names = TRUE)
if (length(fit_files) == 0L) {
  stop(sprintf("No fits matching <setting>-<method>-<dataset>-p<K>.RData under %s", fits_dir),
       call. = FALSE)
}

parse_fit_name <- function(path) {
  base <- sub("\\.RData$", "", basename(path))
  m <- regmatches(base, regexec("^(\\d+)-(\\d+)-(\\d+)-[pP](\\d+)$", base))[[1]]
  if (length(m) != 5L) return(NULL)
  data.frame(
    file       = path,
    setting    = as.integer(m[2]),
    method_id  = as.integer(m[3]),
    dataset    = as.integer(m[4]),
    prop_k     = as.integer(m[5]),
    stringsAsFactors = FALSE
  )
}

index <- do.call(rbind, lapply(fit_files, parse_fit_name))
if (!is.na(filter_setting)) index <- index[index$setting == filter_setting, , drop = FALSE]
if (!is.na(filter_method))  index <- index[index$method_id == filter_method, , drop = FALSE]
if (!is.na(filter_dataset)) index <- index[index$dataset == filter_dataset, , drop = FALSE]
if (nrow(index) == 0L) {
  stop("No fits left after filtering on --setting/--method/--dataset.", call. = FALSE)
}
index <- index[order(index$setting, index$method_id, index$dataset, index$prop_k), , drop = FALSE]

# ---- shared scoring set-up ------------------------------------------
catalog <- get_prop_method_catalog()
obs <- c(rep(TRUE, nrow(y) - ntest), rep(FALSE, ntest))
probs <- c(0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995)
classify_band <- function(q) {
  b <- rep("other", length(q))
  b[q >= 0.90 & q <= 0.95] <- "bulk"
  b[q >= 0.98] <- "tail"
  b
}
bands <- classify_band(probs)

load_yp <- function(path) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (exists("fit.1", envir = env, inherits = FALSE)) {
    obj <- get("fit.1", envir = env)
  } else if (exists("fit", envir = env, inherits = FALSE)) {
    obj <- get("fit", envir = env)
  } else {
    return(NULL)
  }
  if (is.null(obj$yp)) NULL else obj$yp
}

# ---- score each fit --------------------------------------------------
long_rows <- vector("list", nrow(index))
elapsed_lookup <- rep(NA_real_, nrow(index))

for (i in seq_len(nrow(index))) {
  row <- index[i, , drop = FALSE]
  cat(sprintf("[%d/%d] scoring %s\n", i, nrow(index), basename(row$file)))

  yp <- load_yp(row$file)
  if (is.null(yp)) {
    cat("  -> no yp found, skipping\n"); next
  }

  validate   <- y[!obs, , row$dataset, row$setting]
  thresholds <- quantile(y[, , row$dataset, row$setting], probs = probs, na.rm = TRUE)

  brier <- BrierScore(yp, thresholds, validate)
  quant <- QuantScore(yp, probs, validate)

  # pull elapsed time if present
  env <- new.env(parent = emptyenv()); load(row$file, envir = env)
  if (exists("runtime_info", envir = env, inherits = FALSE)) {
    rt <- get("runtime_info", envir = env)
    if (!is.null(rt$elapsed_sec)) elapsed_lookup[i] <- as.numeric(rt$elapsed_sec)
  }

  cat_row <- catalog[catalog$method_id == row$method_id, , drop = FALSE]
  method_label <- if (nrow(cat_row) == 1L) cat_row$label else as.character(row$method_id)

  long_rows[[i]] <- data.frame(
    setting     = row$setting,
    method_id   = row$method_id,
    method_label = method_label,
    dataset     = row$dataset,
    prop_k      = row$prop_k,
    quantile    = probs,
    band        = bands,
    threshold   = unname(thresholds),
    brier       = brier,
    quant_score = quant,
    elapsed_sec = elapsed_lookup[i],
    stringsAsFactors = FALSE
  )
}

brier_long <- do.call(rbind, long_rows)
if (is.null(brier_long) || nrow(brier_long) == 0L) {
  stop("No fits produced scores.", call. = FALSE)
}

# ---- per-K summaries -------------------------------------------------
safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

aggregate_by <- function(df, extra_cols) {
  key_cols <- c("setting", "method_id", "method_label", "dataset", "prop_k", extra_cols)
  key <- do.call(paste, c(df[key_cols], sep = "\r"))
  out <- lapply(split(df, key), function(piece) {
    base <- piece[1, key_cols, drop = FALSE]
    base$n_quantiles <- nrow(piece)
    base$mean_brier  <- safe_mean(piece$brier)
    base$mean_quant  <- safe_mean(piece$quant_score)
    base$elapsed_sec <- safe_mean(piece$elapsed_sec)
    base
  })
  res <- do.call(rbind, out)
  res[order(res$setting, res$method_id, res$dataset, res$prop_k,
            if ("band" %in% extra_cols) res$band else rep("", nrow(res))), , drop = FALSE]
}

summary_overall <- aggregate_by(brier_long, character(0))
summary_overall$band <- "all"
summary_band <- aggregate_by(
  brier_long[brier_long$band %in% c("bulk", "tail"), , drop = FALSE],
  "band"
)
brier_summary <- rbind(summary_overall, summary_band)
brier_summary <- brier_summary[order(brier_summary$setting, brier_summary$method_id,
                                     brier_summary$dataset, brier_summary$prop_k,
                                     brier_summary$band), , drop = FALSE]
rownames(brier_summary) <- NULL

# wide form: prop_k x quantile -> brier (per setting/method/dataset block)
make_wide <- function(df, value_col) {
  ks  <- sort(unique(df$prop_k))
  qs  <- sort(unique(df$quantile))
  m   <- matrix(NA_real_, length(ks), length(qs),
                dimnames = list(prop_k = as.character(ks), quantile = as.character(qs)))
  for (i in seq_len(nrow(df))) {
    m[as.character(df$prop_k[i]), as.character(df$quantile[i])] <- df[[value_col]][i]
  }
  m
}
brier_wide <- by(
  brier_long,
  list(brier_long$setting, brier_long$method_id, brier_long$dataset),
  function(piece) {
    list(
      setting   = piece$setting[1],
      method_id = piece$method_id[1],
      dataset   = piece$dataset[1],
      brier     = make_wide(piece, "brier"),
      quant     = make_wide(piece, "quant_score")
    )
  },
  simplify = FALSE
)

# ---- write outputs ---------------------------------------------------
write.csv(brier_long,   file.path(tables_dir, "mrts_brier_long.csv"),    row.names = FALSE)
write.csv(brier_summary, file.path(tables_dir, "mrts_brier_summary.csv"), row.names = FALSE)
save(brier_long, brier_summary, brier_wide,
     file = file.path(results_dir, "mrts_brier_analysis.RData"))

# ---- console report --------------------------------------------------
cat("\n=== Per-K Brier / Quantile summary (mean over all quantiles) ===\n")
print_overall <- summary_overall[, c("setting", "method_id", "dataset", "prop_k",
                                     "mean_brier", "mean_quant", "elapsed_sec")]
print(format(print_overall, digits = 5), row.names = FALSE)

cat("\n=== Per-K summary by band (bulk = q[0.90,0.95], tail = q>=0.98) ===\n")
print_band <- summary_band[, c("setting", "method_id", "dataset", "prop_k", "band",
                               "mean_brier", "mean_quant")]
print(format(print_band, digits = 5), row.names = FALSE)

cat("\nWrote:\n")
cat("  ", file.path(tables_dir, "mrts_brier_long.csv"), "\n")
cat("  ", file.path(tables_dir, "mrts_brier_summary.csv"), "\n")
cat("  ", file.path(results_dir, "mrts_brier_analysis.RData"), "\n")
