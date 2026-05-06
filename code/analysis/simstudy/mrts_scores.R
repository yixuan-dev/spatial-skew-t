## ---------------------------------------------------------------------------
## mrts_scores.R
##
## Stage 1: compute Brier + Quantile scores for every fitted result file
## and save them into per-(setting, K) bundles:
##
##   scores{setting}_{K}mrts.RData
##
## where K = 0 means the no-MRTS baseline.  Stage 2 (mrts_helpfulness.R)
## consumes these bundles.
##
## Each bundle contains arrays:
##   quant.score, brier.score   dim = c(prob, dataset, method)
##   elapsed.sec                dim = c(dataset, method)
## plus metadata: probs, methods, datasets, setting, mrts_k.
##
## NA marks a fit that was not run or could not be loaded.
##
## Run from this directory:
##   Rscript mrts_scores.R
## ---------------------------------------------------------------------------

rm(list = ls())

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]),
                               winslash = "/", mustWork = FALSE)
  script_dir <- dirname(script_path)
  if (dir.exists(script_dir)) setwd(script_dir)
}

source("./mrts_cov_helpers.R")
source("../../R/auxfunctions.R")
load("simdata.RData")

results_dir <- trimws(Sys.getenv("SIMSTUDY_MRTS_RESULTS_DIR", unset = "results"))

probs   <- c(0.90, 0.91, 0.92, 0.93, 0.94, 0.95,
             0.96, 0.97, 0.98, 0.99, 0.995)
methods <- 1:5
datasets <- seq_len(dim(y)[3])
obs <- c(rep(TRUE, nrow(y) - ntest), rep(FALSE, ntest))

## (setting, K) jobs to compute.  K = 0 means no-MRTS baseline.
jobs <- list(
  list(setting = 4L, mrts_k = 0L),
  list(setting = 4L, mrts_k = 10L),
  list(setting = 4L, mrts_k = 15L),
  list(setting = 4L, mrts_k = 20L),
  list(setting = 5L, mrts_k = 0L),
  list(setting = 5L, mrts_k = 15L)
)

dim_vec <- c(length(probs), length(datasets), length(methods))
dimnames_list <- list(
  prob    = as.character(probs),
  dataset = as.character(datasets),
  method  = as.character(methods)
)

for (job in jobs) {
  setting <- job$setting
  mrts_k  <- job$mrts_k
  k_arg   <- if (mrts_k == 0L) NA_integer_ else mrts_k
  out     <- sprintf("scores%d_%dmrts.RData", setting, mrts_k)

  cat(sprintf("\n=== setting %d, K = %d -> %s ===\n", setting, mrts_k, out))

  quant.score <- array(NA_real_, dim = dim_vec, dimnames = dimnames_list)
  brier.score <- array(NA_real_, dim = dim_vec, dimnames = dimnames_list)
  elapsed.sec <- array(NA_real_,
                       dim = dim_vec[2:3],
                       dimnames = dimnames_list[c("dataset", "method")])

  for (d_idx in seq_along(datasets)) {
    dataset_id <- datasets[d_idx]
    thresholds <- quantile(y[, , dataset_id, setting],
                           probs = probs, na.rm = TRUE)
    validate   <- y[!obs, , dataset_id, setting]

    for (m_idx in seq_along(methods)) {
      method_id <- methods[m_idx]
      f <- build_simstudy_result_file(results_dir, setting, method_id,
                                      dataset_id, mrts_k = k_arg)
      if (!file.exists(f)) next

      e <- new.env()
      load(f, envir = e)
      fit_obj <- if (exists("fit.1", envir = e))
                   get("fit.1", envir = e)
                 else if (exists("fit", envir = e))
                   get("fit", envir = e)
                 else NULL
      if (is.null(fit_obj) || is.null(fit_obj$yp)) next

      pred <- fit_obj$yp
      quant.score[, d_idx, m_idx] <- QuantScore(pred, probs, validate)
      brier.score[, d_idx, m_idx] <- BrierScore(pred, thresholds, validate)

      if (exists("runtime_info", envir = e)) {
        rt <- get("runtime_info", envir = e)
        if (!is.null(rt$elapsed_sec)) {
          elapsed.sec[d_idx, m_idx] <- as.numeric(rt$elapsed_sec)
        }
      }
    }

    cat(sprintf("  dataset %2d / %d\n", dataset_id, length(datasets)))
  }

  save(quant.score, brier.score, elapsed.sec,
       probs, methods, datasets, setting, mrts_k,
       file = out)
  cat(" written:", out, "\n")
}
