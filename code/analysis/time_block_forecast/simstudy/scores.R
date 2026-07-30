#########################################################################
# scores.R - Stage 1 of the time-block forecast post-fit pipeline.
#
# Loads the per-block predictive samples written by run-settings.R
# (results/<setting>-<method>-<dataset>.RData) and computes the
# lead-time scoring protocol of Section 4.3 of ar2_rethink.tex:
#
#   - lead-time curves:   univariate proper scores evaluated separately
#                         at each lead h = 1..H (CRPS, and Brier at a
#                         grid of high thresholds);
#   - joint summary:      energy score and variogram score over the
#                         per-site length-H time vector -- multivariate
#                         proper scores that credit temporal dependence.
#
# Output: output/results/scores<setting><suffix>.RData, or --out=<path>.
#
# Also carried through from the fit files, so that they survive the fit
# deletion performed by the chunked driver: the per-block fit diagnostics
# (fit.diag) and posterior summaries (post.summary), plus the hn_prior
# flag recording which lambda prior produced the fits.
#
# Usage:
#   Rscript scores.R --setting=<id> [--data=<path>]
#                    [--methods=<spec>]   default 1:2
#                    [--datasets=<spec>]  default 1..nsets
#                    [--es_draws=<n>]     default all draws
#                    [--out=<path>]       default output/results/scores<S>.RData
#########################################################################

rm(list = ls())

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0L) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]),
    winslash = "/", mustWork = FALSE)
  if (dir.exists(dirname(script_path))) setwd(dirname(script_path))
}

source("./time_block_helpers.R")
if (!requireNamespace("scoringRules", quietly = TRUE)) {
  stop("scores.R requires the 'scoringRules' package.", call. = FALSE)
}

# ---- CLI parsing ------------------------------------------------------
cli_args <- commandArgs(trailingOnly = TRUE)
parsed <- extract_leading_flags(cli_args,
  c("data", "setting", "methods", "datasets", "es_draws", "out"))
flags <- parsed$values

if (is.null(flags$setting) || !nzchar(flags$setting)) {
  stop("scores.R: --setting=<id> is required.", call. = FALSE)
}

data_path <- resolve_simstudy_data_path(flags$data)
load(data_path)
data_suffix <- derive_data_suffix(data_path)

setting <- as.integer(parse_index_expr(flags$setting, "setting"))
if (length(setting) != 1L || setting < 1L || setting > dim(y)[4]) {
  stop(sprintf("--setting must be a single integer in 1..%d", dim(y)[4]), call. = FALSE)
}

results_dir <- derive_results_dir(data_path, "results")
if (!dir.exists(results_dir)) {
  stop(sprintf("results directory not found: %s (run run-settings.R first)", results_dir),
    call. = FALSE)
}

methods <- if (!is.null(flags$methods) && nzchar(flags$methods)) {
  parse_index_expr(flags$methods, "methods")
} else {
  1:2
}
datasets <- if (!is.null(flags$datasets) && nzchar(flags$datasets)) {
  parse_index_expr(flags$datasets, "datasets")
} else {
  seq_len(as.integer(dim(y)[3]))
}

blocks <- tbf_blocks(block_seams, block_H, nt)
n_blocks <- length(blocks)
H <- block_H

# threshold grid for the lead-time Brier score (high quantiles of the
# full data, matching the Morris study's probs grid).
probs <- c(0.90, 0.92, 0.94, 0.95, 0.96, 0.98, 0.99)

# The energy score is O(m^2) in the number of predictive draws m, and the
# fits carry m = 1e4. At full m it costs ~7 min per (dataset, method) cell,
# which on the chunked driver sits directly in front of the fit deletion.
# Thinning the iteration axis to es_max_draws evenly spaced draws affects
# only es_sample/vs_sample -- CRPS and Brier keep every draw. Same mechanism
# and default as code/analysis/simstudy/scores.R:170-175.
es_max_draws <- if (!is.null(flags$es_draws) && nzchar(flags$es_draws)) {
  as.integer(flags$es_draws)
} else {
  NA_integer_
}
if (!is.na(es_max_draws) && es_max_draws < 100L) {
  stop("--es_draws must be at least 100", call. = FALSE)
}

out_dir <- "output/results"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# An explicit --out lets one process score one dataset into its own chunk
# cache; the fixed default path cannot be written concurrently, and on a
# resumed run it would clobber an already-merged final cache.
out_file <- if (!is.null(flags$out) && nzchar(flags$out)) {
  flags$out
} else {
  file.path(out_dir, sprintf("scores%d%s.RData", setting, data_suffix))
}
if (!dir.exists(dirname(out_file))) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
}

cat(sprintf(
  "scores: setting=%d data=%s results_dir=%s\n  methods=%s datasets=%s blocks=%d H=%d\n  es_draws=%s out=%s\n",
  setting, data_path, results_dir,
  paste(methods, collapse = ","), paste(range(datasets), collapse = ".."),
  n_blocks, H,
  if (is.na(es_max_draws)) "all" else format(es_max_draws), out_file
))

# ---- preallocate ------------------------------------------------------
nsets <- length(datasets)
nmeth <- length(methods)
nprob <- length(probs)

dn_lead <- list(lead = as.character(seq_len(H)),
                dataset = as.character(datasets),
                method = as.character(methods),
                block = as.character(seq_len(n_blocks)))
dn_blk <- list(dataset = as.character(datasets),
               method = as.character(methods),
               block = as.character(seq_len(n_blocks)))

crps.lead <- array(NA_real_, dim = c(H, nsets, nmeth, n_blocks), dimnames = dn_lead)
brier.lead <- array(NA_real_, dim = c(H, nprob, nsets, nmeth, n_blocks),
  dimnames = c(list(lead = as.character(seq_len(H)),
                    quantile = as.character(probs)), dn_lead[-1]))
energy.score <- array(NA_real_, dim = c(nsets, nmeth, n_blocks), dimnames = dn_blk)
vario.score  <- array(NA_real_, dim = c(nsets, nmeth, n_blocks), dimnames = dn_blk)
elapsed_sec  <- array(NA_real_, dim = c(nsets, nmeth), dimnames = dn_blk[1:2])

# Diagnostics and posterior summaries travel with the scores. They are
# numeric arrays with a NAMED "dataset" dimension, not data.frames: that is
# what makes DESKTOP-61SBCCI/merge_score_caches.R bind them along dataset
# instead of demanding they be identical across chunk caches.
diag_stats <- c("lambda", "beta0", "zbar", "z_pred", "z_ratio", "sdz_ratio",
                "phi1.z", "phi2.z", "phi1.tau", "phi2.tau",
                "sd_lead1", "sd_lead_max", "mu_recon", "data_mean",
                "A_zconsist", "Aprime_sdz", "B_truth", "C_spread")
fit.diag <- array(NA_real_, dim = c(length(diag_stats), nsets, nmeth, n_blocks),
  dimnames = c(list(stat = diag_stats), dn_blk))

post_params <- tbf_post_params()
post.summary <- array(NA_real_,
  dim = c(length(TBF_POST_STATS), length(post_params), nsets, nmeth, n_blocks),
  dimnames = c(list(pstat = TBF_POST_STATS, param = post_params), dn_blk))

# Which lambda prior produced these fits; NA if the cells disagree, which
# the chunk gate treats as fatal (an accidentally mixed cache).
hn_prior_seen <- logical(0)

# ---- per-block scoring helpers ---------------------------------------
score_block <- function(yhat, y_val, es_draws = NA_integer_) {
  # yhat: iters x ns x H predictive sample; y_val: ns x H truth.
  ns <- dim(yhat)[2]
  Hh <- dim(yhat)[3]
  iters_n <- dim(yhat)[1]
  es_idx <- if (!is.na(es_draws) && iters_n > es_draws) {
    unique(round(seq(1, iters_n, length.out = es_draws)))
  } else {
    seq_len(iters_n)
  }

  # lead-time CRPS: average crps_sample over sites at each lead.
  crps_h <- numeric(Hh)
  for (h in seq_len(Hh)) {
    # crps_sample is vectorised: y length ns, dat ns x iters.
    crps_h[h] <- mean(scoringRules::crps_sample(
      y = y_val[, h], dat = t(yhat[, , h])), na.rm = TRUE)
  }

  # lead-time Brier at each threshold quantile.
  thr <- quantile(y_val, probs = probs, na.rm = TRUE)
  brier_h <- matrix(NA_real_, Hh, length(probs))
  for (h in seq_len(Hh)) {
    for (q in seq_along(probs)) {
      phat <- colMeans(yhat[, , h] > thr[q])      # ns predicted exceed prob
      ind <- as.numeric(y_val[, h] > thr[q])
      brier_h[h, q] <- mean((ind - phat)^2, na.rm = TRUE)
    }
  }

  # joint-structure scores: per-site length-H time vector.
  es <- vs <- numeric(ns)
  for (i in seq_len(ns)) {
    dat_i <- t(yhat[es_idx, i, ])                 # H x length(es_idx)
    es[i] <- scoringRules::es_sample(y = y_val[i, ], dat = dat_i)
    vs[i] <- scoringRules::vs_sample(y = y_val[i, ], dat = dat_i, p = 0.5)
  }

  list(crps = crps_h, brier = brier_h,
       energy = mean(es, na.rm = TRUE), vario = mean(vs, na.rm = TRUE))
}

# ---- score loop -------------------------------------------------------
for (di in seq_along(datasets)) {
  set <- datasets[di]
  for (mi in seq_along(methods)) {
    method <- methods[mi]
    f <- build_tbf_result_file(results_dir, setting, method, set)
    if (!file.exists(f)) {
      cat("missing: ", f, "\n", sep = "")
      next
    }
    env <- new.env(parent = emptyenv())
    load(f, envir = env)
    fc <- env$forecast
    if (is.null(fc)) { cat("no forecast: ", f, "\n", sep = ""); next }

    for (b in seq_along(fc$blocks)) {
      blk <- fc$blocks[[b]]
      sc <- score_block(blk$yhat, blk$y_val, es_draws = es_max_draws)
      crps.lead[, di, mi, b] <- sc$crps
      brier.lead[, , di, mi, b] <- sc$brier
      energy.score[di, mi, b] <- sc$energy
      vario.score[di, mi, b] <- sc$vario
    }
    if (exists("runtime_info", envir = env, inherits = FALSE)) {
      rt <- env$runtime_info
      if (!is.null(rt$elapsed_sec)) elapsed_sec[di, mi] <- as.numeric(rt$elapsed_sec)
      if (!is.null(rt$control$lambda_positive)) {
        hn_prior_seen <- c(hn_prior_seen, isTRUE(rt$control$lambda_positive))
      }
    }
    # Carry the inline diagnostics / posterior summaries into the cache.
    # Fits written before these existed simply leave the arrays NA.
    if (exists("fit_diag", envir = env, inherits = FALSE)) {
      fd <- env$fit_diag
      for (b in seq_len(nrow(fd))) {
        bi <- as.integer(fd$block[b])
        fit.diag[, di, mi, bi] <- as.numeric(unlist(fd[b, diag_stats]))
      }
    }
    if (exists("post_summary", envir = env, inherits = FALSE)) {
      ps <- env$post_summary
      for (bi in unique(as.integer(ps$block))) {
        sub <- ps[as.integer(ps$block) == bi, , drop = FALSE]
        idx <- match(post_params, sub$param)
        post.summary[, , di, mi, bi] <-
          t(as.matrix(sub[idx, TBF_POST_STATS, drop = FALSE]))
      }
    }
    rm(env, fc)
    cat(sprintf("dataset %d  method %d done\n", set, method))
  }
}

hn_prior <- if (length(unique(hn_prior_seen)) == 1L) unique(hn_prior_seen) else NA

fits_dir <- results_dir
save(crps.lead, brier.lead, energy.score, vario.score, elapsed_sec,
  fit.diag, post.summary, hn_prior, es_max_draws,
  probs, datasets, methods, setting, block_H, block_seams,
  data_path, data_suffix, fits_dir,
  file = out_file)
cat("\nWrote ", out_file, "\n", sep = "")
