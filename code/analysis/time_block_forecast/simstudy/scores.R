#########################################################################
# scores.R - Stage 1 of the time-block forecast post-fit pipeline.
#
# Loads the per-block predictive samples written by run-settings.R
# (results/<setting>-<method>-<dataset>.RData) and computes the
# lead-time scoring protocol of Section 4.3 of ar2_rethink.tex:
#
#   - lead-time curves:   univariate proper scores evaluated separately
#                         at each lead h = 1..H (CRPS, and Brier at a
#                         grid of high FULL-SERIES threshold quantiles;
#                         the block-quantile rule is retained as a
#                         diagnostic array, see the probs block below);
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
#                    [--hn | --prior=hn|n]  which lambda-prior arm to score;
#                                           selects results_<tag>/ and the
#                                           cache name. Default n.
#                    [--methods=<spec>]   default 1:2
#                    [--datasets=<spec>]  default 1..nsets
#                    [--es_draws=<n>]     default all draws
#                    [--out=<path>]       default
#                                         output/results/scores<S>_<tag>.RData
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
# extract_leading_flags() demands a value after a bare "--flag", so pull the
# valueless form of --hn out first; --hn=TRUE and --prior=hn still parse.
hn_bare <- any(cli_args == "--hn")
cli_args <- cli_args[cli_args != "--hn"]
parsed <- extract_leading_flags(cli_args,
  c("data", "setting", "methods", "datasets", "es_draws", "out", "hn", "prior"))
flags <- parsed$values

# Which arm of the lambda-prior toggle to score. This selects the results
# directory, so scoring the wrong arm is a missing-file error rather than a
# silently mixed cache.
lambda_positive <- hn_bare || tbf_parse_prior(cli_args, flags)
prior_tag <- tbf_prior_tag(lambda_positive)

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

results_dir <- derive_results_dir(data_path, "results", prior_tag)
if (!dir.exists(results_dir)) {
  stop(sprintf(
    "results directory not found: %s (run run-settings.R%s first)",
    results_dir, if (lambda_positive) " --hn" else ""), call. = FALSE)
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

# Threshold grid for the lead-time Brier score (the Morris study's probs
# grid). Two threshold RULES are scored at every cell:
#   brier.lead        thresholds = quantile(y[, , set, setting], probs) --
#                     one vector per (dataset, setting) from the FULL
#                     series, identical to code/analysis/simstudy/
#                     scores.R:333. This is the PRIMARY rule.
#   brier.lead.blockq thresholds = quantile(y_val, probs) -- the block's
#                     own held-out window. Retained for continuity with
#                     caches written before 2026-07-31 ONLY: it is
#                     look-ahead, and it pins the exceedance base rate at
#                     exactly 1-p in every block, so it cannot detect a
#                     level error (the near-unit-root excursion of
#                     setting 5 scores the same as a perfect forecast).
#                     Diagnostic; never promote it to a table.
probs <- c(0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995)

# Provenance tripwire. This is NOT a dataset-array, so
# DESKTOP-61SBCCI/merge_score_caches.R:117-126 requires it to be identical
# across every chunk cache -- which makes merging an old-rule cache with a
# new-rule one a hard error rather than a silent average of two
# definitions of "the Brier score".
brier_threshold_basis <- "full_series"

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
  tbf_score_cache_file(setting, data_suffix, prior_tag, dir = out_dir)
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

# Both Brier rules share this shape; as.character(probs) (giving "0.9",
# never quantile()'s "90%") keeps the quantile dimnames byte-identical
# across chunk caches, which bind_dataset_arrays() requires.
dn_brier <- c(list(lead = as.character(seq_len(H)),
                   quantile = as.character(probs)), dn_lead[-1])
brier.lead <- array(NA_real_, dim = c(H, nprob, nsets, nmeth, n_blocks),
  dimnames = dn_brier)
brier.lead.blockq <- array(NA_real_, dim = c(H, nprob, nsets, nmeth, n_blocks),
  dimnames = dn_brier)
# Site-mean predicted exceedance probability under the full-series rule.
# Together with exceed.rate.lead it is the level-error ledger: predicted
# far below realised means the predictive never reaches the threshold and
# the Brier there has degenerated onto climatology.
pexceed.mean <- array(NA_real_, dim = c(H, nprob, nsets, nmeth, n_blocks),
  dimnames = dn_brier)

# The full-series thresholds themselves, and the realised exceedance base
# rate they induce. Both are functions of the DATA only -- populated even
# for a (dataset, block) whose fit is missing, so do not read anyNA() on
# them as "was this dataset scored". Both carry a named "dataset"
# dimension so merge_score_caches.R binds them along dataset instead of
# demanding byte-identity.
brier.thresholds <- array(NA_real_, dim = c(nprob, nsets),
  dimnames = list(quantile = as.character(probs),
                  dataset = as.character(datasets)))
exceed.rate.lead <- array(NA_real_, dim = c(H, nprob, nsets, n_blocks),
  dimnames = list(lead = as.character(seq_len(H)),
                  quantile = as.character(probs),
                  dataset = as.character(datasets),
                  block = as.character(seq_len(n_blocks))))

energy.score <- array(NA_real_, dim = c(nsets, nmeth, n_blocks), dimnames = dn_blk)
vario.score  <- array(NA_real_, dim = c(nsets, nmeth, n_blocks), dimnames = dn_blk)
elapsed_sec  <- array(NA_real_, dim = c(nsets, nmeth), dimnames = dn_blk[1:2])

# Diagnostics and posterior summaries travel with the scores. They are
# numeric arrays with a NAMED "dataset" dimension, not data.frames: that is
# what makes DESKTOP-61SBCCI/merge_score_caches.R bind them along dataset
# instead of demanding they be identical across chunk caches.
# The A/A'/B/C guard flags are gone since 2026-08-03: the HN prior removed
# the ridge they guarded, and run-settings.R stopped recording them.
diag_stats <- c("lambda", "beta0", "zbar", "z_pred", "z_ratio", "sdz_ratio",
                "phi1.z", "phi2.z", "phi1.tau", "phi2.tau",
                "sd_lead1", "sd_lead_max", "mu_recon", "data_mean")
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
score_block <- function(yhat, y_val, thr_full, es_draws = NA_integer_) {
  # yhat: iters x ns x H predictive sample; y_val: ns x H truth;
  # thr_full: the full-series thresholds for this dataset (mandatory and
  # positional BEFORE es_draws, so an un-updated call site errors instead
  # of silently recycling es_draws into the threshold slot).
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

  # lead-time Brier at each threshold quantile, under BOTH rules.
  #   thr_full : full-series quantiles of y[, , set, setting]  (PRIMARY)
  #   thr_blk  : quantiles of this block's own held-out window (diagnostic)
  thr_blk <- quantile(y_val, probs = probs, na.rm = TRUE, names = FALSE)
  brier_h  <- matrix(NA_real_, Hh, length(probs))
  brier_hb <- matrix(NA_real_, Hh, length(probs))
  pex_h    <- matrix(NA_real_, Hh, length(probs))
  for (h in seq_len(Hh)) {
    s <- yhat[, , h]                              # iters x ns, sliced once
    for (q in seq_along(probs)) {
      phat <- colMeans(s > thr_full[q])           # ns predicted exceed prob
      ind <- as.numeric(y_val[, h] > thr_full[q])
      brier_h[h, q] <- mean((ind - phat)^2, na.rm = TRUE)
      pex_h[h, q] <- mean(phat, na.rm = TRUE)
      phb <- colMeans(s > thr_blk[q])
      inb <- as.numeric(y_val[, h] > thr_blk[q])
      brier_hb[h, q] <- mean((inb - phb)^2, na.rm = TRUE)
    }
  }

  # joint-structure scores: per-site length-H time vector.
  es <- vs <- numeric(ns)
  for (i in seq_len(ns)) {
    dat_i <- t(yhat[es_idx, i, ])                 # H x length(es_idx)
    es[i] <- scoringRules::es_sample(y = y_val[i, ], dat = dat_i)
    vs[i] <- scoringRules::vs_sample(y = y_val[i, ], dat = dat_i, p = 0.5)
  }

  list(crps = crps_h, brier = brier_h, brier_blockq = brier_hb,
       pexceed = pex_h,
       energy = mean(es, na.rm = TRUE), vario = mean(vs, na.rm = TRUE))
}

# ---- score loop -------------------------------------------------------
for (di in seq_along(datasets)) {
  set <- datasets[di]

  # Full-series thresholds for this dataset -- the sibling study's rule
  # (code/analysis/simstudy/scores.R:333) verbatim. `set` is the dataset
  # ID, not `di`. names = FALSE keeps quantile()'s "90%" labels out of the
  # cache; the array dimnames come from as.character(probs) only.
  thr_full <- quantile(y[, , set, setting], probs = probs,
                       na.rm = TRUE, names = FALSE)
  brier.thresholds[, di] <- thr_full

  # Descriptive coverage ledger, computed from the DATA (not a fit), so it
  # is complete even where a fit is missing: the realised exceedance base
  # rate at every (lead, quantile, block) under the full-series rule.
  for (b in seq_along(blocks)) {
    yv_b <- y[, blocks[[b]]$test_times, set, setting]     # ns x H
    for (q in seq_along(probs)) {
      exceed.rate.lead[, q, di, b] <- colMeans(yv_b > thr_full[q], na.rm = TRUE)
    }
  }

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
      # Tie the fit file to the data the thresholds came from. The ledger
      # above indexes y by blocks[[b]]$test_times while the score uses the
      # fit's own blk$y_val; if those ever diverge, or the fit was run on
      # a different data file, the full-series thresholds would produce
      # plausible-looking but wrong scores. Assert, don't assume.
      if (!identical(as.integer(blk$test_times),
                     as.integer(blocks[[b]]$test_times))) {
        stop(sprintf("%s block %d: test_times disagree with tbf_blocks()", f, b),
             call. = FALSE)
      }
      if (b == 1L && !isTRUE(all.equal(as.numeric(blk$y_val),
            as.numeric(y[, blk$test_times, set, setting]), tolerance = 0))) {
        stop(sprintf("%s: y_val does not match %s -- wrong data file?",
                     f, data_path), call. = FALSE)
      }
      sc <- score_block(blk$yhat, blk$y_val, thr_full, es_draws = es_max_draws)
      crps.lead[, di, mi, b] <- sc$crps
      brier.lead[, , di, mi, b] <- sc$brier
      brier.lead.blockq[, , di, mi, b] <- sc$brier_blockq
      pexceed.mean[, , di, mi, b] <- sc$pexceed
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

# The fits found in results_<tag>/ must actually BE that arm. The tagged
# path makes a collision impossible going forward, but a directory
# populated before tagging -- or by hand -- could still disagree with the
# flag, and that would mislabel the whole cache. Assert rather than trust.
if (length(hn_prior_seen) > 0L) {
  if (is.na(hn_prior)) {
    stop(sprintf(
      "%s holds a MIX of lambda priors (%d HN, %d N fits) -- the arms must not share a directory",
      results_dir, sum(hn_prior_seen), sum(!hn_prior_seen)), call. = FALSE)
  }
  if (!identical(isTRUE(hn_prior), isTRUE(lambda_positive))) {
    stop(sprintf(
      "%s holds %s fits but the run asked for the %s arm (--prior=%s)",
      results_dir, if (isTRUE(hn_prior)) "HN" else "N",
      if (lambda_positive) "HN" else "N", prior_tag), call. = FALSE)
  }
}

fits_dir <- results_dir
prior_arm <- prior_tag
save(crps.lead, brier.lead, brier.lead.blockq, pexceed.mean,
  brier.thresholds, exceed.rate.lead, brier_threshold_basis,
  energy.score, vario.score, elapsed_sec,
  fit.diag, post.summary, hn_prior, prior_arm, es_max_draws,
  probs, datasets, methods, setting, block_H, block_seams,
  data_path, data_suffix, fits_dir,
  file = out_file)
cat("\nWrote ", out_file, "\n", sep = "")
