#########################################################################
# chunk_sanity_tbf.R -- validity gate for one time-block score cache. The
# driver deletes a dataset's ~2.4 GB of fits only after this passes, so
# every assertion here exists to stop an unrecoverable deletion.
#
# Usage (from code/analysis/time_block_forecast/simstudy):
#   Rscript DESKTOP-61SBCCI/chunk_sanity_tbf.R <cache.RData> <setting> <datasets> <methods> [es_draws]
#   e.g. Rscript DESKTOP-61SBCCI/chunk_sanity_tbf.R \
#          output/results/scores5_d3.RData 5 3 "c(1,2,4)" 1000
#
# The A/A'/B/C guard flags are no longer recorded (2026-08-03): the
# lambda ~ HN(0, 20) prior is the protection against the reflected ridge,
# and fit.diag carries the 14 numeric summaries only.
#########################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4 || length(args) > 6) {
  stop(paste("usage: chunk_sanity_tbf.R <cache.RData> <setting>",
             "<datasets> <methods> [es_draws] [min_elapsed_sec]"), call. = FALSE)
}
cache_file <- args[1]
setting_expect <- as.integer(args[2])
datasets_expect <- sort(as.integer(eval(parse(text = args[3]))))
methods_expect <- as.integer(eval(parse(text = args[4])))
es_expect <- if (length(args) >= 5) as.integer(args[5]) else NA_integer_
# Production cells take 1.5-2 h; anything far under that means mcmc() gave
# up early. Lowered only by the smoke test, which fits at toy iteration
# counts -- the driver always leaves it at the default.
min_elapsed <- if (length(args) >= 6) as.numeric(args[6]) else 1800

fail <- function(...) {
  cat("[SANITY FAILED] ", sprintf(...), "\n", sep = "")
  quit(status = 1)
}

if (!file.exists(cache_file)) fail("cache not found: %s", cache_file)
e <- new.env()
load(cache_file, envir = e)

need <- c("crps.lead", "brier.lead", "brier.lead.blockq", "pexceed.mean",
          "brier.thresholds", "exceed.rate.lead", "brier_threshold_basis",
          "energy.score", "vario.score",
          "elapsed_sec", "fit.diag", "post.summary", "hn_prior",
          "es_max_draws", "probs", "datasets", "methods", "setting",
          "block_H", "block_seams", "data_path", "data_suffix", "fits_dir")
for (nm in need) {
  if (!exists(nm, envir = e, inherits = FALSE)) fail("object `%s` missing", nm)
}

nD <- length(datasets_expect)
nM <- length(methods_expect)
H <- 15L
NB <- 5L
NPROB <- 11L

# ---- identity of the cache --------------------------------------------
if (!identical(as.integer(e$setting), setting_expect)) {
  fail("setting is %s, expected %d", e$setting, setting_expect)
}
# scores.R defaults --datasets to ALL datasets in simdata (1:50) and fills
# the absent ones with NA. Gating such a cache would delete the fits of a
# chunk that was, for the most part, never scored.
if (!identical(sort(as.integer(e$datasets)), datasets_expect)) {
  fail("datasets are [%s], expected [%s] -- was --datasets passed?",
       paste(e$datasets, collapse = ","), paste(datasets_expect, collapse = ","))
}
if (!identical(dimnames(e$crps.lead)[["dataset"]],
               as.character(sort(as.integer(e$datasets))))) {
  fail("crps.lead dataset dimnames disagree with the `datasets` vector")
}
# scores.R defaults --methods to 1:2. Missing method 4 here would mean its
# fits are deleted having never been scored: ~40 core-hours, unrecoverable.
if (!identical(as.integer(e$methods), methods_expect)) {
  fail("methods are [%s], expected [%s] -- was --methods passed?",
       paste(e$methods, collapse = ","), paste(methods_expect, collapse = ","))
}

# ---- the run's premise -------------------------------------------------
if (!isTRUE(e$hn_prior)) {
  fail("hn_prior is %s, expected TRUE -- was --hn passed to run-settings.R?",
       format(e$hn_prior))
}
if (!identical(as.integer(e$es_max_draws), es_expect)) {
  fail("es_max_draws is %s, expected %s",
       format(e$es_max_draws), format(es_expect))
}

# ---- geometry ----------------------------------------------------------
if (!identical(as.integer(e$block_H), H)) fail("block_H is %s, expected %d", e$block_H, H)
if (!identical(as.integer(e$block_seams), c(50L, 80L, 110L, 140L, 170L))) {
  fail("block_seams are [%s], expected 50,80,110,140,170",
       paste(e$block_seams, collapse = ","))
}
if (!isTRUE(all.equal(as.numeric(e$probs),
                      c(0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97,
                        0.98, 0.99, 0.995)))) {
  fail("probs grid is [%s], expected the 11-value high-quantile grid of scores.R",
       paste(e$probs, collapse = ","))
}
# The probs grid is identical under both threshold rules, so it cannot
# detect a rule change; the basis scalar is what does.
if (!identical(e$brier_threshold_basis, "full_series")) {
  fail(paste("brier_threshold_basis is '%s', expected 'full_series' --",
             "this cache was written by the pre-2026-07-31 scorer"),
       format(e$brier_threshold_basis))
}

chk_dim <- function(nm, expect_dim, expect_names) {
  o <- get(nm, envir = e)
  if (!identical(dim(o), expect_dim)) {
    fail("`%s` dim is [%s], expected [%s]", nm,
         paste(dim(o), collapse = "x"), paste(expect_dim, collapse = "x"))
  }
  if (!identical(names(dimnames(o)), expect_names)) {
    fail("`%s` dimnames names are [%s], expected [%s]", nm,
         paste(names(dimnames(o)), collapse = ","),
         paste(expect_names, collapse = ","))
  }
}
chk_dim("crps.lead", c(H, nD, nM, NB), c("lead", "dataset", "method", "block"))
chk_dim("brier.lead", c(H, NPROB, nD, nM, NB),
        c("lead", "quantile", "dataset", "method", "block"))
chk_dim("brier.lead.blockq", c(H, NPROB, nD, nM, NB),
        c("lead", "quantile", "dataset", "method", "block"))
chk_dim("pexceed.mean", c(H, NPROB, nD, nM, NB),
        c("lead", "quantile", "dataset", "method", "block"))
chk_dim("brier.thresholds", c(NPROB, nD), c("quantile", "dataset"))
chk_dim("exceed.rate.lead", c(H, NPROB, nD, NB),
        c("lead", "quantile", "dataset", "block"))
# brier.thresholds' dataset axis is the one the merge binds along.
if (!identical(dimnames(e$brier.thresholds)[["dataset"]],
               as.character(sort(as.integer(e$datasets))))) {
  fail("brier.thresholds dataset dimnames disagree with the `datasets` vector")
}
chk_dim("energy.score", c(nD, nM, NB), c("dataset", "method", "block"))
chk_dim("vario.score", c(nD, nM, NB), c("dataset", "method", "block"))
chk_dim("elapsed_sec", c(nD, nM), c("dataset", "method"))
chk_dim("fit.diag", c(14L, nD, nM, NB), c("stat", "dataset", "method", "block"))
chk_dim("post.summary", c(5L, 15L, nD, nM, NB),
        c("pstat", "param", "dataset", "method", "block"))

# ---- completeness: no NA is legitimate in this study -------------------
# scores.R silently `next`s past a missing or forecast-less fit, leaving
# NA behind. This is the assertion that catches it.
for (nm in c("crps.lead", "energy.score", "vario.score", "elapsed_sec")) {
  o <- get(nm, envir = e)
  if (anyNA(o)) fail("`%s` has %d NA cells -- a fit was not scored", nm, sum(is.na(o)))
  if (!all(is.finite(o))) fail("`%s` has non-finite cells", nm)
  if (any(o <= 0)) fail("`%s` has non-positive cells", nm)
}

# The Brier arrays are mean((ind - phat)^2) over sites: squared probability
# errors on [0, 1], NOT strictly positive scores. An exact 0 is legitimate
# (measured on setting 5 dataset 2 method 1: 5 of 525 cells under the
# full-series rule, 3 of 525 under the block rule). The old
# `mean(brier.lead == 0) > 0.5` heuristic would not have fired under either
# rule; it is replaced by the threshold-provenance recomputation below,
# which asserts what the score actually IS rather than guessing from shape.
for (nm in c("brier.lead", "brier.lead.blockq", "pexceed.mean")) {
  o <- get(nm, envir = e)
  if (anyNA(o)) fail("`%s` has %d NA cells -- a fit was not scored", nm, sum(is.na(o)))
  if (!all(is.finite(o))) fail("`%s` has non-finite cells", nm)
  if (any(o < 0 | o > 1)) fail("`%s` is outside [0, 1] (min %g, max %g)",
                               nm, min(o), max(o))
}
# Wiring check: the two rules cannot produce identical numbers unless both
# arms of score_block() were handed the same thresholds.
if (identical(e$brier.lead, e$brier.lead.blockq)) {
  fail("`brier.lead` and `brier.lead.blockq` are identical -- both arms got the same thresholds")
}
# These two derive from the data, not from a fit; any NA means the dataset
# loop never reached them. (They are NOT a was-this-scored proxy.)
if (anyNA(e$brier.thresholds)) fail("`brier.thresholds` has NA cells")
if (anyNA(e$exceed.rate.lead)) fail("`exceed.rate.lead` has NA cells")

# ---- diagnostics / posterior summaries must have arrived ---------------
for (st in c("lambda", "beta0", "sd_lead_max")) {
  if (anyNA(e$fit.diag[st, , , , drop = FALSE])) {
    fail("`fit.diag` row `%s` has NA -- inline diagnostics did not reach the cache", st)
  }
}
if (anyNA(e$post.summary["mean", "beta0", , , , drop = FALSE])) {
  fail("`post.summary` has NA for beta0 -- posterior summaries did not reach the cache")
}
# (The per-cell output/diag mirrors are gone since 2026-08-03; the two
# assertions above ARE the durability check -- the cache itself carries
# the diagnostics and posterior summaries past the fit deletion.)

# ---- provenance of the data --------------------------------------------
if (!identical(e$data_suffix, "")) fail("data_suffix is '%s', expected ''", e$data_suffix)
if (!identical(e$fits_dir, "results")) fail("fits_dir is '%s', expected 'results'", e$fits_dir)
if (!identical(basename(e$data_path), "simdata.RData")) {
  fail("data_path is '%s', expected .../simdata.RData", e$data_path)
}

# ---- provenance of the Brier thresholds --------------------------------
# This gate authorises an irreversible 2.4 GB fit deletion, and the ONE
# thing a shape check cannot see is whether the thresholds really are the
# full-series ones. So recompute them, and the base rates they induce,
# from the data file the cache names, and require bit-exact agreement.
# Costs ~0.6 s; everything (probs, seams, H, path) is read from the cache,
# nothing hardcoded, so a toy smoke cache exercises this too. Fail closed
# on a missing data file -- a gate that cannot verify must not pass.
if (!file.exists(e$data_path)) {
  fail("data_path does not exist: %s -- cannot verify threshold provenance",
       e$data_path)
}
d <- new.env(parent = emptyenv())
load(e$data_path, envir = d)
if (!exists("y", envir = d, inherits = FALSE)) fail("`y` not found in %s", e$data_path)
pr <- as.numeric(e$probs)
ds_dn <- dimnames(e$brier.thresholds)[["dataset"]]
for (di in seq_along(datasets_expect)) {
  set <- datasets_expect[di]
  if (!identical(ds_dn[di], as.character(set))) {
    fail("brier.thresholds dataset dimname [%d] is '%s', expected '%d'",
         di, ds_dn[di], set)
  }
  thr <- quantile(d$y[, , set, setting_expect], probs = pr,
                  na.rm = TRUE, names = FALSE)
  got <- as.numeric(e$brier.thresholds[, di])
  if (!isTRUE(all.equal(got, thr, tolerance = 0))) {
    fail("dataset %d thresholds are NOT the full-series quantiles (max abs diff %g)",
         set, max(abs(got - thr)))
  }
  if (is.unsorted(thr, strictly = TRUE)) {
    fail("dataset %d thresholds are not strictly increasing", set)
  }
  for (b in seq_len(NB)) {
    tt <- as.integer(e$block_seams[b]) + seq_len(as.integer(e$block_H))
    yv <- d$y[, tt, set, setting_expect]
    for (q in seq_along(pr)) {
      er <- colMeans(yv > thr[q], na.rm = TRUE)
      if (!isTRUE(all.equal(as.numeric(e$exceed.rate.lead[, q, di, b]),
                            as.numeric(er), tolerance = 0))) {
        fail("dataset %d block %d prob %g: exceed.rate.lead does not match the data",
             set, b, pr[q])
      }
    }
  }
}
rm(d); invisible(gc(FALSE))

# ---- method identity (partial cover for the ar1 forecast wiring) -------
mdn <- dimnames(e$fit.diag)[["method"]]
has_phi1 <- function(m) all(!is.na(e$fit.diag["phi1.z", , as.character(m), ]))
has_phi2 <- function(m) all(!is.na(e$fit.diag["phi2.z", , as.character(m), ]))
if ("1" %in% mdn && (has_phi1(1) || has_phi2(1))) {
  fail("method 1 carries phi chains -- it should be i.i.d. in time")
}
if ("2" %in% mdn && !(has_phi1(2) && has_phi2(2))) {
  fail("method 2 is missing a phi lag -- it should be AR(2)")
}
if ("4" %in% mdn && !(has_phi1(4) && !has_phi2(4))) {
  fail("method 4 phi structure is not single-lag -- it should be AR(1)")
}

# ---- runtime plausibility (last, so a toy smoke cache can exercise -----
# ---- every assertion above before tripping on this one) ----------------
# A production cell takes 1.5-2 h. Far under that means mcmc() gave up
# early -- options(warn = 2) makes any warning fatal inside the sampler,
# and the cell would still have saved whatever it had.
if (any(e$elapsed_sec < min_elapsed)) {
  fail("elapsed_sec has cells under %.0f s (min %.0f s) -- did a fit abort?",
       min_elapsed, min(e$elapsed_sec))
}
if (any(e$elapsed_sec > 6 * 3600)) {
  fail("elapsed_sec has cells over 6 h (max %.0f s)", max(e$elapsed_sec))
}

# ---- descriptive only: reported, never fatal ---------------------------
lam_neg <- sum(e$fit.diag["lambda", , , ] < 0, na.rm = TRUE)

# Brier coverage ledger under the full-series rule. Both resolutions are
# printed so nobody quotes the flattering block-level number: brier.lead's
# unit is the lead, and at high thresholds most leads of a quiet block have
# no exceedance at all. The level check is the diagnostic that matters --
# predicted far below realised means the predictive never reaches the
# threshold and the Brier there is climatology, not forecast skill.
er <- e$exceed.rate.lead
q95i <- which(abs(as.numeric(e$probs) - 0.95) < 1e-12)
er_blk <- apply(er, c(2, 3, 4), mean)              # (quantile, dataset, block)
pex95 <- apply(e$pexceed.mean[, q95i, , , , drop = FALSE], 4, mean)
obs95 <- mean(er[, q95i, , ])

cat(sprintf(paste("[SANITY OK] %s: setting %d, datasets %s, methods [%s],",
                  "%d blocks, HN prior, es_draws %s\n"),
            cache_file, setting_expect,
            paste(range(datasets_expect), collapse = ".."),
            paste(methods_expect, collapse = ","), NB, format(es_expect)))
cat(sprintf("            diagnostics (descriptive): lambda<0 %d of %d block-cells\n",
            lam_neg, nD * nM * NB))
cat(sprintf("            brier thresholds (full series): q90 %.2f .. q99 %.2f\n",
            min(e$brier.thresholds[1, ]), max(e$brier.thresholds[NPROB, ])))
cat(sprintf("            exceedance ledger: base rate 0 in %.1f%% of lead-cells, %.1f%% of block-cells\n",
            100 * mean(er == 0), 100 * mean(er_blk == 0)))
cat(sprintf("            level check q95: realised %.4f vs predicted [%s] by method [%s]\n",
            obs95, paste(sprintf("%.4f", pex95), collapse = ", "),
            paste(dimnames(e$pexceed.mean)[["method"]], collapse = ", ")))
