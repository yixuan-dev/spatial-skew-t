# brier_bootstrap.R
# ---------------------------------------------------------------------------
# Site-clustered bootstrap confidence intervals and paired tests for held-out
# Brier-score DIFFERENCES between two fitted settings.
#
# WHY THIS EXISTS
#   autoselect.R ranks candidate settings by the pooled Brier score, which
#   collapses all validation site-days into a single number. Those rankings are
#   often decided at the 4th decimal (e.g. Rank1..Rank4 all read 0.030), i.e.
#   inside the noise. This script attaches an uncertainty statement to each
#   pairwise Brier difference.
#
# THE ESTIMAND
#   For two settings A and B we report the PAIRED difference of their pooled
#   Brier scores, Delta = BS_A - BS_B (BS_A is exactly what autoselect prints).
#   Delta < 0 means A is better.
#
# THE UNCERTAINTY
#   The unit of spatial independence is the VALIDATION SITE, not the site-day
#   (days within a site are correlated). We therefore bootstrap by resampling
#   sites with replacement, keeping every day of a resampled site together
#   (a cluster / block bootstrap). Because A and B are scored on the SAME
#   held-out cells, we form the per-site difference FIRST and resample that,
#   so the large shared "how hard is this site" component cancels (pairing).
#
#   Pooled Brier is a ratio  sum_sites SSE / sum_sites ndays , so the bootstrap
#   statistic is a ratio estimator: resample sites, recompute the ratio.
#
# OUTPUT
#   One row per (pair x exceedance level): BS_A, BS_B, Delta, a 95% percentile
#   CI for Delta, a two-sided bootstrap p-value, and n_exceed (how many
#   exceedance events actually stand behind the score at that level -- the
#   honest driver of why high thresholds are unstable).
# ---------------------------------------------------------------------------

rm(list = ls())

# ---- locate script dir (so it runs from anywhere) ----
.this <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), winslash = "/", mustWork = FALSE)) else "."
})
setwd(.this)

# ---- source BrierScoreSite (per-site Brier), preferring the ar2 copy ----
.r_root <- normalizePath(file.path(.this, "../../../R"), winslash = "/", mustWork = FALSE)
source(local({
  ar2 <- file.path(.r_root, "ar2", "auxfunctions.R")
  if (file.exists(ar2)) ar2 else file.path(.r_root, "auxfunctions.R")
}))

# ===========================================================================
# CONFIG  (edit here, or override the split via env var US_ALL_VAL_RESULTS_DIR)
# ===========================================================================
setup_file  <- Sys.getenv("US_ALL_SETUP_FILE",      unset = "us-all-setup-auto-200-200-400.RData")
fits_dir    <- Sys.getenv("US_ALL_VAL_RESULTS_DIR", unset = "fits-200-200")
out_tag     <- Sys.getenv("US_ALL_BOOT_TAG",        unset = "200-200")  # suffix for the CSV
target_probs <- c(0.90, 0.95, 0.98, 0.99, 0.995)
B           <- 2000L          # bootstrap resamples
seed        <- 2024L

# Matched pairs: each isolates ONE extension by holding everything else fixed.
#   c(A, B, label): reports Delta = BS_A - BS_B; A is the "treatment" model.
pairs <- list(
  c(111, 61,  "AR(2) vs AR(1)  | skew-t, K=7, T=50"),   # AR(2) temporal effect
  c(204, 51,  "MRTS(10) vs none | skew-t, K=1, T=0, TS") # MRTS spatial-mean effect
)
# ===========================================================================

# ---- load data split ----
se <- new.env(); load(setup_file, envir = se)
Y  <- get("Y", se); val_sites <- get("split.lst", se)$val
Y_val <- Y[val_sites, , drop = FALSE]              # np x nt observed
np <- nrow(Y_val); nt <- ncol(Y_val)
nday <- rowSums(!is.na(Y_val))                     # non-NA days per site
keep <- nday > 0                                   # drop all-missing sites

thresholds <- quantile(Y, probs = target_probs, na.rm = TRUE)  # same as autoselect

cat(sprintf("Setup           : %s\n", setup_file))
cat(sprintf("Fits dir        : %s\n", fits_dir))
cat(sprintf("Validation sites: %d (%d with >=1 obs day)\n", np, sum(keep)))
cat(sprintf("Days            : %d   |   NA cells: %.1f%%\n", nt, 100 * mean(is.na(Y_val))))
cat("Thresholds (full Y quantiles):\n")
for (k in seq_along(target_probs))
  cat(sprintf("  p=%.3f -> %.3f   (n_exceed=%d)\n",
              target_probs[k], thresholds[k], sum(Y_val > thresholds[k], na.rm = TRUE)))
cat("\n")

# ---- per-site SSE for one setting: SSE_site = (mean brier over its days) * ndays ----
site_sse <- function(sid) {
  fp <- file.path(fits_dir, sprintf("val-%d.RData", sid))
  if (!file.exists(fp)) stop("missing fit: ", fp)
  le <- new.env(); load(fp, envir = le); fit <- get("fit", le)
  yp <- fit$yp                                      # draws x np x nt
  cat(sprintf("  loaded setting %-4d  draws=%d\n", sid, dim(yp)[1]))
  per_site_mean <- BrierScoreSite(preds = yp, thresholds = thresholds, validate = Y_val)
  per_site_mean * nday                              # np x nthresh  (SSE, not mean)
}

# drop any pair whose two fits are not both present in this split
has_fit <- function(sid) file.exists(file.path(fits_dir, sprintf("val-%d.RData", sid)))
keep_pair <- vapply(pairs, function(p) has_fit(as.integer(p[1])) && has_fit(as.integer(p[2])), logical(1))
if (any(!keep_pair)) for (p in pairs[!keep_pair])
  cat(sprintf("  [skip] pair '%s' -- missing fit for %s or %s in %s\n",
              p[3], p[1], p[2], fits_dir))
pairs <- pairs[keep_pair]
if (length(pairs) == 0L) stop("No requested pair has both fits present in ", fits_dir)

need <- unique(unlist(lapply(pairs, function(p) as.integer(p[1:2]))))
cat("Scoring", length(need), "settings (per-site)...\n")
SSE <- setNames(lapply(need, site_sse), as.character(need))
cat("\n")

# ---- the cluster bootstrap on a paired site-difference ----
boot_pair <- function(sidA, sidB, k) {
  sseA <- SSE[[as.character(sidA)]][keep, k]
  sseB <- SSE[[as.character(sidB)]][keep, k]
  d    <- sseA - sseB                               # per-site SSE difference
  w    <- nday[keep]                                # per-site day counts (denominator)
  ratio <- function(idx) sum(d[idx]) / sum(w[idx])  # pooled Brier difference
  point <- ratio(seq_along(d))

  set.seed(seed)
  bootv <- replicate(B, {
    idx <- sample.int(length(d), replace = TRUE)     # resample SITES (clusters)
    ratio(idx)
  })
  ci <- quantile(bootv, c(0.025, 0.975), names = FALSE)
  p  <- min(1, 2 * min(mean(bootv <= 0), mean(bootv >= 0)))  # two-sided percentile p

  data.frame(
    level    = target_probs[k],
    BS_A     = sum(sseA) / sum(w),
    BS_B     = sum(sseB) / sum(w),
    Delta    = point,
    CI_lo    = ci[1], CI_hi = ci[2],
    boot_p   = p,
    n_exceed = sum(Y_val[keep, ] > thresholds[k], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

# ---- run every pair over every level ----
out <- do.call(rbind, lapply(pairs, function(p) {
  A <- as.integer(p[1]); Bx <- as.integer(p[2]); lab <- p[3]
  rows <- do.call(rbind, lapply(seq_along(target_probs), function(k) boot_pair(A, Bx, k)))
  cbind(pair = lab, A = A, B = Bx, rows)
}))

# ---- report ----
fmt <- out
for (col in c("BS_A", "BS_B", "Delta", "CI_lo", "CI_hi"))
  fmt[[col]] <- sprintf("%+.5f", out[[col]])
fmt$boot_p <- sprintf("%.3f", out$boot_p)
fmt$sig    <- ifelse(out$CI_lo > 0 | out$CI_hi < 0, "*", "")   # CI excludes 0

cat("=== Paired Brier differences (Delta = BS_A - BS_B; negative => A better) ===\n")
cat("=== 95% site-clustered bootstrap CI; '*' = CI excludes 0 ===\n\n")
print(fmt[, c("pair", "level", "BS_A", "BS_B", "Delta", "CI_lo", "CI_hi", "boot_p", "n_exceed", "sig")],
      row.names = FALSE)

dir.create("output/us-all-auto/tables", recursive = TRUE, showWarnings = FALSE)
outfile <- sprintf("output/us-all-auto/tables/brier_bootstrap_pairs_%s.csv", out_tag)
write.csv(out, outfile, row.names = FALSE)
cat(sprintf("\nWritten: %s\n", outfile))
