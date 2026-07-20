# crps_bootstrap.R
# ---------------------------------------------------------------------------
# CRPS companion to brier_bootstrap.R. Same site-clustered (block) bootstrap on
# the paired per-site score difference, but the per-site score is the CRPS of
# the posterior predictive distribution (whole-distribution proper score),
# not the single-threshold Brier. CRPS is level-free, so there is one number
# per pair (no exceedance-level loop).
#
#   CRPS(F,y) = integral (F(L) - 1{y<=L})^2 dL,  i.e. Brier integrated over all
#   thresholds -- more information per forecast, hence more power than Brier at
#   any single L (see the time-block result).
#
# Unit of independence = validation SITE (all its days kept together). Estimand
# = paired difference of pooled CRPS, Delta = CRPS_A - CRPS_B (negative => A
# better). Pooled CRPS = sum_site(sum_days crps) / sum_site(ndays), a ratio, so
# the bootstrap statistic is the same ratio recomputed on resampled sites.
#
# Writes: output/us-all-auto/tables/crps_bootstrap_pairs_<split>.csv
# ---------------------------------------------------------------------------

rm(list = ls())
.this <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), winslash = "/", mustWork = FALSE)) else "."
})
setwd(.this)
stopifnot(requireNamespace("scoringRules", quietly = TRUE))

# ---- CONFIG (mirrors brier_bootstrap.R) ----
setup_file <- Sys.getenv("US_ALL_SETUP_FILE",      unset = "us-all-setup-auto-200-200-400.RData")
fits_dir   <- Sys.getenv("US_ALL_VAL_RESULTS_DIR", unset = "fits-200-200")
out_tag    <- Sys.getenv("US_ALL_BOOT_TAG",        unset = "200-200")
B          <- 2000L
seed       <- 2024L
pairs <- list(
  c(111, 61,  "AR(2) vs AR(1)  | skew-t, K=7, T=50"),
  c(204, 51,  "MRTS(10) vs none | skew-t, K=1, T=0, TS")
)

# ---- data split ----
se <- new.env(); load(setup_file, envir = se)
Y  <- get("Y", se); val_sites <- get("split.lst", se)$val
Y_val <- Y[val_sites, , drop = FALSE]
np <- nrow(Y_val); nt <- ncol(Y_val)
nday <- rowSums(!is.na(Y_val)); keep <- nday > 0

cat(sprintf("Setup: %s | fits: %s | val sites: %d | days: %d\n",
            setup_file, fits_dir, np, nt))

# drop pairs missing a fit in this split
has_fit <- function(sid) file.exists(file.path(fits_dir, sprintf("val-%d.RData", sid)))
keep_pair <- vapply(pairs, function(p) has_fit(as.integer(p[1])) && has_fit(as.integer(p[2])), logical(1))
for (p in pairs[!keep_pair]) cat(sprintf("  [skip] '%s' -- missing fit in %s\n", p[3], fits_dir))
pairs <- pairs[keep_pair]
if (!length(pairs)) stop("no pair has both fits present in ", fits_dir)

# ---- per-site CRPS sum: sum over a site's days of crps_sample(y, draws) ----
site_scrps <- function(sid) {
  le <- new.env(); load(file.path(fits_dir, sprintf("val-%d.RData", sid)), envir = le)
  yp <- get("fit", le)$yp                          # draws x np x nt
  d <- dim(yp)
  dat <- matrix(aperm(yp, c(2, 3, 1)), nrow = d[2] * d[3], ncol = d[1])  # (np*nt) x draws
  yv <- as.vector(Y_val)                           # site-fastest, matches dat rows
  crps_cell <- rep(NA_real_, length(yv))           # crps_sample rejects NA obs
  ok <- !is.na(yv)
  crps_cell[ok] <- scoringRules::crps_sample(y = yv[ok], dat = dat[ok, , drop = FALSE])
  crps_mat <- matrix(crps_cell, np, nt)            # site x day
  cat(sprintf("  scored setting %-4d draws=%d\n", sid, d[1]))
  rowSums(crps_mat, na.rm = TRUE)                  # per-site CRPS sum (length np)
}

need <- unique(unlist(lapply(pairs, function(p) as.integer(p[1:2]))))
SCR <- setNames(lapply(need, site_scrps), as.character(need))

# ---- site-clustered bootstrap on the paired difference ----
boot_pair <- function(sidA, sidB) {
  a <- SCR[[as.character(sidA)]][keep]; b <- SCR[[as.character(sidB)]][keep]
  d <- a - b; w <- nday[keep]
  ratio <- function(idx) sum(d[idx]) / sum(w[idx])
  point <- ratio(seq_along(d))
  set.seed(seed)
  bootv <- replicate(B, ratio(sample.int(length(d), replace = TRUE)))
  ci <- quantile(bootv, c(0.025, 0.975), names = FALSE)
  p  <- min(1, 2 * min(mean(bootv <= 0), mean(bootv >= 0)))
  data.frame(CRPS_A = sum(a) / sum(w), CRPS_B = sum(b) / sum(w),
             Delta = point, CI_lo = ci[1], CI_hi = ci[2], boot_p = p)
}

out <- do.call(rbind, lapply(pairs, function(p) {
  A <- as.integer(p[1]); Bx <- as.integer(p[2])
  cbind(pair = p[3], A = A, B = Bx, boot_pair(A, Bx))
}))

fmt <- out
for (c0 in c("CRPS_A", "CRPS_B", "Delta", "CI_lo", "CI_hi")) fmt[[c0]] <- sprintf("%+.4f", out[[c0]])
fmt$boot_p <- sprintf("%.3f", out$boot_p)
fmt$sig <- ifelse(out$CI_lo > 0 | out$CI_hi < 0, "*", "")
cat("\n=== Paired CRPS differences (Delta = CRPS_A - CRPS_B; negative => A better) ===\n")
cat("=== 95% site-clustered bootstrap CI; '*' = CI excludes 0 ===\n\n")
print(fmt[, c("pair", "CRPS_A", "CRPS_B", "Delta", "CI_lo", "CI_hi", "boot_p", "sig")], row.names = FALSE)

dir.create("output/us-all-auto/tables", recursive = TRUE, showWarnings = FALSE)
outfile <- sprintf("output/us-all-auto/tables/crps_bootstrap_pairs_%s.csv", out_tag)
write.csv(out, outfile, row.names = FALSE)
cat(sprintf("\nWritten: %s\n", outfile))
