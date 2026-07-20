# twcrps_bootstrap.R
# ---------------------------------------------------------------------------
# Threshold-weighted CRPS over the TAIL BAND, via the exact identity
#   twCRPS(F,y) = integral_band ( F(z) - 1{y<=z} )^2 dz = integral_band BS(z) dz,
# i.e. the per-cell Brier score integrated over a band of thresholds z. Here the
# band is [q(0.90), q(0.995)] with uniform weight w(z)=1 (a tail-restricted
# CRPS), evaluated on a fine z-grid and integrated by the trapezoid rule.
#
# This sits BETWEEN single-threshold Brier (one point of the integrand -> high
# variance, low power in the far tail) and plain CRPS (whole distribution ->
# powerful but bulk-dominated). It answers: is plain CRPS's "AR(2) worse for
# ozone interpolation" a TAIL effect (goal-relevant) or a bulk artifact?
#
# Same site-clustered bootstrap as brier_bootstrap.R / crps_bootstrap.R:
# per-site score, resample validation SITES, paired difference, ratio estimator.
#
# Writes: output/us-all-auto/tables/twcrps_bootstrap_pairs_<split>.csv
# ---------------------------------------------------------------------------

rm(list = ls())
.this <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), winslash = "/", mustWork = FALSE)) else "."
})
setwd(.this)

# ---- CONFIG ----
setup_file <- Sys.getenv("US_ALL_SETUP_FILE",      unset = "us-all-setup-auto-200-200-400.RData")
fits_dir   <- Sys.getenv("US_ALL_VAL_RESULTS_DIR", unset = "fits-200-200")
out_tag    <- Sys.getenv("US_ALL_BOOT_TAG",        unset = "200-200")
band_probs <- c(0.90, 0.995)     # tail band (in quantile units)
G          <- 41L                # z-grid resolution over the band
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
yv <- as.vector(Y_val)                          # site-fastest

band <- quantile(Y, probs = band_probs, na.rm = TRUE)
zgrid <- seq(band[1], band[2], length.out = G)  # ppb thresholds across the tail band
cat(sprintf("Setup: %s | fits: %s | sites: %d | band [%.1f, %.1f] ppb, G=%d\n",
            setup_file, fits_dir, np, band[1], band[2], G))

has_fit <- function(sid) file.exists(file.path(fits_dir, sprintf("val-%d.RData", sid)))
keep_pair <- vapply(pairs, function(p) has_fit(as.integer(p[1])) && has_fit(as.integer(p[2])), logical(1))
for (p in pairs[!keep_pair]) cat(sprintf("  [skip] '%s' -- missing fit in %s\n", p[3], fits_dir))
pairs <- pairs[keep_pair]
if (!length(pairs)) stop("no pair has both fits present in ", fits_dir)

# ---- per-site twCRPS: integral over the band of per-cell Brier, summed over days ----
site_twcrps <- function(sid) {
  le <- new.env(); load(file.path(fits_dir, sprintf("val-%d.RData", sid)), envir = le)
  yp <- get("fit", le)$yp                        # draws x np x nt
  ypm <- matrix(yp, nrow = dim(yp)[1])           # draws x (np*nt), cols site-fastest (match yv)
  acc <- numeric(length(yv)); prevBS <- NULL
  for (k in seq_along(zgrid)) {
    z <- zgrid[k]
    phat <- colMeans(ypm > z)                    # P(Y>z) per cell
    ind  <- as.numeric(yv > z)                   # NA where yv NA -> BS NA -> cell drops out
    BS   <- (phat - ind)^2
    if (k > 1) acc <- acc + 0.5 * (BS + prevBS) * (zgrid[k] - zgrid[k - 1])  # trapezoid
    prevBS <- BS
  }
  cat(sprintf("  scored setting %-4d draws=%d\n", sid, dim(yp)[1]))
  rowSums(matrix(acc, np, nt), na.rm = TRUE)     # per-site twCRPS sum (length np)
}

need <- unique(unlist(lapply(pairs, function(p) as.integer(p[1:2]))))
TW <- setNames(lapply(need, site_twcrps), as.character(need))

boot_pair <- function(sidA, sidB) {
  a <- TW[[as.character(sidA)]][keep]; b <- TW[[as.character(sidB)]][keep]
  d <- a - b; w <- nday[keep]
  ratio <- function(idx) sum(d[idx]) / sum(w[idx])
  point <- ratio(seq_along(d))
  set.seed(seed)
  bootv <- replicate(B, ratio(sample.int(length(d), replace = TRUE)))
  ci <- quantile(bootv, c(0.025, 0.975), names = FALSE)
  p  <- min(1, 2 * min(mean(bootv <= 0), mean(bootv >= 0)))
  data.frame(twCRPS_A = sum(a) / sum(w), twCRPS_B = sum(b) / sum(w),
             Delta = point, CI_lo = ci[1], CI_hi = ci[2], boot_p = p)
}

out <- do.call(rbind, lapply(pairs, function(p) {
  cbind(pair = p[3], A = as.integer(p[1]), B = as.integer(p[2]),
        boot_pair(as.integer(p[1]), as.integer(p[2])))
}))

fmt <- out
for (c0 in c("twCRPS_A", "twCRPS_B", "Delta", "CI_lo", "CI_hi")) fmt[[c0]] <- sprintf("%+.4f", out[[c0]])
fmt$boot_p <- sprintf("%.3f", out$boot_p)
fmt$sig <- ifelse(out$CI_lo > 0 | out$CI_hi < 0, "*", "")
cat("\n=== Paired tail-band twCRPS differences (Delta = A - B; negative => A better) ===\n")
cat("=== 95% site-clustered bootstrap CI; '*' = CI excludes 0 ===\n\n")
print(fmt[, c("pair", "twCRPS_A", "twCRPS_B", "Delta", "CI_lo", "CI_hi", "boot_p", "sig")], row.names = FALSE)

dir.create("output/us-all-auto/tables", recursive = TRUE, showWarnings = FALSE)
outfile <- sprintf("output/us-all-auto/tables/twcrps_bootstrap_pairs_%s.csv", out_tag)
write.csv(out, outfile, row.names = FALSE)
cat(sprintf("\nWritten: %s\n", outfile))
