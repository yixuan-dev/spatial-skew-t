# block_bootstrap_sensitivity.R
# ---------------------------------------------------------------------------
# Sensitivity of the paired-score CIs to BETWEEN-SITE spatial correlation.
#
# The site-clustered bootstrap in brier/crps/twcrps_bootstrap.R resamples
# individual validation SITES, assuming they are approximately independent
# clusters. Nearby ozone monitors are spatially correlated (the skew-t process
# is built to make them jointly extreme), so the effective number of clusters
# is < n_site and the site-bootstrap CI is anticonservative (too narrow).
#
# This script re-runs the paired difference for the ozone AR(2) and MRTS pairs
# under a spatial BLOCK bootstrap: the domain is gridded into contiguous blocks,
# and whole blocks (all their sites and days) are resampled together, so
# spatially-correlated neighbours stay together. Coarser blocks absorb longer-
# range correlation and widen the CI. We report site vs block (coarse/med/fine)
# side by side for three scores (Brier q0.95, tail-band twCRPS, plain CRPS).
#
# Reading: for a NULL, an anticonservative site CI is safe (a wider block CI
# still covers 0). The one at-risk cell is the significant plain-CRPS
# "AR(2) worse"; the test is whether it survives coarser blocks.
#
# Writes: output/us-all-auto/tables/block_bootstrap_sensitivity_<split>.csv
# ---------------------------------------------------------------------------

rm(list = ls())
.this <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), winslash = "/", mustWork = FALSE)) else "."
})
setwd(.this)
.r_root <- normalizePath(file.path(.this, "../../../R"), winslash = "/", mustWork = FALSE)
source(local({ ar2 <- file.path(.r_root, "ar2", "auxfunctions.R"); if (file.exists(ar2)) ar2 else file.path(.r_root, "auxfunctions.R") }))
stopifnot(requireNamespace("scoringRules", quietly = TRUE))

setup_file <- Sys.getenv("US_ALL_SETUP_FILE",      unset = "us-all-setup-auto-200-200-400.RData")
fits_dir   <- Sys.getenv("US_ALL_VAL_RESULTS_DIR", unset = "fits-200-200")
out_tag    <- Sys.getenv("US_ALL_BOOT_TAG",        unset = "200-200")
B <- 2000L; seed <- 2024L
probs <- c(0.90, 0.95, 0.98, 0.99, 0.995); q95_idx <- 2L
band_probs <- c(0.90, 0.995); G <- 41L
pairs <- list(c(111, 61, "AR(2) vs AR(1)"), c(204, 51, "MRTS vs none"))
grids <- list(site = NULL, block_coarse = c(4, 3), block_med = c(6, 4), block_fine = c(8, 6))

# ---- data ----
se <- new.env(); load(setup_file, envir = se)
Y <- get("Y", se); S <- get("S", se); val <- get("split.lst", se)$val
Y_val <- Y[val, , drop = FALSE]; np <- nrow(Y_val); nt <- ncol(Y_val)
nday <- rowSums(!is.na(Y_val)); keep <- nday > 0
thresholds <- quantile(Y, probs = probs, na.rm = TRUE)
band <- quantile(Y, probs = band_probs, na.rm = TRUE); zgrid <- seq(band[1], band[2], length.out = G)
coords <- S[val, , drop = FALSE][keep, , drop = FALSE]
w <- nday[keep]

has_fit <- function(sid) file.exists(file.path(fits_dir, sprintf("val-%d.RData", sid)))
pairs <- pairs[vapply(pairs, function(p) has_fit(as.integer(p[1])) && has_fit(as.integer(p[2])), logical(1))]

# ---- per-site score sums (SSE-analogue: sum over a site's days), all three scores ----
per_site_all <- function(sid) {
  yp <- get("fit", local({ le <- new.env(); load(file.path(fits_dir, sprintf("val-%d.RData", sid)), envir = le); le }))$yp
  d <- dim(yp)
  brier_sse <- BrierScoreSite(yp, thresholds, Y_val)[, q95_idx] * nday      # Brier @ q0.95
  ypm <- matrix(yp, nrow = d[1]); yv <- as.vector(Y_val); ok <- !is.na(yv)
  # plain CRPS per cell
  crps_cell <- rep(NA_real_, length(yv))
  crps_cell[ok] <- scoringRules::crps_sample(y = yv[ok], dat = t(ypm[, ok, drop = FALSE]))
  crps_sse <- rowSums(matrix(crps_cell, np, nt), na.rm = TRUE)
  # tail-band twCRPS per cell (trapezoid integral of Brier over zgrid)
  acc <- numeric(length(yv)); prevBS <- NULL
  for (k in seq_along(zgrid)) {
    BS <- (colMeans(ypm > zgrid[k]) - as.numeric(yv > zgrid[k]))^2
    if (k > 1) acc <- acc + 0.5 * (BS + prevBS) * (zgrid[k] - zgrid[k - 1]); prevBS <- BS
  }
  twcrps_sse <- rowSums(matrix(acc, np, nt), na.rm = TRUE)
  cat(sprintf("  scored %-4d\n", sid))
  cbind(brier = brier_sse, crps = crps_sse, twcrps = twcrps_sse)[keep, ]
}

need <- unique(unlist(lapply(pairs, function(p) as.integer(p[1:2]))))
SSE <- setNames(lapply(need, per_site_all), as.character(need))

# ---- block assignment ----
assign_blocks <- function(co, nx, ny) {
  xr <- range(co[, 1]); yr <- range(co[, 2])
  xb <- pmin(nx, 1L + floor(nx * (co[, 1] - xr[1]) / (diff(xr) + 1e-9)))
  yb <- pmin(ny, 1L + floor(ny * (co[, 2] - yr[1]) / (diff(yr) + 1e-9)))
  as.integer((yb - 1L) * nx + xb)
}

boot_ratio <- function(d, w, grp) {
  units <- unique(grp); idx <- split(seq_along(d), grp)
  set.seed(seed)
  bootv <- replicate(B, { drawn <- sample(units, replace = TRUE)
    s <- unlist(idx[as.character(drawn)], use.names = FALSE); sum(d[s]) / sum(w[s]) })
  ci <- quantile(bootv, c(.025, .975), names = FALSE)
  list(Delta = sum(d) / sum(w), lo = ci[1], hi = ci[2],
       p = min(1, 2 * min(mean(bootv <= 0), mean(bootv >= 0))), n_units = length(units))
}

rows <- list()
for (p in pairs) {
  A <- as.character(as.integer(p[1])); Bx <- as.character(as.integer(p[2]))
  for (sc in c("brier", "crps", "twcrps")) {
    d <- SSE[[A]][, sc] - SSE[[Bx]][, sc]
    for (gname in names(grids)) {
      grp <- if (is.null(grids[[gname]])) seq_along(d) else assign_blocks(coords, grids[[gname]][1], grids[[gname]][2])
      r <- boot_ratio(d, w, grp)
      rows[[length(rows) + 1]] <- data.frame(
        pair = p[3], score = sc, scheme = gname, n_units = r$n_units,
        Delta = r$Delta, CI_lo = r$lo, CI_hi = r$hi, boot_p = r$p,
        sig = if (r$lo > 0 || r$hi < 0) "*" else "", stringsAsFactors = FALSE)
    }
  }
}
out <- do.call(rbind, rows)

fmt <- out
for (c0 in c("Delta", "CI_lo", "CI_hi")) fmt[[c0]] <- formatC(out[[c0]], format = "fg", digits = 2, flag = "+")
fmt$boot_p <- sprintf("%.3f", out$boot_p)
cat("\n=== site vs spatial-block bootstrap (ozone ", out_tag, ") ===\n", sep = "")
cat("=== Delta = A - B (neg => A better); '*' CI excludes 0 ===\n\n")
print(fmt[, c("pair", "score", "scheme", "n_units", "Delta", "CI_lo", "CI_hi", "boot_p", "sig")], row.names = FALSE)

dir.create("output/us-all-auto/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(out, sprintf("output/us-all-auto/tables/block_bootstrap_sensitivity_%s.csv", out_tag), row.names = FALSE)
cat(sprintf("\nWritten: output/us-all-auto/tables/block_bootstrap_sensitivity_%s.csv\n", out_tag))
