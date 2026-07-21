# mc_error_check.R
# ---------------------------------------------------------------------------
# Measure (not assert) the posterior Monte-Carlo error in the pooled scores and
# their paired difference, to substantiate folder 01's claim that MC error is
# negligible next to the evaluation (bootstrap) uncertainty.
#
# Batch-means estimator: split the M posterior draws into G contiguous batches,
# recompute the pooled score from each batch's P-hat, take the across-batch SD;
# the MC s.e. of the full-sample pooled score is SD_batch / sqrt(G). The two
# models are independent MCMC runs, so MC_SE(Delta) = sqrt(SE_A^2 + SE_B^2).
# Compare to the site-bootstrap CI half-width for the same pair/score.
#
# Ozone AR(2) pair (111 vs 61), 200/200 split, scores: Brier q0.95, CRPS, twCRPS.
# Writes: output/us-all-auto/tables/mc_error_check.csv
# ---------------------------------------------------------------------------

rm(list = ls())
.this <- local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), winslash = "/", mustWork = FALSE)) else "." })
setwd(.this)
.r_root <- normalizePath(file.path(.this, "../../../R"), winslash = "/", mustWork = FALSE)
source(local({ ar2 <- file.path(.r_root, "ar2", "auxfunctions.R"); if (file.exists(ar2)) ar2 else file.path(.r_root, "auxfunctions.R") }))
stopifnot(requireNamespace("scoringRules", quietly = TRUE))

fits_dir <- "fits-200-200"; G <- 10L
probs <- c(0.90, 0.95, 0.98, 0.99, 0.995); q95 <- 2L
band <- NULL; zG <- 41L
pairA <- 111L; pairB <- 61L

se <- new.env(); load("us-all-setup-auto-200-200-400.RData", envir = se)
Y <- get("Y", se); val <- get("split.lst", se)$val
Y_val <- Y[val, , drop = FALSE]; np <- nrow(Y_val); nt <- ncol(Y_val)
nday <- rowSums(!is.na(Y_val)); tot_days <- sum(nday)
thr <- quantile(Y, probs = probs, na.rm = TRUE)
bandv <- quantile(Y, probs = c(0.90, 0.995), na.rm = TRUE); zgrid <- seq(bandv[1], bandv[2], length.out = zG)
yv <- as.vector(Y_val); ok <- !is.na(yv)

# pooled scores from a given set of draw indices (a batch)
pooled_scores <- function(ypm, idx) {
  sub <- ypm[idx, , drop = FALSE]                       # (batch draws) x cells
  # Brier q0.95
  phat <- colMeans(sub > thr[q95]); ind <- as.numeric(yv > thr[q95])
  brier <- mean((phat - ind)^2, na.rm = TRUE)
  # plain CRPS
  cc <- rep(NA_real_, length(yv)); cc[ok] <- scoringRules::crps_sample(y = yv[ok], dat = t(sub[, ok, drop = FALSE]))
  crps <- mean(cc, na.rm = TRUE)
  # tail-band twCRPS
  acc <- numeric(length(yv)); prev <- NULL
  for (k in seq_along(zgrid)) { BS <- (colMeans(sub > zgrid[k]) - as.numeric(yv > zgrid[k]))^2
    if (k > 1) acc <- acc + 0.5 * (BS + prev) * (zgrid[k] - zgrid[k - 1]); prev <- BS }
  twcrps <- mean(acc, na.rm = TRUE)
  c(brier = brier, crps = crps, twcrps = twcrps)
}

batch_se <- function(sid) {
  yp <- get("fit", local({ le <- new.env(); load(file.path(fits_dir, sprintf("val-%d.RData", sid)), envir = le); le }))$yp
  ypm <- matrix(yp, nrow = dim(yp)[1])
  M <- nrow(ypm); bat <- split(seq_len(M), cut(seq_len(M), G, labels = FALSE))
  sc <- t(vapply(bat, function(ix) pooled_scores(ypm, ix), numeric(3)))  # G x 3
  full <- pooled_scores(ypm, seq_len(M))
  cat(sprintf("  %-4d done\n", sid))
  list(full = full, mc_se = apply(sc, 2, sd) / sqrt(G))
}

cat("batch-means MC error (G =", G, "batches)\n")
A <- batch_se(pairA); B <- batch_se(pairB)

# bootstrap CI half-widths for the same pair (from existing CSVs)
halfwidth <- function(f, score_row) {
  d <- read.csv(file.path("output/us-all-auto/tables", f), stringsAsFactors = FALSE)
  r <- score_row(d); (r$CI_hi - r$CI_lo) / 2
}
hw <- c(
  brier = halfwidth("brier_bootstrap_pairs_200-200.csv",
                    function(d) d[startsWith(d$pair, "AR(2)") & abs(d$level - 0.95) < 1e-6, ]),
  crps  = halfwidth("crps_bootstrap_pairs_200-200.csv",  function(d) d[startsWith(d$pair, "AR(2)"), ]),
  twcrps= halfwidth("twcrps_bootstrap_pairs_200-200.csv",function(d) d[startsWith(d$pair, "AR(2)"), ]))

sc_names <- c("brier", "crps", "twcrps")
Delta <- A$full - B$full
mc_se_Delta <- sqrt(A$mc_se^2 + B$mc_se^2)
out <- data.frame(
  score = sc_names, Delta = Delta[sc_names],
  mc_se_Delta = mc_se_Delta[sc_names], boot_CI_halfwidth = hw[sc_names],
  ratio_mc_to_boot = mc_se_Delta[sc_names] / hw[sc_names], row.names = NULL)

cat("\n=== MC error of Delta vs bootstrap CI half-width (ozone AR(2), 200/200) ===\n\n")
print(within(out, {
  Delta <- signif(Delta, 3); mc_se_Delta <- signif(mc_se_Delta, 3)
  boot_CI_halfwidth <- signif(boot_CI_halfwidth, 3); ratio_mc_to_boot <- signif(ratio_mc_to_boot, 2)
}), row.names = FALSE)

dir.create("output/us-all-auto/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(out, "output/us-all-auto/tables/mc_error_check.csv", row.names = FALSE)
cat("\nWritten: output/us-all-auto/tables/mc_error_check.csv\n")
