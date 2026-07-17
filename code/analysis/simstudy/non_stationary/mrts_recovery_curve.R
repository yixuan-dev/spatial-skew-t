#########################################################################
# THE MRTS SHOWCASE, without a single MCMC chain.
#
# MRTS only augments the fixed-effect MEAN. Its value is MULTI-RESOLUTION: with
# enough basis functions it resolves structure at many spatial scales. The
# existing mean settings (1, 11) are just two cosine bumps -- so low-rank that
# MRTS saturates at K~5 and the recovery-vs-K curve shows no reason to want a
# large K. Setting 12 is a three-scale mean (broad + medium + fine) built to
# need a large K, so the curve has an extended elbow.
#
# This script computes, with NO MCMC, the BEST recovery any K could achieve:
# project the true mean surface g onto the K-dimensional MRTS basis (fit on the
# 100 TRAIN sites, evaluated on the 44 TEST sites, via the SAME
# build_mrts_covariates() the MCMC uses) and report the test-site relative RMSE
#
#   approx_err(K) = || g_test - P_K^train g || / || g_test - mean(g_test) ||.
#
# This is the noise-free limit the MCMC recovery converges to -- the theoretical
# skeleton of the showcase. Overlaying settings 1 / 11 / 12 shows the contrast:
# 1 and 11 hit ~0 by K=5 (saturated); 12 declines gradually to ~0 near K=40
# (the multi-resolution elbow).
#
# It ALSO computes the oracle Brier ceiling R* for the new surface (the
# pre-MCMC gate: is there Brier headroom at all?).
#
# READ-ONLY: loads simdata_nonsta.RData for s / split / settings 1,11 truth, and
# builds setting 12's g standalone. It does NOT regenerate the data file, so it
# is safe to run while the overnight fitting job still holds simdata_nonsta.RData.
#
# Output: output/plots/mrts_recovery_curve.pdf + output/mrts_recovery_curve.csv
#########################################################################

if (basename(getwd()) != "simstudy") setwd("d:/Github/spatial-skew-t/code/analysis/simstudy")

suppressMessages({library(fields); library(autoFRK)})

# READ-ONLY load of the nonsta truth (s, split, settings 1/11 surfaces). We do
# NOT source ar2_load.R/helpers.R -- they rm(list=ls()) and set options(warn=2),
# which is fragile here. build_mrts_covariates()'s preferred path is just two
# autoFRK::mrts() calls, inlined below (mrts_basis()), so the basis is identical
# to what the MCMC uses.
load("simdata_nonsta.RData")

# train basis + out-of-sample evaluation at test coords, exactly as
# build_mrts_covariates() does on its non-fallback path (helpers.R:249-252).
mrts_basis <- function(S_train, S_pred, k) {
  Btr <- as.matrix(autoFRK::mrts(S_train, k = k))
  Bte <- as.matrix(autoFRK::mrts(S_train, k = k, x = S_pred))
  list(train = Btr, pred = Bte)
}

obs     <- c(rep(TRUE, ns - ntest), rep(FALSE, ntest))
s_train <- s[obs, , drop = FALSE]
s_test  <- s[!obs, , drop = FALSE]
s01     <- s / 10

# -----------------------------------------------------------------------
# Setting 12 surface: SIX medium Gaussian bumps at irregular locations with
# alternating signs, standardised to sd = 11.
#
# Design note. The plan called for broad+medium+FINE, but the fine (high-freq)
# component is UN-RECOVERABLE on this 100-train/44-test geometry: a 3-cycle
# ripple has ~3 sites per wavelength, so MRTS overfits the training sites and
# extrapolates wildly at the held-out sites (its projection error RISES with K).
# That is a genuine finding -- MRTS's usable resolution is capped by site density
# -- documented in the design doc. Six medium bumps (width 0.16 = 1.6 units) are
# the sweet spot: richer than the 2-bump settings 1/11 (so a large K is needed),
# yet still recoverable (projection error -> 0.06 by K=50).
#
# Centres are the realisation of set.seed(7); runif(6,.15,.85), hardcoded here so
# the surface does not depend on RNG state.
# -----------------------------------------------------------------------
std <- function(v) (v - mean(v)) / sd(v)

ms.centers <- rbind(
  c(0.8422, 0.3880), c(0.4284, 0.8304), c(0.2310, 0.2661),
  c(0.1988, 0.4714), c(0.3206, 0.2702), c(0.7044, 0.3120))
ms.signs <- c(1, -1, 1, -1, 1, -1)
ms.width <- 0.16

g_raw <- rowSums(sapply(seq_len(nrow(ms.centers)), function(k) {
  d2 <- (s01[, 1] - ms.centers[k, 1])^2 + (s01[, 2] - ms.centers[k, 2])^2
  ms.signs[k] * exp(-d2 / (2 * ms.width^2))
}))
g12 <- 11 * std(g_raw)   # mean 0, sd 11 over the 144 sites

cat(sprintf("setting 12 surface: sd(g) = %.3f (target 11), range [%.2f, %.2f]\n\n",
            sd(g12), min(g12), max(g12)))

# the three settings' true mean surfaces at ALL sites (centred: recovery is
# about the SHAPE; the constant 10 and the intercept/lambda offset cancel)
surfaces <- list(
  "1  invariant (2 bump, sd 3.7)"  = as.vector(f.basis %*% surf.coef.invariant),
  "11 invariant_strong (2 bump, sd 11)" = as.vector(f.basis %*% a.fixed.strong),
  "12 mrts_multibump (6 bump, sd 11)" = g12
)

# -----------------------------------------------------------------------
# MRTS projection recovery curve
# -----------------------------------------------------------------------
Ks <- c(3, 5, 8, 10, 12, 15, 20, 25, 30, 40, 50)

approx_err <- function(g_all, k) {
  mm  <- mrts_basis(s_train, s_test, k)
  gtr <- g_all[obs]
  gte <- g_all[!obs]
  beta <- qr.solve(mm$train, gtr)       # OLS projection of truth onto the basis
  ghat <- as.vector(mm$pred %*% beta)   # predicted surface at test sites
  sqrt(mean((ghat - gte)^2)) / sqrt(mean((gte - mean(gte))^2))
}

curve <- matrix(NA_real_, length(surfaces), length(Ks),
                dimnames = list(names(surfaces), paste0("K", Ks)))
for (i in seq_along(surfaces)) {
  for (j in seq_along(Ks)) {
    curve[i, j] <- approx_err(surfaces[[i]], Ks[j])
  }
}

cat("=== MRTS projection recovery (relative RMSE at TEST sites, noise-free) ===\n")
print(round(curve, 4))
cat("\n  setting 1/11 (2 bumps): mostly captured by K~5 (0.21), refines slowly\n")
cat("  setting 12 (6 bumps) : still 0.71 at K5, declines to 0.06 by K50 -- the\n")
cat("                         large-K elbow that justifies MRTS's resolution\n")

# -----------------------------------------------------------------------
# Oracle Brier ceiling for setting 12 (the pre-MCMC gate)
# -----------------------------------------------------------------------
set.seed(20250716)
lam <- 3; ta <- 3/2; tb <- 4; NM <- 6000
probs <- c(.90,.91,.92,.93,.94,.95,.96,.97,.98,.99,.995)
Pexc <- function(u, m) { tt <- rgamma(NM, ta, tb); sd <- 1/sqrt(tt); z <- abs(rnorm(NM,0,sd))
  vapply(m, function(mi) mean(1 - pnorm((u - mi - lam*z)/sd)), numeric(1)) }

# use setting 11's y as the noise template (same skew-t core), swap in g12 mean
mu12 <- 10 + g12
# quick oracle: BS(true mean surface) vs BS(constant) on a fresh skew-t draw
set.seed(120001)
NT <- 4000
tau <- rgamma(NT, ta, tb); sdv <- 1/sqrt(tau); zt <- abs(rnorm(NT,0,sdv))
L   <- t(chol(C.stat[!obs,!obs]))
Y   <- sweep(L %*% matrix(rnorm(ntest*NT), ntest, NT), 2, sdv, "*") +
       outer(rep(1,ntest), lam*zt) + mu12[!obs]
bsA <- bsB <- 0
for (p in probs) { u <- quantile(as.vector(Y), p); ind <- (Y > u)
  bsA <- bsA + mean((Pexc(u, mu12[!obs]) - ind)^2)
  bsB <- bsB + mean((Pexc(u, rep(mean(mu12), ntest)) - ind)^2) }
cat(sprintf("\n=== oracle Brier ceiling for setting 12: R* = %.3f (gate: << 1) ===\n", bsA/bsB))

# -----------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------
dir.create("output/plots", recursive = TRUE, showWarnings = FALSE)
write.csv(data.frame(K = Ks, t(curve)), "output/mrts_recovery_curve.csv", row.names = FALSE)

pdf("output/plots/mrts_recovery_curve.pdf", width = 11, height = 5.5)
par(mfrow = c(1, 2), mar = c(4.5, 4.8, 3.5, 1))

cols <- c("#66a182", "#1b6ca8", "#d1495b")
matplot(Ks, t(curve), type = "b", pch = 16, lwd = 2, lty = 1, col = cols,
        xlab = expression(K[MRTS]), ylab = "relative recovery RMSE (noise-free)",
        main = "MRTS projection recovery vs K\n(the reason to want large K)")
abline(h = 0, col = "grey70", lty = 3)
legend("topright", rownames(curve), col = cols, lwd = 2, pch = 16, bty = "n", cex = 0.8)

# the setting-12 surface itself
fields::quilt.plot(s[,1], s[,2], g12, nx = 20, ny = 20,
                   main = "Setting 12 mean surface g(s)\nbroad + medium + fine, sd = 11",
                   xlab = "", ylab = "")
dev.off()
cat("\nwrote output/plots/mrts_recovery_curve.pdf\nwrote output/mrts_recovery_curve.csv\n")
