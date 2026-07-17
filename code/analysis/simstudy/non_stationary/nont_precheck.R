#########################################################################
# Pre-MCMC gates for the non-t settings 13-21, run BEFORE setup_nonsta.R
# is extended. Two cheap checks, no MCMC:
#
# (a) noise-free MRTS projection recovery vs K for each of the nine new
#     mean surfaces (same construction as mrts_recovery_curve.R): confirms
#     every surface is recoverable on the 100-train/44-test geometry, that
#     setting 13's fine dip and 16's front do not blow up in extrapolation,
#     and that 21 plateaus (expected: C^0 creases vs C^inf basis).
#
# (b) oracle Brier ceiling R*(sigma_G) on a grid sigma_G in {3, 5, 7.5}
#     for the Gaussian-error settings: pick the shared sigma_G that puts
#     R* in the 0.4-0.5 band (comparable to settings 11/12), to be
#     hardcoded in setup_nonsta.R.
#
# READ-ONLY: loads simdata_nonsta.RData only for the site geometry and
# C.stat; does not write any data file. Safe to run while fitting jobs
# hold simdata_nonsta.RData.
#
# Output: output/nont_precheck.csv (+ _oracle.csv) and
#         output/plots/nont_precheck.pdf
#########################################################################

if (basename(getwd()) != "simstudy") setwd("d:/Github/spatial-skew-t/code/analysis/simstudy")

suppressMessages({library(fields); library(autoFRK)})

load("simdata_nonsta.RData")          # s, ns, ntest, C.stat, ...
source("non_stationary/nont_builders.R")

obs     <- c(rep(TRUE, ns - ntest), rep(FALSE, ntest))
s_train <- s[obs, , drop = FALSE]
s_test  <- s[!obs, , drop = FALSE]
s01     <- s / 10

# -----------------------------------------------------------------------
# The nine surfaces (settings 13-21; 18 shares 13's surface). All pushed
# through ortho_std() exactly as setup_nonsta.R will.
# -----------------------------------------------------------------------
surf <- list(
  "13 franke"        = ortho_std(franke_fn(s01),        s),
  "14 ridge_curved"  = ortho_std(ridge_fn(s01),         s),
  "15 annulus"       = ortho_std(annulus_fn(s, width = 1.5), s),
  "16 sigmoid_front" = ortho_std(front_fn(s),           s),
  "19 poly2"         = ortho_std(poly2_fn(s01),         s),
  "20 transform_mix" = ortho_std(transform_mix_fn(s),   s),
  "21 nonsmooth"     = ortho_std(nonsmooth_fn(s),       s)
)

# setting 17: three representative GP draws (the setup will use seeds
# 860000 + set; use the first three here)
C.gpmean <- matern_cor(s, range = 2, smoothness = 1.5)
L.gpmean <- t(chol(C.gpmean))
for (k in 1:3) {
  set.seed(860000 + k)
  surf[[sprintf("17 gp_mean (draw %d)", k)]] <- ortho_std(as.vector(L.gpmean %*% rnorm(ns)), s)
}

# -----------------------------------------------------------------------
# (a) MRTS projection recovery curves
# -----------------------------------------------------------------------
mrts_basis <- function(S_train, S_pred, k) {
  Btr <- as.matrix(autoFRK::mrts(S_train, k = k))
  Bte <- as.matrix(autoFRK::mrts(S_train, k = k, x = S_pred))
  list(train = Btr, pred = Bte)
}

Ks <- c(3, 5, 8, 10, 12, 15, 20, 25, 30, 40, 50)

approx_err <- function(g_all, k) {
  mm  <- mrts_basis(s_train, s_test, k)
  beta <- qr.solve(mm$train, g_all[obs])
  ghat <- as.vector(mm$pred %*% beta)
  gte  <- g_all[!obs]
  sqrt(mean((ghat - gte)^2)) / sqrt(mean((gte - mean(gte))^2))
}

curve <- matrix(NA_real_, length(surf), length(Ks),
                dimnames = list(names(surf), paste0("K", Ks)))
for (i in seq_along(surf)) {
  for (j in seq_along(Ks)) {
    curve[i, j] <- approx_err(surf[[i]], Ks[j])
  }
}

cat("=== (a) MRTS projection recovery (relative RMSE at TEST sites, noise-free) ===\n")
print(round(curve, 4))

# -----------------------------------------------------------------------
# (b) oracle Brier ceiling R* vs sigma_G, Gaussian error with correlation
# C.stat. Predictive law at test site i is N(10 + g_i, sigma^2); thresholds
# are pooled test-sample quantiles, as in scores.R.
# -----------------------------------------------------------------------
probs   <- c(.90, .91, .92, .93, .94, .95, .96, .97, .98, .99, .995)
sigmas  <- c(3, 4, 5, 7.5)
L_test  <- t(chol(C.stat[!obs, !obs]))
NT      <- 4000

oracle_Rstar <- function(g_all, sigma) {
  mu <- 10 + g_all[!obs]
  set.seed(130001)
  Y  <- mu + sigma * (L_test %*% matrix(rnorm(ntest * NT), ntest, NT))
  bsA <- bsB <- 0
  for (p in probs) {
    u   <- quantile(as.vector(Y), p)
    ind <- (Y > u)
    pA  <- 1 - pnorm((u - mu) / sigma)             # true mean surface
    pB  <- 1 - pnorm((u - mean(mu)) / sigma)       # spatial constant
    bsA <- bsA + mean((pA - ind)^2)
    bsB <- bsB + mean((pB - ind)^2)
  }
  bsA / bsB
}

oracle <- sapply(sigmas, function(sg) sapply(surf, oracle_Rstar, sigma = sg))
colnames(oracle) <- paste0("sigma", sigmas)

cat("\n=== (b) oracle Brier ceiling R* (target band 0.4-0.5) ===\n")
print(round(oracle, 3))
cat(sprintf("\ncolumn medians: %s\n",
            paste(sprintf("sigma=%.1f -> %.3f", sigmas, apply(oracle, 2, median)),
                  collapse = ", ")))

# -----------------------------------------------------------------------
# Output
# -----------------------------------------------------------------------
dir.create("output/plots", recursive = TRUE, showWarnings = FALSE)
write.csv(data.frame(surface = rownames(curve), round(curve, 4)),
          "output/nont_precheck.csv", row.names = FALSE)
write.csv(data.frame(surface = rownames(oracle), round(oracle, 4)),
          "output/nont_precheck_oracle.csv", row.names = FALSE)

pdf("output/plots/nont_precheck.pdf", width = 11, height = 8)
par(mfrow = c(1, 1), mar = c(4.5, 4.8, 3.5, 1))
cols <- rep(c("#1b6ca8", "#d1495b", "#66a182", "#edae49", "#775b9f",
              "#2e4057", "#00798c", "#9a4c50", "#c17fa0", "#5c821a"),
            length.out = nrow(curve))
matplot(Ks, t(curve), type = "b", pch = 16, lwd = 2, lty = 1, col = cols,
        xlab = expression(K[MRTS]), ylab = "relative recovery RMSE (noise-free)",
        main = "Settings 13-21: MRTS projection recovery vs K")
abline(h = 0, col = "grey70", lty = 3)
legend("topright", rownames(curve), col = cols, lwd = 2, pch = 16, bty = "n", cex = 0.75)

par(mfrow = c(2, 3), mar = c(2, 2, 3, 1))
for (nm in names(surf)) {
  fields::quilt.plot(s[, 1], s[, 2], surf[[nm]], nx = 20, ny = 20, main = nm,
                     xlab = "", ylab = "")
}
dev.off()

cat("\nwrote output/nont_precheck.csv, output/nont_precheck_oracle.csv,\n")
cat("      output/plots/nont_precheck.pdf\n")
