#########################################################################
# Diagnostics for simdata_nonsta.RData -- run BEFORE spending any MCMC.
#
# Each setting claims to be non-stationary in a particular moment.  This script
# checks the claim is true in the generated data, that the strength is
# detectable but not absurd, and that the settings which only mean to change the
# DEPENDENCE really did leave the marginal alone.
#
# TWO RESIDUALS, because the diagnostics want opposite things.
#
# The model is  y_t(s) = 10 + mean_t(s) + lambda(s) z_t + sigma_t eps_t(s) + u_t(s),
# with z_t a per-day scalar (nknots = 1) and u_t whatever additive structure the
# setting injects.
#
#  R.dep   = y - 10 - mean_t(s) - lambda(s) z_t
#            Removes the intercept, the mean surface and the skew term EXACTLY,
#            using the saved truth (z.t, W.varying, lambda.field).  What is left
#            is sigma_t eps_t + u_t: the dependence signal, and nothing else.
#            Removing lambda(s) z_t matters -- it is common to every site on day
#            t, so it puts a large distance-free positive correlation between
#            every pair and would swamp the variogram.
#
#  R.shape = y - truemean.field
#            Removes only the TIME-AVERAGED true mean, so the skew term SURVIVES.
#            This is the residual whose per-site skewness/kurtosis carries the
#            higher-moment structure of settings 9 and 10.
#            (Subtracting the true mean, rather than the per-day spatial mean, is
#            what keeps setting 1's mean surface from interacting with the
#            per-day sigma_t normalisation and faking a skewness pattern.)
#
# Both are then centred per day and, where a shape/correlation statistic is
# wanted, divided by their per-day spatial sd -- which removes sigma_t.  sigma_t
# is heavy-tailed (tau ~ Gamma(3/2, 4) => df 3) and would otherwise let a handful
# of days decide every average.
#
# Outputs: output/nonsta_diagnostics.pdf and output/nonsta_diagnostics.csv
#########################################################################

if (basename(getwd()) != "simstudy") setwd("d:/Github/spatial-skew-t/code/analysis/simstudy")

library(fields)
library(SpatialTools)

source("../../R/ar2/auxfunctions.R")
source("non_stationary/ns_cov_builders.R")

mem <- function(s, knots) {
  d <- fields::rdist(s, knots)
  apply(d, 1, which.min)
}

load("simdata_nonsta.RData")

dir.create("output", showWarnings = FALSE)
pdf_fp <- "output/nonsta_diagnostics.pdf"
csv_fp <- "output/nonsta_diagnostics.csv"

nsettings <- nrow(settings.nonsta)
stype     <- settings.nonsta$surf_type
d.mat     <- as.matrix(dist(s))
NREF      <- nsettings + 1L                # index of the stationary reference

# -----------------------------------------------------------------------
# A stationary reference realisation: the same skew-t core with NO structure at
# all.  Anything the diagnostics call "non-stationary" must look different from
# THIS; anything this reference also shows is an artefact of the skew-t core or
# of the finite sample, not of the data-generating process.
# -----------------------------------------------------------------------
cat("generating the stationary reference...\n")
y.ref <- array(NA, dim = c(ns, nt, nsets))
z.ref <- matrix(NA, nt, nsets)
for (set in seq_len(nsets)) {
  set.seed(990000 + set)
  ref <- rpotspatTS(
    nt = nt, x = x, s = s, beta = c(10, 0, 0), gamma = 0.9, nu = 0.5, rho = 1,
    phi.z = 0, phi.w = 0, phi.tau = 0, lambda = 3,
    tau.alpha = 3, tau.beta = 8, nknots = 1, dist = "t", cov.type = "matern"
  )
  y.ref[, , set] <- ref$y
  z.ref[, set]   <- ref$z[1, ]
}

# -----------------------------------------------------------------------
# Residuals
# -----------------------------------------------------------------------
centre_days <- function(R) sweep(R, 2, colMeans(R), "-")
std_days    <- function(R) sweep(R, 2, apply(R, 2, stats::sd), "/")

# R.dep: intercept, mean surface and skew term all removed using the truth.
dep_resid <- function(k) {
  R <- matrix(NA_real_, ns, nt * nsets)
  for (set in seq_len(nsets)) {
    ya <- if (k == NREF) y.ref[, , set] else y[, , set, k]
    zt <- if (k == NREF) z.ref[, set]   else z.t[[k]][1, , set]

    # lambda(s): a field only in setting 10, a constant 3 everywhere else
    lam <- if (k <= nsettings && stype[k] == "lambda_varying") {
      lambda.field
    } else {
      rep(lambda.par$lambda0, ns)
    }

    # the fixed/varying mean surface (settings 1 and 2 only)
    M <- if (k <= nsettings && stype[k] == "invariant") {
      matrix(as.vector(f.basis %*% surf.coef.invariant), ns, nt)
    } else if (k <= nsettings && stype[k] == "varying") {
      f.basis %*% t(W.varying[[set]])
    } else {
      0
    }

    R[, (set - 1) * nt + seq_len(nt)] <- ya - 10 - M - outer(lam, zt)
  }
  centre_days(R)
}

# R.shape: only the time-averaged true mean removed, so the skew term survives.
shape_resid <- function(k) {
  R <- matrix(NA_real_, ns, nt * nsets)
  for (set in seq_len(nsets)) {
    ya <- if (k == NREF) y.ref[, , set] else y[, , set, k]
    tm <- if (k == NREF) rep(10, ns)    else truemean.field[[k]][, set]
    R[, (set - 1) * nt + seq_len(nt)] <- ya - tm
  }
  std_days(centre_days(R))
}

skew_fn <- function(v) mean((v - mean(v))^3) / stats::sd(v)^3
kurt_fn <- function(v) mean((v - mean(v))^4) / stats::sd(v)^4 - 3

cat("building residuals...\n")
Rdep <- lapply(seq_len(NREF), dep_resid)
Rshp <- lapply(seq_len(NREF), shape_resid)
lab  <- c(sprintf("%d %s", seq_len(nsettings), stype), "ref stationary")

pw.sd   <- sapply(Rdep, function(R) apply(R, 1, stats::sd))       # dependence resid
pw.skew <- sapply(Rshp, function(R) apply(R, 1, skew_fn))         # shape resid
pw.kurt <- sapply(Rshp, function(R) apply(R, 1, kurt_fn))

# -----------------------------------------------------------------------
# Variogram / correlogram of the day-standardised dependence residual.
# -----------------------------------------------------------------------
vgm_of <- function(Rs) {
  V <- tcrossprod(Rs) / ncol(Rs)
  0.5 * (outer(diag(V), diag(V), "+") - 2 * V)
}

quad <- 1 + (s[, 1] >= 5) + 2 * (s[, 2] >= 5)          # 1=SW 2=SE 3=NW 4=NE
brk  <- seq(0, 7, by = 0.5)
mid  <- head(brk, -1) + diff(brk) / 2
ut   <- upper.tri(d.mat)

bin_mean <- function(M, keep) {
  b <- cut(d.mat[ut & keep], brk, include.lowest = TRUE)
  as.numeric(tapply(M[ut & keep], b, mean))
}

cat("computing variograms and correlograms...\n")
Rdep.s <- lapply(Rdep, std_days)
qv <- lapply(Rdep.s, function(Rs) {
  G <- vgm_of(Rs)
  sapply(1:4, function(q) bin_mean(G, outer(quad == q, quad == q, "&")))
})
pc <- sapply(Rdep.s, function(Rs) bin_mean(stats::cor(t(Rs)), matrix(TRUE, ns, ns)))

# -----------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------
pdf(pdf_fp, width = 13, height = 9)
qcol <- c("#1b6ca8", "#d1495b", "#66a182", "#e8a33d")

map_panel <- function(v, ttl, zlim = NULL) {
  fields::quilt.plot(s[, 1], s[, 2], v, nx = 16, ny = 16, zlim = zlim,
                     main = ttl, xlab = "", ylab = "", cex.main = 0.9)
}

# 0. the truth fields, to compare against the empirical panels below
par(mfrow = c(2, 3), mar = c(2, 2, 3, 4), oma = c(0, 0, 2, 0))
map_panel(f.basis[, 1], "f1(s)   <- mean in setting 1")
map_panel(f.basis[, 2], "f2(s)")
map_panel(rho.field,    "rho(s) = range   <- setting 6, SAME f1")
map_panel(g.field,      "g(s) = skewness   <- setting 9")
map_panel(h.field,      "h(s) = tail weight   <- setting 9")
map_panel(lambda.field, "lambda(s)   <- setting 10")
mtext("Ground truth. Settings 1 and 6 are driven by the SAME f1 -- once in the mean, once in the range.",
      outer = TRUE, cex = 0.95)

# 1. pointwise sd of the dependence residual: replace_C (4/6/8) and 9 stay FLAT
par(mfrow = c(3, 4), mar = c(2, 2, 3, 4), oma = c(0, 0, 2, 0))
zl <- quantile(pw.sd, c(0.01, 0.99))
for (k in seq_len(NREF)) map_panel(pw.sd[, k], lab[k], zlim = zl)
mtext("Pointwise sd of the dependence residual -- flat = homoscedastic (4/6/8/9 by design; 3/5/7 vary)",
      outer = TRUE, cex = 0.95)

# 2. quadrant variograms -- THE KEY FIGURE
par(mfrow = c(3, 4), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
for (k in seq_len(NREF)) {
  matplot(mid, qv[[k]], type = "l", lty = 1, lwd = 2, col = qcol,
          xlab = "distance", ylab = "semivariance", ylim = c(0, 1.7),
          main = lab[k], cex.main = 0.95)
  if (k == 1) legend("bottomright", c("SW", "SE", "NW", "NE"), lty = 1, lwd = 2,
                     col = qcol, bty = "n", cex = 0.8)
}
mtext("Variogram by quadrant -- curves SEPARATE <=> non-stationary dependence (compare the flat reference)",
      outer = TRUE, cex = 0.95)

# 3. correlogram -- setting 3 must go NEGATIVE at long range (a Matern cannot)
par(mfrow = c(1, 1), mar = c(5, 5, 5, 2), oma = c(0, 0, 0, 0))
cols <- c(rainbow(nsettings, end = 0.8), "black")
matplot(mid, pc, type = "l", lty = c(rep(1, nsettings), 2), lwd = 2, col = cols,
        xlab = "distance", ylab = "correlation",
        main = "Pooled correlation vs distance\nsetting 3's low-rank cosine basis goes NEGATIVE at long range -- a monotone Matern cannot")
abline(h = 0, col = "grey50", lty = 3)
legend("topright", lab, lty = c(rep(1, nsettings), 2), lwd = 2, col = cols,
       bty = "n", cex = 0.8)

# 4. pointwise skewness -- 9 tracks g(s) (sign change), 10 tracks lambda(s)
par(mfrow = c(3, 4), mar = c(2, 2, 3, 4), oma = c(0, 0, 2, 0))
zl <- quantile(pw.skew, c(0.01, 0.99))
for (k in seq_len(NREF)) map_panel(pw.skew[, k], lab[k], zlim = zl)
mtext("Pointwise skewness of the shape residual -- 9 should track g(s) (SIGN CHANGE), 10 should track lambda(s)",
      outer = TRUE, cex = 0.95)

# 5. pointwise excess kurtosis -- setting 9 should track h(s)
par(mfrow = c(3, 4), mar = c(2, 2, 3, 4), oma = c(0, 0, 2, 0))
zl <- quantile(pw.kurt, c(0.01, 0.99))
for (k in seq_len(NREF)) map_panel(pw.kurt[, k], lab[k], zlim = zl)
mtext("Pointwise excess kurtosis -- setting 9 should track h(s)", outer = TRUE, cex = 0.95)

# 6. marginal QQ -- the replace_C settings must MATCH the stationary reference
par(mfrow = c(1, 3), mar = c(5, 5, 4, 2), oma = c(0, 0, 2, 0))
qs   <- seq(0.001, 0.999, length.out = 400)
qref <- quantile(as.vector(y.ref), qs)
for (k in c(4, 6, 8)) {
  plot(qref, quantile(as.vector(y[, , , k]), qs), type = "l", lwd = 2, col = "#1b6ca8",
       xlab = "stationary reference quantile", ylab = sprintf("setting %d quantile", k),
       main = lab[k])
  abline(0, 1, col = "#d1495b", lty = 2)
}
mtext("replace_C settings changed the DEPENDENCE only -- the marginal must sit on the diagonal",
      outer = TRUE, cex = 0.95)

dev.off()

# -----------------------------------------------------------------------
# Summary table
# -----------------------------------------------------------------------
lag_i <- which.min(abs(mid - 1.5))     # a representative short lag

summ <- data.frame(
  setting   = c(seq_len(nsettings), NA),
  surf_type = c(stype, "ref_stationary"),
  moment    = c(settings.nonsta$moment, "none"),
  route     = c(settings.nonsta$route,  "none"),
  # spatial MEAN structure carried by the truth = what MRTS could ever recover
  truemean_sd = c(sapply(truemean.field, function(m) sd(m[, 1])), 0),
  y_sd        = c(sapply(seq_len(nsettings), function(k) sd(y[, , , k])), sd(y.ref)),
  # variance non-stationarity
  pw_sd_ratio = apply(pw.sd, 2, function(v) max(v) / min(v)),
  # dependence non-stationarity: spread of the quadrant variograms at lag ~1.5
  vgm_quad_spread = sapply(qv, function(g) {
    v <- g[lag_i, ]; (max(v) - min(v)) / mean(v)
  }),
  # the low-rank / Tzeng-vs-Matern signature
  min_correlation = apply(pc, 2, min),
  # higher moments
  skew_range = apply(pw.skew, 2, function(v) diff(range(v))),
  kurt_range = apply(pw.kurt, 2, function(v) diff(range(v))),
  # do the shape fields land where the truth says they should?
  cor_skew_g      = apply(pw.skew, 2, function(v) cor(v, g.field)),
  cor_skew_lambda = apply(pw.skew, 2, function(v) cor(v, lambda.field)),
  cor_kurt_h      = apply(pw.kurt, 2, function(v) cor(v, h.field)),
  stringsAsFactors = FALSE
)
write.csv(summ, csv_fp, row.names = FALSE)

vgm_ref <- summ$vgm_quad_spread[NREF]
summ$vgm_vs_ref <- summ$vgm_quad_spread / vgm_ref     # >> 1 = real non-stationarity

cat("\n============ DEPENDENCE / VARIANCE ============\n")
print(format(summ[, c("setting", "surf_type", "moment", "truemean_sd",
                      "pw_sd_ratio", "vgm_quad_spread", "vgm_vs_ref",
                      "min_correlation")], digits = 3), row.names = FALSE)
cat("\n============ HIGHER MOMENTS ============\n")
print(format(summ[, c("setting", "surf_type", "skew_range", "kurt_range",
                      "cor_skew_g", "cor_skew_lambda", "cor_kurt_h")], digits = 3),
      row.names = FALSE)

cat(sprintf("\nreference noise floor: quadrant-variogram spread = %.3f, min correlation = %.3f\n",
            vgm_ref, summ$min_correlation[NREF]))
cat(sprintf("wrote %s\nwrote %s\n", pdf_fp, csv_fp))
