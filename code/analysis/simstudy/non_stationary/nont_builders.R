#########################################################################
# Builders for the NON-t settings 13-21 of setup_nonsta.R: nine
# nonlinear mean surfaces plus the non-Gaussian error transform.
#
# Every surface builder returns the RAW ns-vector; setup_nonsta.R pushes it
# through ortho_std() so that (a) the component the K=0 design [1, s1, s2]
# already spans is removed (otherwise the oracle Brier headroom overstates
# what MRTS can add) and (b) the remainder is scaled to a common sd.
# Composite surfaces standardise each raw term to unit sd BEFORE summing,
# so no term dominates by accident of units.
#
# References:
#   Franke (1979), 'A critical comparison of some methods for interpolation
#     of scattered data' -- franke_fn()
#   Ricker wavelet (Mexican hat) -- annulus_fn()
#   Jones & Pewsey (2009), Biometrika -- sas_transform(), sas_moments()
#########################################################################

# -----------------------------------------------------------------------
# ortho_std(): remove the [1, s1, s2] projection (OLS on the ns sites),
# then scale the residual to target_sd. The result is exactly the part of
# the surface the K=0 baseline design cannot represent.
# -----------------------------------------------------------------------
ortho_std <- function(g_raw, s, target_sd = 11) {
  r <- stats::resid(stats::lm(g_raw ~ s[, 1] + s[, 2]))
  as.vector(target_sd * r / stats::sd(r))
}

# unit-sd standardisation for composite terms (mean removed for symmetry,
# though ortho_std() would absorb it anyway)
unit_sd <- function(v) (v - mean(v)) / stats::sd(v)

# -----------------------------------------------------------------------
# Setting 13/18: Franke's function on [0,1]^2. The canonical TPS test
# surface: two broad bumps, one medium bump, one small dip -- multiscale by
# construction. dip_widen (default 1) scales the width of the sharpest term
# (the dip, sd ~ 0.79 units on [0,10]); the precheck decides whether it
# needs widening to stay inside the site-density resolution cap.
# -----------------------------------------------------------------------
franke_fn <- function(s01, dip_widen = 1) {
  x <- s01[, 1]; y <- s01[, 2]
  0.75 * exp(-((9 * x - 2)^2 + (9 * y - 2)^2) / 4) +
  0.75 * exp(-(9 * x + 1)^2 / 49 - (9 * y + 1) / 10) +
  0.50 * exp(-((9 * x - 7)^2 + (9 * y - 3)^2) / 4) -
  0.20 * exp(-((9 * x - 4)^2 + (9 * y - 7)^2) / dip_widen^2)
}

# -----------------------------------------------------------------------
# Setting 14: curved (banana) ridge. A Gaussian ridge whose crest follows
# y01 = 0.5 + 0.25 sin(2 pi x01): the local orientation of the feature
# rotates across the domain, which no single radial bump can mimic.
# -----------------------------------------------------------------------
ridge_fn <- function(s01, width = 0.12) {
  exp(-(s01[, 2] - 0.5 - 0.25 * sin(2 * pi * s01[, 1]))^2 / (2 * width^2))
}

# -----------------------------------------------------------------------
# Setting 15: Ricker (Mexican-hat) annulus in the radial coordinate.
# psi(u) = (1 - u^2) exp(-u^2 / 2) applied to u = (r - radius)/width:
# positive on the ring r = radius, negative side lobes inside and outside
# -- a non-monotone radial mean, on [0,10] coordinates.
# -----------------------------------------------------------------------
annulus_fn <- function(s, center = c(5, 5), radius = 3, width = 1.2) {
  r <- sqrt((s[, 1] - center[1])^2 + (s[, 2] - center[2])^2)
  u <- (r - radius) / width
  (1 - u^2) * exp(-u^2 / 2)
}

# -----------------------------------------------------------------------
# Setting 16: sigmoid front -- a smoothed step along a circular arc, the
# steepest structure the site density plausibly supports (transition width
# ~ 2*w units). On [0,10] coordinates.
# -----------------------------------------------------------------------
front_fn <- function(s, center = c(2, 2), radius = 4.5, width = 1.0) {
  r <- sqrt((s[, 1] - center[1])^2 + (s[, 2] - center[2])^2)
  tanh((r - radius) / width)
}

# -----------------------------------------------------------------------
# Setting 19: full second-order polynomial surface = bowl + saddle,
# the textbook regression nonlinearity (quadratic main effects plus the
# interaction term).
# -----------------------------------------------------------------------
poly2_fn <- function(s01) {
  bowl   <- (s01[, 1] - 0.5)^2 + (s01[, 2] - 0.5)^2
  saddle <- (s01[, 1] - 0.5) * (s01[, 2] - 0.5)
  unit_sd(bowl) + unit_sd(saddle)
}

# -----------------------------------------------------------------------
# Setting 20: monotone-transform mix -- log, exponential decay, and a
# ratio-type interaction, i.e. the nonlinear terms an applied modeller
# actually writes down. Curvature concentrates near the domain boundary,
# where training sites are sparsest. On [0,10] coordinates.
# -----------------------------------------------------------------------
transform_mix_fn <- function(s) {
  unit_sd(log(1 + s[, 1])) +
  unit_sd(exp(-s[, 1] / 2)) +
  unit_sd(s[, 1] / (1 + s[, 2]))
}

# -----------------------------------------------------------------------
# Setting 21: non-smooth terms -- a hinge (broken stick) in s1 and a
# V-shape in s2, kink locations deliberately offset. C^0 creases cannot be
# represented exactly by the C^inf TPS eigenfunctions, so the projection
# residual is expected to PLATEAU above zero: the mean-side counterpart of
# the dependence-side "structurally cannot" settings 3-8.
# -----------------------------------------------------------------------
nonsmooth_fn <- function(s, hinge_at = 5, vee_at = 4) {
  unit_sd(pmax(0, s[, 1] - hinge_at)) +
  unit_sd(abs(s[, 2] - vee_at))
}

# -----------------------------------------------------------------------
# Setting 17: Matern correlation for the random mean surface
# u(s) ~ GP(0, Matern(range, nu = 1.5)). Deliberately SMOOTHER (nu = 1.5)
# and longer-ranged than the error's exponential, so the mean and the error
# live at different smoothness classes. Returns the correlation matrix.
# -----------------------------------------------------------------------
matern_cor <- function(s, range = 2, smoothness = 1.5) {
  d <- as.matrix(stats::dist(s))
  C <- fields::Matern(d, range = range, smoothness = smoothness)
  C <- (C + t(C)) / 2
  diag(C) <- 1
  C
}

# -----------------------------------------------------------------------
# sas_transform(): sinh-arcsinh transform of Jones & Pewsey (2009),
#   X = sinh(delta * asinh(Z) + eps),  Z ~ N(0,1).
# eps controls skewness (eps > 0 -> right skew), delta the tail weight
# (delta = 1 keeps tails LIGHTER than or equal to Gaussian-ish -- crucially,
# all moments exist for every eps, unlike the g-and-h). With delta = 1 this
# is X = sinh(asinh(Z) + eps): skewed, light-tailed, exactly what "non-t
# marginal" needs.
# -----------------------------------------------------------------------
sas_transform <- function(z, eps = 0.8, delta = 1) {
  sinh(delta * asinh(z) + eps)
}

# -----------------------------------------------------------------------
# sas_moments(): mean and sd of sas_transform(Z), Z ~ N(0,1), by Monte
# Carlo (same pattern as gh_moments() in ns_cov_builders.R: one large
# reusable sample, RNG state saved and restored). Closed forms exist but
# the MC keeps the standardisation consistent with the g-and-h setting.
# -----------------------------------------------------------------------
sas_moments <- function(eps = 0.8, delta = 1, nsim = 2e5, seed = 20240102) {
  old <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
  set.seed(seed)
  z <- stats::rnorm(nsim)
  if (!is.null(old)) assign(".Random.seed", old, .GlobalEnv)

  x <- sas_transform(z, eps, delta)
  list(mean = mean(x), sd = stats::sd(x))
}
