#########################################################################
# Builders for the non-stationary structures used by setup_nonsta.R
# (settings 4-10) and by the diagnostics script.
#
# Every covariance builder here returns a raw ns x ns covariance; run it
# through to_correlation() before handing it to rpotspatTS(), which expects
# a unit-diagonal correlation matrix (see CorFx / CorFxDef / CorFxNL in
# R/ar2/auxfunctions.R -- all three end with `diag(cor) <- 1`).
#
# References:
#   Paciorek & Schervish (2006), Environmetrics -- ps_cov(), ps_cov_iso()
#   Schmidt, Guttorp & O'Hagan (2011); Risser & Calder (2015) -- ps_cov_iso()
#   Fuentes (2001, 2002) -- fuentes_cov()
#   Xu & Genton (2017), JASA -- gh_transform(), gh_moments()
#########################################################################

# -----------------------------------------------------------------------
# to_correlation(): raw covariance -> correlation, apply the gamma nugget,
# force a unit diagonal, and assert positive definiteness.
#
# gamma has the same meaning as in CorFx(): it scales every entry and then
# the diagonal is put back to 1, i.e. gamma acts as (1 - nugget).
# -----------------------------------------------------------------------
to_correlation <- function(C, gamma, tol = 1e-8, label = "C") {
  if (!is.matrix(C) || nrow(C) != ncol(C)) {
    stop(sprintf("%s must be a square matrix", label))
  }
  C <- (C + t(C)) / 2                 # kill numerical asymmetry
  cor <- stats::cov2cor(C)
  cor <- gamma * cor
  diag(cor) <- 1

  ev <- min(eigen(cor, symmetric = TRUE, only.values = TRUE)$values)
  if (ev <= tol) {
    stop(sprintf("%s is not positive definite (min eigenvalue = %.3e)", label, ev))
  }
  attr(cor, "min_eigenvalue") <- ev
  cor
}

# -----------------------------------------------------------------------
# ps_kernel_aniso(): per-site 2x2 kernel matrices for settings 4/5.
#
#   Sigma(s) = R(theta(s)) diag(a1(s)^2, a2(s)^2) R(theta(s))'
#   a1(s) = rho0 * exp( kappa_scale*(s2/10 - 0.5) + kappa_aspect*(s1/10 - 0.5))
#   a2(s) = rho0 * exp( kappa_scale*(s2/10 - 0.5) - kappa_aspect*(s1/10 - 0.5))
#   theta(s) = (pi/2) * (s2/10)
#
# Three things vary independently across the domain:
#   SIZE        sqrt|Sigma(s)| = rho0^2 exp(2 kappa_scale (s2/10 - 0.5))  -- with s2
#   ASPECT      a1/a2 = exp(2 kappa_aspect (s1/10 - 0.5))                 -- with s1
#   ORIENTATION theta(s)                                                  -- with s2
#
# The SIZE gradient matters and is easy to get wrong: with kappa_scale = 0 the
# determinant is constant, the correlation ellipse only changes shape, and an
# OMNIDIRECTIONAL variogram (or an isotropic fitted model) averages the
# anisotropy away and sees a stationary process. The scale gradient is what
# makes the non-stationarity visible without conditioning on direction.
#
# Returns a list of ns 2x2 matrices.
# -----------------------------------------------------------------------
ps_kernel_aniso <- function(s, rho0 = 1, kappa_scale = 0.8, kappa_aspect = 1,
                            domain = 10) {
  u1 <- s[, 1] / domain - 0.5
  u2 <- s[, 2] / domain - 0.5
  a1 <- rho0 * exp(kappa_scale * u2 + kappa_aspect * u1)
  a2 <- rho0 * exp(kappa_scale * u2 - kappa_aspect * u1)
  th <- (pi / 2) * (s[, 2] / domain)

  lapply(seq_len(nrow(s)), function(i) {
    R <- matrix(c(cos(th[i]), sin(th[i]), -sin(th[i]), cos(th[i])), 2, 2)
    R %*% diag(c(a1[i]^2, a2[i]^2)) %*% t(R)
  })
}

# -----------------------------------------------------------------------
# ps_cov(): the general Paciorek-Schervish non-stationary covariance.
#
#   C(s,s') = |Sig(s)|^{1/4} |Sig(s')|^{1/4} |Sigbar|^{-1/2} * M_nu(sqrt(Q))
#   Sigbar  = (Sig(s) + Sig(s')) / 2
#   Q       = (s-s')' Sigbar^{-1} (s-s')
#
# The normalising constants make C(s,s) = 1 exactly, so the marginal variance
# is untouched and only the DEPENDENCE becomes non-stationary. nu = 0.5 gives
# the exponential, matching the rest of the simulation study.
#
# ns = 144 here, so the O(ns^2) double loop with a 2x2 solve per pair is
# instant; no need for anything cleverer.
# -----------------------------------------------------------------------
ps_cov <- function(s, Sigma_list, nu = 0.5) {
  if (!isTRUE(all.equal(nu, 0.5))) {
    stop("ps_cov() only implements nu = 0.5 (exponential); the study uses nu = 0.5 throughout")
  }
  ns <- nrow(s)
  det4 <- vapply(Sigma_list, function(S) det(S)^(1 / 4), numeric(1))

  C <- matrix(1, ns, ns)
  for (i in seq_len(ns - 1)) {
    Si <- Sigma_list[[i]]
    for (j in (i + 1):ns) {
      Sbar <- (Si + Sigma_list[[j]]) / 2
      dij  <- s[i, ] - s[j, ]
      Q    <- as.numeric(crossprod(dij, solve(Sbar, dij)))
      val  <- det4[i] * det4[j] / sqrt(det(Sbar)) * exp(-sqrt(Q))
      C[i, j] <- val
      C[j, i] <- val
    }
  }
  C
}

# -----------------------------------------------------------------------
# ps_cov_iso(): the isotropic-kernel special case, Sigma(s) = rho(s)^2 I_2,
# used by settings 6/7 (covariate-driven range).
#
# Substituting into the PS formula and simplifying:
#
#   C(s,s') = 2 rho(s) rho(s') / (rho(s)^2 + rho(s')^2)
#             * exp( -||s-s'|| * sqrt( 2 / (rho(s)^2 + rho(s')^2) ) )
#
# Unit diagonal is immediate (set rho(s') = rho(s)). Positive definiteness is
# inherited from the general PS class. Fully vectorised.
# -----------------------------------------------------------------------
ps_cov_iso <- function(s, rho_field) {
  ns <- nrow(s)
  stopifnot(length(rho_field) == ns, all(rho_field > 0))

  r2  <- rho_field^2
  ssq <- outer(r2, r2, "+")                       # rho(s)^2 + rho(s')^2
  pre <- 2 * outer(rho_field, rho_field, "*") / ssq
  d   <- as.matrix(stats::dist(s))

  pre * exp(-d * sqrt(2 / ssq))
}

# -----------------------------------------------------------------------
# fuentes_cov(): Fuentes (2001) weighted mixture of stationary GPs.
#
#   Y(s) = sum_k w_k(s) Z_k(s),  Z_k stationary exponential with range rho_k
#   C(s,s') = sum_k w_k(s) w_k(s') exp(-||s-s'|| / rho_k)
#
# Weights are Gaussian kernels around the mixture centres, normalised in L2
# (sum_k w_k(s)^2 = 1) so that the diagonal is exactly 1.
#
# C = sum_k diag(w_k) C_k diag(w_k) is a sum of PSD matrices, hence PSD; the
# gamma nugget in to_correlation() makes it strictly PD.
#
# Unlike PS (a smooth gradient), this gives REGIME-like behaviour: one corner
# of the domain is rough, the opposite corner smooth.
# -----------------------------------------------------------------------
fuentes_cov <- function(s, centers, ranges, h) {
  ns <- nrow(s)
  K  <- nrow(centers)
  stopifnot(length(ranges) == K, h > 0)

  w <- sapply(seq_len(K), function(k) {
    d2 <- rowSums(sweep(s, 2, centers[k, ], "-")^2)
    exp(-d2 / (2 * h^2))
  })                                              # ns x K
  w <- w / sqrt(rowSums(w^2))                     # L2-normalise -> unit diagonal

  d <- as.matrix(stats::dist(s))
  C <- matrix(0, ns, ns)
  for (k in seq_len(K)) {
    C <- C + outer(w[, k], w[, k], "*") * exp(-d / ranges[k])
  }
  attr(C, "weights") <- w
  C
}

# -----------------------------------------------------------------------
# gh_transform(): the Tukey g-and-h transform (Xu & Genton 2017).
#
#   tau_{g,h}(z) = (exp(g z) - 1) / g * exp(h z^2 / 2)     (g != 0)
#   tau_{0,h}(z) = z * exp(h z^2 / 2)                      (g  = 0)
#
# g controls skewness (sign of g = direction of skew), h controls tail weight
# (h = 0 recovers the Gaussian). Moments of order r exist iff h < 1/r, so
# h < 0.5 is required for a finite variance.
#
# z may be a vector or an ns x nt matrix; g and h may be scalars or per-site
# vectors (recycled down the columns of z, i.e. one (g, h) per site).
#
# g and h are broadcast to z's shape FIRST: ifelse() takes its result shape
# from the test, so a scalar g would otherwise collapse the output to length 1.
# expm1 keeps the g -> 0 limit (expm1(gz)/g -> z) numerically clean.
# -----------------------------------------------------------------------
gh_transform <- function(z, g, h) {
  gg <- z * 0 + g                     # broadcast to dim(z) via arithmetic recycling
  hh <- z * 0 + h
  near0 <- abs(gg) < 1e-10
  core <- ifelse(near0, z, expm1(gg * z) / ifelse(near0, 1, gg))
  core * exp(hh * z^2 / 2)
}

# -----------------------------------------------------------------------
# gh_moments(): per-site mean and sd of tau_{g(s),h(s)}(Z), Z ~ N(0,1).
#
# Closed forms exist but are fiddly; a single large Monte Carlo sample reused
# across sites is simpler, accurate enough, and keeps the code honest. The
# SAME z sample is used for every site, so the standardisation is consistent
# across the domain.
#
# Returns a list(mean = ns-vector, sd = ns-vector).
# -----------------------------------------------------------------------
gh_moments <- function(g, h, nsim = 2e5, seed = 20240101) {
  stopifnot(length(g) == length(h), all(h < 0.5))   # h >= 0.5 -> infinite variance

  old <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
  set.seed(seed)
  z <- stats::rnorm(nsim)
  if (!is.null(old)) assign(".Random.seed", old, .GlobalEnv)

  m <- numeric(length(g))
  v <- numeric(length(g))
  for (i in seq_along(g)) {
    ti   <- gh_transform(z, g[i], h[i])
    m[i] <- mean(ti)
    v[i] <- stats::sd(ti)
  }
  list(mean = m, sd = v)
}
