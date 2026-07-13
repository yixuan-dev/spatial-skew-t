#########################################################################
# Non-stationary simulation study.  All 10 settings share the same skew-t-1
# core (lambda = 3, dist = "t", nknots = 1); they differ in WHERE the
# non-stationarity lives and by WHICH mechanism it is generated.
#
# The point of the grid is to turn a single data point ("MRTS cannot recover
# the setting-3 random effect") into a line:
#
#   non-stationarity in the 1st moment (MEAN)        -> MRTS recovers it
#   non-stationarity in the 2nd moment (DEPENDENCE)  -> MRTS structurally cannot,
#                                                       whatever the mechanism
#   non-stationarity in higher moments               -> depends on whether an
#                                                       induced mean appears
#
#  ID  surf_type        mechanism                          route       moment
#  --  ---------------  ---------------------------------  ----------  ---------
#   1  invariant        fixed cosine-bump mean surface     mean        1st
#   2  varying          time-varying cosine-bump mean      mean        1st
#   3  ns_dependence    low-rank cosine random effect      additive    2nd
#   4  ps_cov           Paciorek-Schervish local aniso.    replace C   2nd
#   5  ps_add           same PS field, additive            additive    2nd + marg
#   6  covdep_cov       rho(s) = rho0 exp(beta f1(s))      replace C   2nd
#   7  covdep_add       same field, additive               additive    2nd + marg
#   8  fuentes_cov      4-regime weighted GP mixture       replace C   2nd
#   9  gh_marginal      Tukey g-and-h, g(s) and h(s)       additive    3rd + 4th
#  10  lambda_varying   lambda(s) = 3 + 3 f2(s)            skew term   3rd
#
# THE HEADLINE PAIR IS 1 vs 6.  Both are driven by the SAME spatial pattern
# f1(s).  Setting 1 puts f1 in the mean, where MRTS recovers it exactly;
# setting 6 puts f1 in the correlation range, where no amount of K_MRTS helps.
# That rules out "the model lacked the information" and isolates the real
# cause: mean-only augmentation is looking at the wrong moment.
#
# Settings 4/5 and 6/7 are A/B pairs: the same non-stationary structure, once
# as the error correlation (marginal skew-t preserved, dependence changed) and
# once as an additive mean-zero random effect (marginal changed too).  Each
# pair SHARES A BASE SEED, so the underlying skew-t draw (tau, z, knots and the
# standard-normal innovations) is identical and the injection route is the only
# difference.
#
# BACKWARD COMPATIBILITY: settings 1-3 must reproduce bit-for-bit, because the
# fits already in results_nonsta/ are keyed by setting id.  Their code path and
# seeds are untouched.
#
# References: see non_stationary/nonsta_settings_design.md
#########################################################################

rm(list = ls())

library(fields)
library(SpatialTools)

# Exposes rpotspatTS(), CorFx(), mem(), etc.
source("../../R/ar2/auxfunctions.R")
# Exposes ps_cov(), ps_cov_iso(), fuentes_cov(), gh_transform(), to_correlation()
source("non_stationary/ns_cov_builders.R")

# Pure-R override of mem() to avoid the Rcpp build of g.Rcpp.
# With nknots = 1 every site maps to knot 1, but keep the generic form.
mem <- function(s, knots) {
  d <- fields::rdist(s, knots)
  apply(d, 1, which.min)
}

# -----------------------------------------------------------------------
# Common skew-t parameters (match setup.R / setup_def.R)
# -----------------------------------------------------------------------
beta.t      <- c(10, 0, 0)  # intercept only; coords carry NO linear trend
nu.t        <- 0.5
gamma.t     <- 0.9
rho.t       <- 1
tau.alpha.t <- 3
tau.beta.t  <- 8

# All settings are skew-t-1
surf.type <- c("invariant", "varying", "ns_dependence",
               "ps_cov", "ps_add", "covdep_cov", "covdep_add",
               "fuentes_cov", "gh_marginal", "lambda_varying")
nsettings <- length(surf.type)

dist.nonsta   <- rep("t", nsettings)
nknots.nonsta <- rep(1L,  nsettings)
lambda.nonsta <- rep(3,   nsettings)   # setting 10 varies lambda(s) around this

# Which moment each setting perturbs, and how it is injected.
moment.nonsta <- c("mean", "mean", "dependence",
                   "dependence", "dependence", "dependence", "dependence",
                   "dependence", "tails", "skewness")
route.nonsta  <- c("mean", "mean", "additive",
                   "replace_C", "additive", "replace_C", "additive",
                   "replace_C", "additive", "skew_term")
ref.nonsta    <- c("Tzeng & Huang (2018)", "Tzeng & Huang (2018)",
                   "Cressie & Johannesson (2008)",
                   "Paciorek & Schervish (2006)", "Paciorek & Schervish (2006)",
                   "Schmidt et al. (2011); Risser & Calder (2015)",
                   "Schmidt et al. (2011); Risser & Calder (2015)",
                   "Fuentes (2001, 2002)", "Xu & Genton (2017)",
                   "Morris et al. (2017) assumption violated")

# A/B pairs share a base seed so the underlying skew-t draw is identical:
# setting 5 borrows setting 4's, setting 7 borrows setting 6's.  Settings 1-3
# map to themselves, which reproduces the original seeds exactly.
pair.seed <- c(1, 2, 3, 4, 4, 6, 6, 8, 9, 10)

# -----------------------------------------------------------------------
# Non-stationary surface controls (Tzeng & Huang 2018, Scenario 1)
#   f1(s) = cos(   pi * || s01 - c1 ||)
#   f2(s) = cos( 2*pi * || s01 - c2 ||)
# Sites live on [0,10]^2; Tzeng's bumps are defined on [0,1]^2, so we
# evaluate the basis on the rescaled coordinates s01 = s / 10.
#
# surf.scale globally scales the surface amplitude relative to the skew-t
# noise (sd ~ 1/sqrt(tau), tau ~ Gamma(3/2, 4) -> sd ~ 1.6).  Increase it
# to make the spatial structure easier to detect, decrease for a harder
# (lower-SNR) recovery problem.
# -----------------------------------------------------------------------
c1        <- c(0,    1)      # f1 centre (in [0,1]^2 coordinates)
c2        <- c(0.75, 0.25)   # f2 centre
surf.scale <- 1.0
a.fixed   <- c(5, 3) * surf.scale   # time-invariant coefficients (setting 1)
M.var     <- diag(c(25, 9)) * surf.scale^2  # Var(w1, w2) (setting 2)

# -----------------------------------------------------------------------
# Sites & covariates (same RNG seed and grid as setup.R / setup_def.R)
# -----------------------------------------------------------------------
set.seed(20)
s     <- cbind(runif(144, 0, 10), runif(144, 0, 10))
ns    <- nrow(s)
nt    <- 50
nsets <- 50
ntest <- 44

x <- array(1, c(ns, nt, 3))
for (t in 1:nt) {
  x[, t, 2] <- s[, 1]
  x[, t, 3] <- s[, 2]
}

# Evaluate the true basis fields once (s is fixed across settings/sets)
s01 <- s / 10
f1  <- cos(    pi * sqrt((s01[, 1] - c1[1])^2 + (s01[, 2] - c1[2])^2))
f2  <- cos(2 * pi * sqrt((s01[, 1] - c2[1])^2 + (s01[, 2] - c2[2])^2))
f.basis <- cbind(f1 = f1, f2 = f2)   # ns x 2  (the ground-truth surface basis)

# -----------------------------------------------------------------------
# Setting 3 control: NON-STATIONARY DEPENDENCE injected as a mean-zero,
# per-day RANDOM EFFECT (basis-function random effect / FRK-style):
#   u_t(s) = sum_j xi_{tj} * cos(kappa_j * ||s - c_j||),  xi_{tj} ~ N(0, tau_j^2),
# redrawn every day t.  Marginalizing over xi gives a low-rank NON-STATIONARY
# covariance Cov(u_t(s), u_t(s')) = sum_j tau_j^2 phi_j(s) phi_j(s'); the
# variance sum_j tau_j^2 phi_j(s)^2 varies with location -> non-stationary.
#
# This is the SAME cosine basis as Tzeng & Huang (2018, Scenario 1) but written
# NATIVELY on the [0,10]^2 domain (no s/10 rescaling): the wavenumber kappa_j is
# specified directly in [0,10] units.  Algebraically identical to the rescaled
# form, since cos(a_j pi ||s/10 - c_j/10||) = cos((a_j pi / 10) ||s - c_j||),
# so kappa_j = a_j pi / 10 with the centre c_j placed in [0,10]^2.
#
# Purpose: a data-generating process whose non-stationarity lives in the
# DEPENDENCE (second moment / random effect), NOT the mean.  The MRTS extension
# acts on the fixed-effect mean with a coefficient pooled over time, so it
# CANNOT capture this structure: the score-vs-K (and recovery) curve should stay
# flat / rise, demonstrating the structural limit of mean-only augmentation.
#
# It is also the one setting whose correlation goes NEGATIVE at long range: a
# low-rank cosine basis can do that, a monotone Matern cannot.
# -----------------------------------------------------------------------
kappa.re   <- c(pi / 10, pi / 5)            # native wavenumbers on [0,10]^2
centers.re <- rbind(c(0, 10), c(7.5, 2.5))  # bump centres in [0,10]^2
tau.re     <- c(5, 3)                       # component sds (Tzeng: tau^2 = 25, 9)
F.re <- sapply(seq_len(nrow(centers.re)), function(j) {
  cos(kappa.re[j] * fields::rdist(s, centers.re[j, , drop = FALSE]))
})                                          # ns x 2 (native cosine basis)
Tcov.re <- diag(tau.re^2)                   # Var(xi_t) = diag(25, 9)
L.re    <- chol(Tcov.re)                    # Tcov.re = L'L

# -----------------------------------------------------------------------
# The stationary reference correlation: the exponential the skew-t error uses
# by default.  Settings 5/7/9/10 keep it (their structure is injected some
# other way), and the diagnostics use it as the "no non-stationarity" baseline.
# -----------------------------------------------------------------------
C.stat <- CorFx(d = as.matrix(dist(s)), gamma = gamma.t, rho = rho.t, nu = nu.t)

# -----------------------------------------------------------------------
# Settings 4/5: Paciorek-Schervish local anisotropy.
# The correlation ellipse changes SIZE with s2 (kappa_scale), changes ASPECT
# with s1 (kappa_aspect), and ROTATES with s2 -- three independent gradients.
# The PS normalising constants give C(s,s) = 1 exactly, so replacing the error
# correlation with C.ps leaves the skew-t MARGINAL untouched and changes only
# the DEPENDENCE.
#
# kappa_scale must be non-zero.  With aspect and rotation alone the determinant
# |Sigma(s)| is constant, so the ellipse only changes shape; an omnidirectional
# variogram -- and, more to the point, an ISOTROPIC fitted model -- averages
# that away and sees a stationary process.  An earlier version of this setting
# had exactly that bug and registered no non-stationarity at all.
# -----------------------------------------------------------------------
ps.par   <- list(rho0 = rho.t, kappa_scale = 0.8, kappa_aspect = 1.0)
Sigma.ps <- ps_kernel_aniso(s, rho0 = ps.par$rho0,
                            kappa_scale = ps.par$kappa_scale,
                            kappa_aspect = ps.par$kappa_aspect)
C.ps     <- to_correlation(ps_cov(s, Sigma.ps, nu = nu.t), gamma.t, label = "C.ps")

# -----------------------------------------------------------------------
# Settings 6/7: covariate-driven dependence.  THE HEADLINE CONTRAST.
# The SAME f1 that drives the mean surface in setting 1 now drives the local
# correlation range instead:  rho(s) = rho0 * exp(beta_rho * f1(s)).
# f1 is a smooth function of location, so the MRTS basis spans it (setting 1
# proves exactly that) -- yet here it is unrecoverable, because the truth put it
# in the second moment.  Isotropic PS kernel Sigma(s) = rho(s)^2 I, closed form.
# -----------------------------------------------------------------------
covdep.par <- list(rho0 = rho.t, beta_rho = 0.7)
rho.field  <- covdep.par$rho0 * exp(covdep.par$beta_rho * f1)   # in [0.50, 2.01]
C.covdep   <- to_correlation(ps_cov_iso(s, rho.field), gamma.t, label = "C.covdep")

# -----------------------------------------------------------------------
# Setting 8: Fuentes weighted mixture of 4 stationary GPs whose ranges span a
# factor of 10.  Where PS gives a smooth gradient, this gives REGIMES: the
# south-west corner is rough (range 0.4), the north-east smooth (range 4).
# -----------------------------------------------------------------------
fuentes.par <- list(
  centers = rbind(c(2.5, 2.5), c(7.5, 2.5), c(2.5, 7.5), c(7.5, 7.5)),
  ranges  = c(0.4, 1.0, 2.0, 4.0),
  h       = 3
)
C.fuentes <- to_correlation(
  fuentes_cov(s, fuentes.par$centers, fuentes.par$ranges, fuentes.par$h),
  gamma.t, label = "C.fuentes"
)

# -----------------------------------------------------------------------
# Setting 9: Tukey g-and-h random field (Xu & Genton 2017).  A stationary
# Gaussian field is pushed through a POINTWISE transform whose parameters vary
# in space: g(s) sets the skewness (it CHANGES SIGN across the domain) and h(s)
# sets the tail weight (Gaussian in the west, heavy-tailed in the east).
#
# h < 0.5 is required for a finite variance; we cap at 0.24 so the fourth moment
# stays finite too -- otherwise a handful of enormous draws dominate the field
# and the SNR calibration stops meaning anything.
#
# The field is standardised pointwise to mean 0 / variance 1, so the mean and
# the variance stay stationary and ONLY the higher moments (and, mildly, the
# correlation) become non-stationary.  Mean-zero => nothing for a pooled fixed
# mean to recover.
#
# g is driven by f1 and h by f2, DELIBERATELY the opposite way round from
# setting 10's lambda(s) = 3 + 3*f2.  If g also rode on f2 the two settings
# would carry the same spatial pattern in their higher moments and no
# diagnostic could tell them apart -- cor(skewness, g) and cor(skewness, lambda)
# would be identical by construction.
# -----------------------------------------------------------------------
gh.par  <- list(g_scale = 0.8, h_scale = 0.12)
g.field <- gh.par$g_scale * f1            # in [-0.80, 0.80]  (sign changes)
h.field <- gh.par$h_scale * (1 + f2)      # in [ 0.00, 0.24]
stopifnot("g-and-h needs h < 0.5 for a finite variance" = max(h.field) < 0.5)
gh.mom  <- gh_moments(g.field, h.field, nsim = 2e5)

# -----------------------------------------------------------------------
# Setting 10: spatially varying skewness.  The fitted model HAS a lambda -- it
# just assumes it is a constant scalar (Morris et al. 2017).  So this is a
# misspecification study of the model's own assumption.
#
# Because nknots = 1, z_t is a single scalar per day and mu_t = x'beta + lambda*z_t,
# so swapping the constant lambda for a field is a pure post-hoc adjustment with
# no change to the generator:
#   y_new[, t] = y_old[, t] + (lambda(s) - lambda0) * z_t
#
# This is the "partial capture" setting: E[Y(s)] = 10 + lambda(s) E[z_t], so
# lambda(s) INDUCES a spatially varying mean that MRTS can pick up, even though
# the heterogeneous skewness of the predictive distribution cannot be.
# -----------------------------------------------------------------------
lambda.par   <- list(lambda0 = 3, lambda1 = 3)
lambda.field <- lambda.par$lambda0 + lambda.par$lambda1 * f2   # in [0, 6]

# -----------------------------------------------------------------------
# Amplitude of the additive mean-zero random effects (settings 5, 7, 9).
# Chosen so their sd is comparable to setting 3's random effect
# (sd = sqrt(sum_j tau_j^2 phi_j(s)^2) ~ 4) and a couple of times the skew-t
# noise sd (~1.6): strong enough to matter, not so strong it swamps everything.
# -----------------------------------------------------------------------
sigma.add <- 3

L.stat   <- t(chol(C.stat))
L.ps     <- t(chol(C.ps))
L.covdep <- t(chol(C.covdep))

# -----------------------------------------------------------------------
# Storage
# -----------------------------------------------------------------------
y       <- array(NA, dim = c(ns, nt, nsets, nsettings))
tau.t   <- vector("list", length = nsettings)
z.t     <- vector("list", length = nsettings)
knots.t <- vector("list", length = nsettings)

# Save the surface building blocks so the true mean can be reconstructed:
#   setting 1: g(s)   = f.basis %*% a.fixed         (same for all sets)
#   setting 2: g_t(s) = f.basis %*% t(W[[set]])     (nt x 2 weights)
surf.coef.invariant <- a.fixed
W.varying           <- vector("list", length = nsets)  # filled below
#   setting 3: u_t(s) = F.re %*% t(W3),  W3 (nt x 2) ~ N(0, diag(25, 9)) per day
W.nsdep             <- vector("list", length = nsets)  # filled below

# THE UNIFIED TRUTH CONTRACT.  truemean.field[[setting]][, set] is the
# time-averaged true mean at every site, for that setting and that dataset.
# scores.R reads this directly instead of re-deriving the mean per surf_type --
# the old if/else chain silently fell through to the "varying" formula for any
# surf_type it did not recognise, which would have produced plausible-looking
# but wrong recovery numbers for every setting added here.
truemean.field <- lapply(seq_len(nsettings), function(i) matrix(NA_real_, ns, nsets))

L.var <- chol(M.var)  # M.var = L'L

for (setting in seq_len(nsettings)) {
  stype  <- surf.type[setting]
  nknots <- nknots.nonsta[setting]
  tau.t.setting   <- array(NA, dim = c(nknots, nt, nsets))
  z.t.setting     <- array(NA, dim = c(nknots, nt, nsets))
  knots.t.setting <- array(NA, dim = c(nknots, nt, 2, nsets))

  # Which correlation the skew-t error itself uses.  The three "replace_C"
  # settings hand rpotspatTS a precomputed non-stationary correlation; everyone
  # else keeps the stationary exponential and injects structure another way.
  C.err <- switch(stype,
    ps_cov      = C.ps,
    covdep_cov  = C.covdep,
    fuentes_cov = C.fuentes,
    NULL
  )
  cov.type.setting <- if (is.null(C.err)) "matern" else "precomputed"

  for (set in seq_len(nsets)) {
    # Reproducible skew-t draw.  A/B pairs share a base seed (pair.seed), so
    # settings 4/5 and 6/7 get an identical underlying draw and differ only in
    # how the non-stationarity is injected.  Settings 1-3 map to themselves,
    # reproducing the original seeds exactly.
    set.seed(pair.seed[setting] * 1000 + set)
    data <- rpotspatTS(
      nt        = nt,
      x         = x,
      s         = s,
      beta      = beta.t,
      gamma     = gamma.t,
      nu        = nu.t,
      rho       = rho.t,
      phi.z     = 0,
      phi.w     = 0,
      phi.tau   = 0,
      lambda    = lambda.nonsta[setting],
      tau.alpha = tau.alpha.t,
      tau.beta  = tau.beta.t,
      nknots    = nknots,
      dist      = dist.nonsta[setting],
      cov.type  = cov.type.setting,
      C.mat     = C.err
    )

    # Inject the non-stationary structure on top of the skew-t draw.  Because
    # the surface / random effect enters additively and is independent of the
    # skew-t error, adding it post-hoc here is identical to injecting it into
    # `mu` inside rpotspatTS().
    if (stype == "invariant") {
      g_s <- as.vector(f.basis %*% a.fixed)        # ns-vector, constant in t
      data$y <- data$y + matrix(g_s, ns, nt)
      truemean.field[[setting]][, set] <- 10 + g_s

    } else if (stype == "varying") {
      # Tzeng time-varying weights, drawn reproducibly per dataset.
      set.seed(900000 + set)
      W <- matrix(rnorm(nt * 2), nt, 2) %*% L.var  # nt x 2
      W.varying[[set]] <- W
      data$y <- data$y + f.basis %*% t(W)          # ns x nt
      truemean.field[[setting]][, set] <- 10 + as.vector(f.basis %*% colMeans(W))

    } else if (stype == "ns_dependence") {
      # u_t(s) = F.re %*% t(W3), W3 (nt x 2) ~ N(0, diag(25, 9)) redrawn per day.
      # Mean-zero in t -> no fixed mean structure for a pooled-coefficient MRTS
      # mean to recover; the non-stationarity is entirely in the dependence.
      set.seed(750000 + set)
      W3 <- matrix(rnorm(nt * 2), nt, 2) %*% L.re  # nt x 2
      W.nsdep[[set]] <- W3
      data$y <- data$y + F.re %*% t(W3)            # ns x nt
      truemean.field[[setting]][, set] <- 10 + as.vector(F.re %*% colMeans(W3))

    } else if (stype %in% c("ps_cov", "covdep_cov", "fuentes_cov")) {
      # The non-stationarity is already inside the error correlation rpotspatTS
      # just used.  The mean was never touched: it is flat at x'beta = 10
      # everywhere.  This is the cleanest possible form of the claim -- there is
      # NO mean structure at all, so no K_MRTS can help.
      truemean.field[[setting]][, set] <- 10

    } else if (stype %in% c("ps_add", "covdep_add")) {
      # The same non-stationary correlation, but as an additive mean-zero
      # Gaussian random effect on top of a STATIONARY skew-t draw.  The marginal
      # changes here (skew-t convolved with a Gaussian), unlike its replace_C
      # twin -- and that difference is exactly what the A/B pair isolates.
      set.seed(if (stype == "ps_add") 810000 + set else 820000 + set)
      L.add <- if (stype == "ps_add") L.ps else L.covdep
      u <- sigma.add * (L.add %*% matrix(rnorm(ns * nt), ns, nt))   # ns x nt
      data$y <- data$y + u
      truemean.field[[setting]][, set] <- 10 + rowMeans(u)

    } else if (stype == "gh_marginal") {
      # Stationary Gaussian field -> pointwise Tukey g-and-h transform with
      # spatially varying (g, h) -> standardise to mean 0 / variance 1 per site.
      set.seed(830000 + set)
      Z <- L.stat %*% matrix(rnorm(ns * nt), ns, nt)   # ns x nt, N(0,1) margins
      u <- sigma.add * (gh_transform(Z, g.field, h.field) - gh.mom$mean) / gh.mom$sd
      data$y <- data$y + u
      truemean.field[[setting]][, set] <- 10 + rowMeans(u)

    } else if (stype == "lambda_varying") {
      # nknots = 1 -> z_t is one scalar per day, so swapping the constant lambda
      # for a field is a pure post-hoc adjustment of the skew term.
      zt <- data$z[1, ]                                       # nt-vector
      data$y <- data$y + outer(lambda.field - lambda.par$lambda0, zt)
      truemean.field[[setting]][, set] <- 10 + lambda.field * mean(zt)

    } else {
      stop(sprintf("unhandled surf_type '%s'", stype))
    }

    y[, , set, setting]        <- data$y
    tau.t.setting[, , set]     <- data$tau
    z.t.setting[, , set]       <- data$z
    knots.t.setting[, , , set] <- data$knots
  }

  cat(sprintf(
    "finished setting %2d (%-14s | %-10s | %s)\n",
    setting, stype, moment.nonsta[setting], route.nonsta[setting]
  ))

  tau.t[[setting]]   <- tau.t.setting
  z.t[[setting]]     <- z.t.setting
  knots.t[[setting]] <- knots.t.setting
}

settings.nonsta <- data.frame(
  setting   = seq_len(nsettings),
  dist      = dist.nonsta,
  nknots    = nknots.nonsta,
  lambda    = lambda.nonsta,
  surf_type = surf.type,
  moment    = moment.nonsta,
  route     = route.nonsta,
  reference = ref.nonsta,
  stringsAsFactors = FALSE
)

save(
  y, tau.t, z.t, knots.t, settings.nonsta,
  ns, nt, s, nsets, ntest, x,
  # unified truth contract read by scores.R
  truemean.field,
  # ground-truth surface building blocks (settings 1-2)
  f.basis, c1, c2, surf.scale, surf.coef.invariant, W.varying, M.var,
  # setting 3 non-stationary-dependence random-effect building blocks
  F.re, W.nsdep, kappa.re, centers.re, tau.re,
  # settings 4-10 building blocks (also read by diagnose_nonsta.R)
  C.stat, sigma.add,
  Sigma.ps, C.ps, ps.par,
  rho.field, C.covdep, covdep.par,
  C.fuentes, fuentes.par,
  g.field, h.field, gh.par, gh.mom,
  lambda.field, lambda.par,
  file = "simdata_nonsta.RData"
)

cat(sprintf(
  "\nsaved simdata_nonsta.RData: %d settings, dim(y) = [%s]\n",
  nsettings, paste(dim(y), collapse = ", ")
))
