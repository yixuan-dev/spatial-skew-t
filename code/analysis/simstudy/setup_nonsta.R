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
# ---------------------------------------------------------------------
# NON-t EXTENSION (settings 13-21).  Settings 1-12 all share the skew-t
# core, so the fitted skew-t is CORRECTLY specified at large K -- a
# confound for every "MRTS helps/hurts" conclusion.  Settings 13-21 drop
# rpotspatTS entirely: nonlinear mean surfaces (orthogonalised against the
# K=0 design [1, s1, s2], standardised to sd 11) plus Gaussian-GP errors
# (sinh-arcsinh transformed for setting 18).  sigma.g = 4 was calibrated by
# non_stationary/nont_precheck.R to put the oracle Brier ceiling R* in the
# 0.4-0.5 band (median 0.446), comparable to settings 11/12.
#
#  ID  surf_type        mean surface                       error
#  --  ---------------  ---------------------------------  --------------
#  13  franke           Franke's function (multiscale)     Gaussian GP
#  14  ridge_curved     curved (banana) Gaussian ridge     Gaussian GP
#  15  annulus          Ricker ring, r0 = 3, w = 1.5       Gaussian GP
#  16  sigmoid_front    tanh arc front                     Gaussian GP
#  17  gp_mean          GP draw per dataset, fixed in t    Gaussian GP
#  18  franke_sas       same surface as 13                 sinh-arcsinh GP
#  19  poly2            bowl + saddle (2nd-order poly)     Gaussian GP
#  20  transform_mix    log + exp-decay + ratio terms      Gaussian GP
#  21  nonsmooth        hinge + |.| (C^0 creases)          Gaussian GP
#
# simdata_nonsta.RData carries all 21 settings (the 13-21 extension was
# generated as simdata_nonsta_ext.RData, verified against the frozen
# original by non_stationary/nont_verify.R -- settings 1-12 bit-identical
# -- and then merged back on 2026-07-17).  Settings 1-12 regenerate
# bit-for-bit, so existing fits in results_nonsta/ keyed by setting id
# remain valid.
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
# Exposes the settings 13-21 surface builders, ortho_std(), sas_transform(),
# sas_moments(), matern_cor()
source("non_stationary/nont_builders.R")

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

# Settings 1-12 are skew-t-1; settings 13-21 are the non-t extension
surf.type <- c("invariant", "varying", "ns_dependence",
               "ps_cov", "ps_add", "covdep_cov", "covdep_add",
               "fuentes_cov", "gh_marginal", "lambda_varying",
               "invariant_strong", "mrts_multibump",
               "franke", "ridge_curved", "annulus", "sigmoid_front",
               "gp_mean", "franke_sas", "poly2", "transform_mix",
               "nonsmooth")
nsettings <- length(surf.type)
nont.first <- 13L                       # settings >= this skip rpotspatTS

dist.nonsta   <- c(rep("t", 12),
                   rep("gauss", 5), "gauss_sas", rep("gauss", 3))
nknots.nonsta <- rep(1L,  nsettings)    # >= 13: dims of the NA tau/z/knots arrays
lambda.nonsta <- c(rep(3, 12),          # setting 10 varies lambda(s) around this
                   rep(NA_real_, 9))    # no skew term in the non-t settings

# Which moment each setting perturbs, and how it is injected.
moment.nonsta <- c("mean", "mean", "dependence",
                   "dependence", "dependence", "dependence", "dependence",
                   "dependence", "tails", "skewness",
                   "mean", "mean",
                   rep("mean", 9))
route.nonsta  <- c("mean", "mean", "additive",
                   "replace_C", "additive", "replace_C", "additive",
                   "replace_C", "additive", "skew_term",
                   "mean", "mean",
                   rep("mean", 9))
ref.nonsta    <- c("Tzeng & Huang (2018)", "Tzeng & Huang (2018)",
                   "Cressie & Johannesson (2008)",
                   "Paciorek & Schervish (2006)", "Paciorek & Schervish (2006)",
                   "Schmidt et al. (2011); Risser & Calder (2015)",
                   "Schmidt et al. (2011); Risser & Calder (2015)",
                   "Fuentes (2001, 2002)", "Xu & Genton (2017)",
                   "Morris et al. (2017) assumption violated",
                   "Tzeng & Huang (2018), high-SNR",
                   "multi-bump mean showcase (data-limited resolution)",
                   "Franke (1979)", "curved ridge (direction-varying)",
                   "Ricker annulus (non-monotone radial)",
                   "tanh arc front (steepest recoverable)",
                   "random smooth mean (GP function class)",
                   "Jones & Pewsey (2009) sinh-arcsinh error",
                   "2nd-order polynomial (regression terms)",
                   "log/exp/ratio transforms (regression terms)",
                   "hinge + |.| creases (C^0 vs C^inf basis)")

# A/B pairs share a base seed so the underlying draw is identical:
# setting 5 borrows setting 4's, setting 7 borrows setting 6's, and setting 18
# borrows setting 13's (same Franke mean, same Gaussian innovations -- only the
# sinh-arcsinh transform differs).  Settings 1-3 map to themselves, which
# reproduces the original seeds exactly.
pair.seed <- c(1, 2, 3, 4, 4, 6, 6, 8, 9, 10, 11, 12,
               13, 14, 15, 16, 17, 13, 19, 20, 21)

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
# Setting 11: the SAME surface as setting 1, scaled up 3x.  This is the BRIER
# POSITIVE CONTROL, and it exists because setting 1 turns out to be incapable of
# being one.
#
# MRTS only touches the fixed-effect mean and beta is pooled over time, so the
# most any K can do is recover the TIME-AVERAGED mean surface.  The best relative
# Brier a mean basis could ever achieve is therefore
#
#   R* = BS(true mean surface) / BS(spatial constant),
#
# computable with no MCMC at all (see the oracle calculation in the design doc).
# At surf.scale = 1 -- i.e. setting 1 -- that ceiling is R* = 0.86, so even a
# PERFECT mean basis buys at most 14% of Brier; we measured 1.5%.  Setting 1 was
# never able to demonstrate that MRTS helps the Brier score, which left every
# other setting's Brier null uninterpretable ("maybe MRTS never helps").
#
# surf.scale.strong = 3 puts sd(mu) = 11.0 against a predictive sd of 7.5
# (SNR 1.7) and drops the ceiling to R* = 0.375.  That is an unambiguous
# positive control.
# -----------------------------------------------------------------------
surf.scale.strong <- 3.0
a.fixed.strong    <- c(5, 3) * surf.scale.strong   # = c(15, 9); sd(g) = 11.04

# -----------------------------------------------------------------------
# Setting 12: a RICH mean surface that showcases where MRTS earns its keep.
# Settings 1/11 are only two cosine bumps -- so low-rank that MRTS's projection
# recovery is 90% done by K=5, showing no reason to want a large basis. Setting
# 12 is SIX medium Gaussian bumps at irregular locations with alternating signs
# (width 0.16 = 1.6 units), standardised to sd = 11 (high SNR, matching 11).
#
# Its noise-free MRTS projection recovery (mrts_recovery_curve.R) declines
# gradually -- 0.71 at K=5, 0.19 at K=15, 0.06 at K=50 -- versus the 2-bump
# surface's 0.21 at K=5. THAT extended elbow is the point: a rich mean needs a
# large K, and MRTS delivers it.
#
# A finer (high-frequency) component was tested and dropped: on this 100-train/
# 44-test geometry a 3-cycle ripple is UN-recoverable (MRTS overfits the train
# sites and extrapolates wildly at held-out sites, projection error rising with
# K). MRTS's usable resolution is capped by site density -- a genuine finding,
# recorded in nonsta_settings_design.md. Six medium bumps sit just inside that
# recoverable ceiling.
#
# Centres are the realisation of set.seed(7); runif(6, .15, .85), hardcoded so
# the surface is RNG-independent. Time-invariant -> pooled beta recovers it.
# The surface itself (g.multibump) is built below, once s01 exists.
# -----------------------------------------------------------------------
ms.par <- list(
  centers = rbind(
    c(0.8422, 0.3880), c(0.4284, 0.8304), c(0.2310, 0.2661),
    c(0.1988, 0.4714), c(0.3206, 0.2702), c(0.7044, 0.3120)),
  signs = c(1, -1, 1, -1, 1, -1),
  width = 0.16,
  target_sd = 11
)

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

# Setting 12 surface (see ms.par above): 6 alternating-sign Gaussian bumps,
# standardised to sd = 11.  Time-invariant, so pooled beta can recover it.
g.multibump.raw <- rowSums(sapply(seq_len(nrow(ms.par$centers)), function(k) {
  d2 <- (s01[, 1] - ms.par$centers[k, 1])^2 + (s01[, 2] - ms.par$centers[k, 2])^2
  ms.par$signs[k] * exp(-d2 / (2 * ms.par$width^2))
}))
g.multibump <- ms.par$target_sd *
  (g.multibump.raw - mean(g.multibump.raw)) / sd(g.multibump.raw)  # ns-vector, sd 11

# -----------------------------------------------------------------------
# Settings 13-21: non-t DGPs.  Every mean surface is orthogonalised against
# the K=0 design [1, s1, s2] (so sd 11 measures only what the baseline
# cannot fit) and standardised to sd 11.  Errors are sigma.g * (C.stat
# Gaussian GP), i.i.d. in time; setting 18 pushes the same innovations
# through the sinh-arcsinh transform (standardised by MC moments, like the
# g-and-h of setting 9).
#
# Constants fixed by non_stationary/nont_precheck.R:
#   sigma.g = 4          -> oracle Brier ceiling R* median 0.446 (band .4-.5)
#   annulus width = 1.5  -> projection error declines monotonically with K
#                           (1.2 was too fine: error ROSE around K = 15)
#   franke dip unwidened -> recoverable (0.07 at K = 50)
# -----------------------------------------------------------------------
sigma.g <- 4

g.nont <- list(
  franke        = ortho_std(franke_fn(s01),             s),
  ridge_curved  = ortho_std(ridge_fn(s01),              s),
  annulus       = ortho_std(annulus_fn(s, width = 1.5), s),
  sigmoid_front = ortho_std(front_fn(s),                s),
  poly2         = ortho_std(poly2_fn(s01),              s),
  transform_mix = ortho_std(transform_mix_fn(s),        s),
  nonsmooth     = ortho_std(nonsmooth_fn(s),            s)
)
g.nont$franke_sas <- g.nont$franke   # setting 18 = setting 13's surface

# Setting 17: one random smooth mean per dataset (fixed across t, so the
# pooled-beta mean channel can recover it).  Matern nu = 1.5, range 2 --
# smoother and longer-ranged than the exponential error.
gpmean.par <- list(range = 2, smoothness = 1.5, seed.base = 860000)
L.gpmean   <- t(chol(matern_cor(s, range = gpmean.par$range,
                                smoothness = gpmean.par$smoothness)))
G.gpmean <- sapply(seq_len(nsets), function(set) {
  set.seed(gpmean.par$seed.base + set)
  ortho_std(as.vector(L.gpmean %*% rnorm(ns)), s)
})                                       # ns x nsets, one surface per column

# Setting 18: sinh-arcsinh error, X = sinh(delta * asinh(Z) + eps).
# eps = 0.8 gives a right skew comparable to the skew-t core's; delta = 1
# keeps the tails light (all moments finite) -- skewed but decisively non-t.
sas.par <- list(eps = 0.8, delta = 1)
sas.mom <- sas_moments(sas.par$eps, sas.par$delta)

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
    # ---- Non-t settings (13-21): no rpotspatTS ----------------------------
    # y = 10 + g(s) + sigma.g * (C.stat Gaussian GP), i.i.d. in time.
    # Setting 18 shares setting 13's base seed, so Z is IDENTICAL and the
    # sinh-arcsinh transform is the only difference (A/B pair).  tau/z/knots
    # stay NA for these settings; scores.R reads only truemean.field.
    if (setting >= nont.first) {
      set.seed(pair.seed[setting] * 1000 + set)
      Z <- L.stat %*% matrix(rnorm(ns * nt), ns, nt)
      err <- if (stype == "franke_sas") {
        sigma.g * (sas_transform(Z, sas.par$eps, sas.par$delta) - sas.mom$mean) / sas.mom$sd
      } else {
        sigma.g * Z
      }
      g_s <- if (stype == "gp_mean") G.gpmean[, set] else g.nont[[stype]]
      y[, , set, setting] <- 10 + g_s + err   # g_s recycles down columns
      truemean.field[[setting]][, set] <- 10 + g_s
      next
    }

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

    } else if (stype == "invariant_strong") {
      # Same surface as setting 1, 3x the amplitude: the Brier positive control.
      g_s <- as.vector(f.basis %*% a.fixed.strong)
      data$y <- data$y + matrix(g_s, ns, nt)
      truemean.field[[setting]][, set] <- 10 + g_s

    } else if (stype == "mrts_multibump") {
      # Six medium Gaussian bumps, sd 11, time-invariant: the MRTS showcase --
      # rich enough that recovery needs a large K (see mrts_recovery_curve.R).
      data$y <- data$y + matrix(g.multibump, ns, nt)
      truemean.field[[setting]][, set] <- 10 + g.multibump

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
  # ground-truth surface building blocks (settings 1-2, and 11)
  f.basis, c1, c2, surf.scale, surf.coef.invariant, W.varying, M.var,
  surf.scale.strong, a.fixed.strong,
  # setting 12 multi-bump mean showcase
  ms.par, g.multibump,
  # setting 3 non-stationary-dependence random-effect building blocks
  F.re, W.nsdep, kappa.re, centers.re, tau.re,
  # settings 4-10 building blocks (also read by diagnose_nonsta.R)
  C.stat, sigma.add,
  Sigma.ps, C.ps, ps.par,
  rho.field, C.covdep, covdep.par,
  C.fuentes, fuentes.par,
  g.field, h.field, gh.par, gh.mom,
  lambda.field, lambda.par,
  # settings 13-21 non-t building blocks
  sigma.g, g.nont, G.gpmean, gpmean.par, sas.par, sas.mom,
  file = "simdata_nonsta.RData"
)

cat(sprintf(
  "\nsaved simdata_nonsta.RData: %d settings, dim(y) = [%s]\n",
  nsettings, paste(dim(y), collapse = ", ")
))
