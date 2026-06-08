#########################################################################
# Simulation study with a NON-STATIONARY spatial mean surface added on top
# of the skew-t-1 data-generating process.  Mirrors setup_def.R in style,
# but instead of deforming the covariance it injects a smooth, location-
# based mean surface g(s) built from the cosine "bump" basis of
# Tzeng & Huang (2018, Scenario 1).
#
# Purpose: give the MRTS covariate extension a data-generating process that
# actually CONTAINS spatial mean structure for the basis to recover, so the
# choice of K_MRTS can be justified by a recovery curve (Brier/quantile
# score vs K) rather than asserted.
#
# Two settings, both skew-t-1 (lambda = 3, dist = "t", nknots = 1):
#   1 - TIME-INVARIANT surface:  g(s) = a1 f1(s) + a2 f2(s), fixed (a1, a2).
#       The truth lies inside the space the Morris model can represent
#       (beta is pooled over time -> a single, time-constant coefficient
#       per covariate; see updateBeta() in R/ar2/update_params.R).  MRTS
#       can therefore recover g(s) exactly in the limit -> the score-vs-K
#       curve has a clean elbow and K can be justified as "past the plateau".
#
#   2 - TIME-VARYING surface:    g_t(s) = w1(t) f1(s) + w2(t) f2(s),
#       (w1(t), w2(t)) ~ N(0, diag(25, 9)) per Tzeng & Huang.  Because the
#       model's beta is time-constant, this structure is NOT fully in the
#       representable space; MRTS can at best capture the time-average.
#       This is the robustness / misspecification companion to setting 1.
#
# The two basis fields f1, f2 are evaluated once at all sites (s is fixed),
# and the building blocks (f.basis, surface coefficients / weights) are
# saved so the recovery analysis can reconstruct the true surface exactly.
#########################################################################

rm(list = ls())

library(fields)
library(SpatialTools)

# Exposes rpotspatTS(), CorFx(), mem(), etc.
source("../../R/ar2/auxfunctions.R")

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
dist.nonsta   <- c("t", "t")
nknots.nonsta <- c(1,   1)
lambda.nonsta <- c(3,   3)
surf.type     <- c("invariant", "varying")
nsettings     <- length(surf.type)

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
# Storage
# -----------------------------------------------------------------------
y       <- array(NA, dim = c(ns, nt, nsets, nsettings))
tau.t   <- vector("list", length = nsettings)
z.t     <- vector("list", length = nsettings)
knots.t <- vector("list", length = nsettings)

# Save the surface building blocks so the true mean can be reconstructed:
#   setting 1: g(s)      = f.basis %*% a.fixed              (same for all sets)
#   setting 2: g_t(s)    = f.basis %*% t(W[[set]])          (nt x 2 weights)
surf.coef.invariant <- a.fixed
W.varying           <- vector("list", length = nsets)  # filled below

L.var <- chol(M.var)  # M.var = L'L

for (setting in seq_len(nsettings)) {
  nknots <- nknots.nonsta[setting]
  tau.t.setting   <- array(NA, dim = c(nknots, nt, nsets))
  z.t.setting     <- array(NA, dim = c(nknots, nt, nsets))
  knots.t.setting <- array(NA, dim = c(nknots, nt, 2, nsets))

  for (set in seq_len(nsets)) {
    # Reproducible skew-t draw (mirror setup_def.R seeding scheme)
    set.seed(setting * 1000 + set)
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
      dist      = dist.nonsta[setting]
    )

    # Add the non-stationary mean surface on top of the skew-t draw.
    # Because the surface enters the mean additively and the skew-t error
    # is independent of the mean, adding it post-hoc here is identical to
    # injecting it into `mu` inside rpotspatTS().
    if (surf.type[setting] == "invariant") {
      g_s <- as.vector(f.basis %*% a.fixed)        # ns-vector, constant in t
      data$y <- data$y + matrix(g_s, ns, nt)
    } else {
      # Tzeng time-varying weights, drawn reproducibly per dataset.
      # Generated only once (during setting 2) and stored for reuse.
      set.seed(900000 + set)
      W <- matrix(rnorm(nt * 2), nt, 2) %*% L.var  # nt x 2
      W.varying[[set]] <- W
      data$y <- data$y + f.basis %*% t(W)          # ns x nt
    }

    y[, , set, setting]        <- data$y
    tau.t.setting[, , set]     <- data$tau
    z.t.setting[, , set]       <- data$z
    knots.t.setting[, , , set] <- data$knots
  }

  cat(sprintf(
    "finished setting %d (%s surface)\n",
    setting, surf.type[setting]
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
  stringsAsFactors = FALSE
)

save(
  y, tau.t, z.t, knots.t, settings.nonsta,
  ns, nt, s, nsets, ntest, x,
  # ground-truth surface building blocks (for the recovery analysis)
  f.basis, c1, c2, surf.scale, surf.coef.invariant, W.varying, M.var,
  file = "simdata_nonsta.RData"
)

cat(sprintf(
  "saved simdata_nonsta.RData with %d non-stationary-mean settings\n",
  nsettings
))
