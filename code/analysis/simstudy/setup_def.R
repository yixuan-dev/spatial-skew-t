#########################################################################
# Simulation study with the geometric-anisotropy ("deformed") exponential
# covariance.  Mirrors setup.R but switches CorFx -> CorFxDef.
#
# Three settings, all skew-t-1 (lambda = 3, dist = "t", nknots = 1),
# varying only the deformation parameters (theta, ratio):
#   1 - theta = 0,    ratio = 1.00   (isotropic baseline)
#   2 - theta = pi/4, ratio = 0.50   (moderate anisotropy at 45 deg)
#   3 - theta = pi/6, ratio = 0.25   (strong anisotropy at 30 deg)
#########################################################################

rm(list = ls())

library(fields)
library(SpatialTools)

# Source prop/auxfunctions.R which now exposes:
#   CorFxDef()   - geometric-anisotropy exponential
#   rpotspatTS() - accepts cov.type = "deformed", theta, ratio
source("../../R/prop/auxfunctions.R")

# Pure-R override of mem() to avoid the Rcpp build of g.Rcpp.
# When nknots = 1 every site simply maps to knot 1, but keep the
# generic form for completeness.
mem <- function(s, knots) {
  d <- fields::rdist(s, knots)
  apply(d, 1, which.min)
}

# Common parameters (match setup.R)
beta.t      <- c(10, 0, 0)
nu.t        <- 0.5      # ignored when cov.type = "deformed"
gamma.t     <- 0.9
rho.t       <- 1
tau.alpha.t <- 3
tau.beta.t  <- 8

# Per-setting deformation parameters
dist.def   <- c("t",  "t",  "t")
nknots.def <- c(1,    1,    1)
lambda.def <- c(3,    3,    3)
theta.def  <- c(0,    pi / 4, pi / 6)
ratio.def  <- c(1.00, 0.50,  0.25)
nsettings  <- length(theta.def)

# Sites & covariates (same RNG seed and grid as setup.R)
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

# Storage
y       <- array(NA, dim = c(ns, nt, nsets, nsettings))
tau.t   <- vector("list", length = nsettings)
z.t     <- vector("list", length = nsettings)
knots.t <- vector("list", length = nsettings)

for (setting in seq_len(nsettings)) {
  nknots <- nknots.def[setting]
  tau.t.setting   <- array(NA, dim = c(nknots, nt, nsets))
  z.t.setting     <- array(NA, dim = c(nknots, nt, nsets))
  knots.t.setting <- array(NA, dim = c(nknots, nt, 2, nsets))

  for (set in seq_len(nsets)) {
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
      lambda    = lambda.def[setting],
      tau.alpha = tau.alpha.t,
      tau.beta  = tau.beta.t,
      nknots    = nknots,
      dist      = dist.def[setting],
      cov.type  = "deformed",
      theta     = theta.def[setting],
      ratio     = ratio.def[setting]
    )

    y[, , set, setting]        <- data$y
    tau.t.setting[, , set]     <- data$tau
    z.t.setting[, , set]       <- data$z
    knots.t.setting[, , , set] <- data$knots
  }

  cat(sprintf(
    "finished setting %d (theta = %.4f, ratio = %.2f)\n",
    setting, theta.def[setting], ratio.def[setting]
  ))

  tau.t[[setting]]   <- tau.t.setting
  z.t[[setting]]     <- z.t.setting
  knots.t[[setting]] <- knots.t.setting
}

settings.def <- data.frame(
  setting = seq_len(nsettings),
  dist    = dist.def,
  nknots  = nknots.def,
  lambda  = lambda.def,
  theta   = theta.def,
  ratio   = ratio.def
)

save(
  y, tau.t, z.t, knots.t, settings.def,
  ns, nt, s, nsets, ntest, x,
  file = "simdata_def.RData"
)

cat(sprintf(
  "saved simdata_def.RData with %d deformed-covariance settings\n",
  nsettings
))
