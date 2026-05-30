#########################################################################
# A simulation study to determine if a thresholded or skew methods improve
# return-level estimation and threshold exceedance prediction over
# standard kriging methods
#
# data settings:
#   1 - Gaussian
#   2 - t-1
#   3 - t-5
#   4 - skew t-1 (lambda = 3)
#   5 - skew t-5 w/partition (lambda = 3)
#   6 - max-stable with mu=1, sig=1, xi=0.2, alpha = 0.5, bw = 1
#   7 - x = setting 4, set T = q(0.80)
#       y = x,              x > T
#       y = T * exp(x - T), x <= T
#   8 - Brown-Resnick with range = 1, smooth = 0.5
#   9 - skew t-1 with fixed AR(2), stronger time dependence: phi1 = 0.80, phi2 = -0.35
#  10 - skew t-1 with fixed AR(2), weaker time dependence: phi1 = 0.12, phi2 = -0.05
#  11 - skew t-5 with fixed AR(2), stronger time dependence: phi1 = 0.80, phi2 = -0.35
#  12 - skew t-5 with fixed AR(2), weaker time dependence:   phi1 = 0.12, phi2 = -0.05
#  13 - skew t-1, AR(2) on z only (phi_z = (0.80, -0.35), phi_tau = phi_w = 0): pure level channel
#  14 - skew t-1, AR(2) on tau only (phi_tau = (0.80, -0.35), phi_z = phi_w = 0): pure volatility channel
#  15 - skew t-5, AR(2) on w only  (phi_w  = (0.80, -0.35), phi_z = phi_tau = 0): pure knot-mixing channel (K=5 needed)
#
# analysis methods:
#  1 - Gaussian
#  2 - skew t-1
#  3 - t-1 (T = 0.80)
#  4 - skew t-5
#  5 - t-5 (T = 0.80)
#  6 - Max-stable (T = 0.80)
#  7 - Skew-t, K = 1 + AR(2) temporal
#  8 - Skew-t, K = 5 + AR(2) temporal
#
#########################################################################
##### For revisions, see below initial generation
# clean out previous variables
rm(list = ls())

# libraries
library(fields)
library(SpatialTools)

# necessary functions
# source("../../R/mcmc_cont_lambda.R", chdir = T)
# source("../../R/auxfunctions.R")
# source("../../biometrics/auxfunctions.R")
source("../../R/ar2/update_params_cpp.R")
source("../../R/ar2/auxfunctions.R")
source("../../R/ar2/auxfunctions_ar2.R")
source("./max-stab/Bayes_GEV.R")

# data settings
beta.t <- c(10, 0, 0)
nu.t <- 0.5
gamma.t <- 0.9
dist.t <- c("gaussian", "t", "t", "t", "t")
nknots.t <- c(1, 1, 5, 1, 5, NA, NA, NA, 1, 1, 5, 5, 1, 1, 5)
rho.t <- c(1, 1, 1, 1, 1)
lambda.t <- c(0, 0, 0, 3, 3)
tau.alpha.t <- 6
tau.beta.t <- 16

# covariate data
set.seed(20)
s <- cbind(runif(144, 0, 10), runif(144, 0, 10))
knots.x <- seq(1, 9, length = 12)
knots.gev <- expand.grid(knots.x, knots.x)
ns <- nrow(s)
nt <- 50
nsets <- 50
nsettings <- 15
ntest <- 44

x <- array(1, c(ns, nt, 3))
for (t in 1:nt) {
  x[, t, 2] <- s[, 1]
  x[, t, 3] <- s[, 2]
}

build_phi_triple <- function(setting) {
  phi_strong <- c(0.80, -0.35)
  phi_weak   <- c(0.12, -0.05)
  phi_zero   <- c(0, 0)

  triple <- switch(as.character(setting),
    "9"  = list(phi.z = phi_strong, phi.tau = phi_strong, phi.w = phi_strong),
    "10" = list(phi.z = phi_weak,   phi.tau = phi_weak,   phi.w = phi_weak),
    "11" = list(phi.z = phi_strong, phi.tau = phi_strong, phi.w = phi_strong),
    "12" = list(phi.z = phi_weak,   phi.tau = phi_weak,   phi.w = phi_weak),
    "13" = list(phi.z = phi_strong, phi.tau = phi_zero,   phi.w = phi_zero),
    "14" = list(phi.z = phi_zero,   phi.tau = phi_strong, phi.w = phi_zero),
    "15" = list(phi.z = phi_zero,   phi.tau = phi_zero,   phi.w = phi_strong),
    stop(sprintf("No phi triple configured for setting %s", setting))
  )

  for (nm in names(triple)) {
    p <- triple[[nm]]
    if (any(p != 0) && !check_ar2_stability(p[1], p[2])) {
      stop(sprintf(
        "Non-stationary AR(2) in %s for setting %s: phi1=%.3f, phi2=%.3f",
        nm, setting, p[1], p[2]
      ))
    }
  }

  triple
}

# Storage for datasets
y <- array(NA, dim = c(ns, nt, nsets, nsettings))
tau.t <- z.t <- knots.t <- vector("list", length = nsettings)
phi.path <- vector("list", length = nsettings)

for (setting in 1:nsettings) {
  if (setting < 6) {
    nknots <- nknots.t[setting]
    tau.t.setting <- array(NA, dim = c(nknots, nt, nsets))
    z.t.setting <- array(NA, dim = c(nknots, nt, nsets))
    knots.t.setting <- array(NA, dim = c(nknots, nt, 2, nsets))
    for (set in 1:nsets) {
      set.seed(setting * 100 + set)
      data <- rpotspatTS(
        nt = nt, x = x, s = s, beta = beta.t, gamma = gamma.t, nu = nu.t,
        rho = rho.t[setting], phi.z = 0, phi.w = 0, phi.tau = 0,
        lambda = lambda.t[setting], tau.alpha = tau.alpha.t,
        tau.beta = tau.beta.t, nknots = nknots.t[setting],
        dist = dist.t[setting]
      )

      y[, , set, setting] <- data$y
      tau.t.setting[, , set] <- data$tau
      z.t.setting[, , set] <- data$z
      knots.t.setting[, , , set] <- data$knots
    }
    print(setting)
    tau.t[[setting]] <- tau.t.setting
    z.t[[setting]] <- z.t.setting
    knots.t[[setting]] <- knots.t.setting
  } else if (setting == 6) {
    for (set in 1:nsets) {
      set.seed(setting * 100 + set)
      y[, , set, setting] <- rgevspatial(nreps = nt, S = s, knots = knots.gev, xi = 0.2)
    }
    print(setting)
  } else if (setting == 7) {
    tau.t.setting <- array(NA, dim = c(nknots.t[4], nt, nsets))
    z.t.setting <- array(NA, dim = c(nknots.t[4], nt, nsets))
    knots.t.setting <- array(NA, dim = c(nknots.t[4], nt, 2, nsets))
    for (set in 1:nsets) {
      set.seed(setting * 100 + set)
      data <- rpotspatTS(
        nt = nt, x = x, s = s, beta = beta.t, gamma = gamma.t, nu = nu.t,
        rho = rho.t[4], phi.z = 0, phi.w = 0, phi.tau = 0,
        lambda = lambda.t[4], tau.alpha = tau.alpha.t,
        tau.beta = tau.beta.t, nknots = nknots.t[4],
        dist = dist.t[4]
      )
      y.set <- data$y
      y.quant <- quantile(y.set, probs = 0.80)
      y.set <- ifelse(
        y.set < y.quant,
        y.quant * exp((y.set - y.quant) / 3),
        y.set
      )
      y[, , set, setting] <- y.set

      tau.t.setting[, , set] <- data$tau
      z.t.setting[, , set] <- data$z
      knots.t.setting[, , , set] <- data$knots
    }
    tau.t[[setting]] <- tau.t.setting
    z.t[[setting]] <- z.t.setting
    knots.t[[setting]] <- knots.t.setting
    print(setting)
  } else if (setting >= 9) {
    phi.triple <- build_phi_triple(setting)
    phi.path[[setting]] <- list(
      type = "channel-specific",
      phi.z   = phi.triple$phi.z,
      phi.tau = phi.triple$phi.tau,
      phi.w   = phi.triple$phi.w
    )
    nknots <- nknots.t[setting]
    tau.t.setting <- array(NA, dim = c(nknots, nt, nsets))
    z.t.setting <- array(NA, dim = c(nknots, nt, nsets))
    knots.t.setting <- array(NA, dim = c(nknots, nt, 2, nsets))
    for (set in 1:nsets) {
      set.seed(setting * 100 + set)
      data <- rpotspatTS_arp(
        nt = nt, x = x, s = s, beta = beta.t, gamma = gamma.t, nu = nu.t,
        rho = rho.t[4], lambda = lambda.t[4], tau.alpha = tau.alpha.t,
        tau.beta = tau.beta.t, nknots = nknots, dist = dist.t[4],
        phi.z = phi.triple$phi.z, phi.w = phi.triple$phi.w, phi.tau = phi.triple$phi.tau
      )
      y[, , set, setting] <- data$y
      tau.t.setting[, , set] <- data$tau
      z.t.setting[, , set] <- data$z
      knots.t.setting[, , , set] <- data$knots
    }
    print(setting)
    tau.t[[setting]] <- tau.t.setting
    z.t[[setting]] <- z.t.setting
    knots.t[[setting]] <- knots.t.setting
  }
}

save(y, tau.t, z.t, knots.t, phi.path, ns, nt, s, nsets, ntest,
  x, # covariate data that should be the same for all datasets
  file = "simdata.RData"
)

# modify simdata for revisions
rm(list = ls())
library(fields)
library(SpatialTools)
library(SpatialExtremes)
source("../../R/mcmc_cont_lambda.R", chdir = T)
source("../../R/auxfunctions.R")
source("./max-stab/Bayes_GEV.R")
load(file = "simdata.RData")
ns <- dim(y)[1]
nt <- dim(y)[2]
nsets <- dim(y)[3]
nsettings <- dim(y)[4]
y.new <- y
setting <- 8
for (set in 1:nsets) {
  set.seed(setting * 100 + set)
  y.new[, , set, 8] <- t(rmaxstab(
    n = nt, coord = s, cov.mod = "brown", range = 1,
    smooth = 0.5
  ))
  print(paste("finished set", set))
}
y <- y.new

save(y, tau.t, z.t, knots.t, phi.path, ns, nt, s, nsets, ntest,
  x, # covariate data that should be the same for all datasets
  file = "simdata.RData"
)

# par(mfrow=c(2, 3))
# quilt.plot(s[, 1], s[, 2], z=y[, 1, 1, 1], nx=20, ny=20)
# quilt.plot(s[, 1], s[, 2], z=y[, 1, 1, 2], nx=20, ny=20)
# quilt.plot(s[, 1], s[, 2], z=y[, 1, 1, 3], nx=20, ny=20)
# quilt.plot(s[, 1], s[, 2], z=y[, 1, 1, 4], nx=20, ny=20)
# quilt.plot(s[, 1], s[, 2], z=y[, 1, 1, 5], nx=20, ny=20)
# quilt.plot(s[, 1], s[, 2], z=y[, 1, 1, 6], nx=20, ny=20)

# # plot the data
# s1 <- s2 <- seq(0, 1, length = 50)
# s <- expand.grid(s1, s2)
# X <- array(1, dim = c(nrow(s), 2, 3))
# X.temp <- cbind(rep(1, nrow(s)), s[, 1], s[, 2])
# X[, 1, ] <- X.temp
# X[, 2, ] <- X.temp
# y.gau <- rpotspat(
#   nt = 2, x = X, s = s, beta = c(10, 0, 0), alpha = alpha.t, nu = nu.t,
#   gau.rho = 0.1, t.rho = 0.1, mixprob = 0, z.alpha = 0, tau.alpha = tau.alpha.t,
#   tau.beta = tau.beta.t, nknots = 1
# )
# par(mfrow = c(1, 2))
# quilt.plot(s[, 1], s[, 2], z = y.gau$y[, 1], nx = 50, ny = 50, main = "Gaussian")
# hist(y.gau$y, main = "Histogram", xlab = "")

# y.t1 <- rpotspat(
#   nt = 2, x = X, s = s, beta = c(10, 0, 0), alpha = alpha.t, nu = nu.t,
#   gau.rho = 0.1, t.rho = 0.1, mixprob = 1, z.alpha = 0, tau.alpha = tau.alpha.t,
#   tau.beta = tau.beta.t, nknots = 1
# )
# par(mfrow = c(1, 2))
# quilt.plot(s[, 1], s[, 2], z = y.t1$y[, 1], nx = 50, ny = 50, main = "t, K=1")
# hist(y.t1$y, main = "Histogram", xlab = "")

# y.t3 <- rpotspat(
#   nt = 2, x = X, s = s, beta = c(10, 0, 0), alpha = alpha.t, nu = nu.t,
#   gau.rho = 0.1, t.rho = 0.1, mixprob = 1, z.alpha = 0, tau.alpha = tau.alpha.t,
#   tau.beta = tau.beta.t, nknots = 3
# )
# par(mfrow = c(1, 2))
# quilt.plot(s[, 1], s[, 2], z = y.t3$y[, 1], nx = 50, ny = 50, main = "t, K=3")
# hist(y.t3$y, main = "Histogram", xlab = "")

# y.st1 <- rpotspat(
#   nt = 2, x = X, s = s, beta = c(10, 0, 0), alpha = alpha.t, nu = nu.t,
#   gau.rho = 0.1, t.rho = 0.1, mixprob = 1, z.alpha = 3, tau.alpha = tau.alpha.t,
#   tau.beta = tau.beta.t, nknots = 1
# )
# par(mfrow = c(1, 2))
# quilt.plot(s[, 1], s[, 2], z = y.st1$y[, 1], nx = 50, ny = 50, main = "skew-t, K=1, alpha=3")
# hist(y.st1$y[, 1], main = "Histogram", xlab = "")

# y.st3 <- rpotspat(
#   nt = 2, x = X, s = s, beta = c(10, 0, 0), alpha = alpha.t, nu = nu.t,
#   gau.rho = 0.1, t.rho = 0.1, mixprob = 1, z.alpha = 1, tau.alpha = tau.alpha.t,
#   tau.beta = tau.beta.t, nknots = 3
# )
# par(mfrow = c(1, 2))
# quilt.plot(s[, 1], s[, 2], z = y.s35$y[, 1], nx = 50, ny = 50, main = "skew-t, K=3, alpha=1")
# hist(y.st5$y[, 1], main = "Histogram", xlab = "")

# par(mfrow = c(1, 2))
# quilt.plot(s[, 1], s[, 2], z = y.gau$y[, 1], nx = 50, ny = 50, main = "Gaussian")
# quilt.plot(s[, 1], s[, 2], z = y.st3$y[, 1], nx = 50, ny = 50, main = "skew-t, K=3, alpha=1")

# par(mfrow = c(1, 2))
# quilt.plot(s[, 1], s[, 2], z = y.gau$y[, 1], nx = 50, ny = 50, main = "Gaussian")
# quilt.plot(s[, 1], s[, 2], z = y.t3$y[, 1], nx = 50, ny = 50, main = "t, K=3")

# # remove simstudy "truth" settings to help diagnose errors.
# rm(list = ls())
