#########################################################################
# setup.R - data generation for the time-block forecast simulation study.
#
# Companion to tex/time_block_strategy/ar2_rethink.tex. The note diagnoses
# why a marginal Brier score under a spatial-only holdout is structurally
# blind to AR(2) temporal pooling (Propositions 1-2) and prescribes a
# time-block holdout that forecasts forward from the seam state
# (Section 4). This script generates the *data* for that redesign.
#
# Key difference from code/analysis/simstudy/setup.R:
#   - The Morris baseline varies the *observation family* across settings
#     and holds out spatial sites (ntest).
#   - Here every setting is the same skew-t, K = 1 family; settings vary
#     only the *strength of AR(2) temporal dependence* on the latent
#     processes (tau, z, w). The holdout is over *time blocks*, so the
#     record is long (nt = 200) and no sites are withheld.
#
# data settings (AR(2) coefficients phi = (phi1, phi2) shared by the
#                latent tau*, z*, w* Gaussian processes; eq. (2) of the note):
#   1 - i.i.d. in time      phi = (0, 0)        null control
#   2 - weak AR(2)          phi = (0.12, -0.05) small spectral radius
#   3 - moderate AR(2)      phi = (0.60, -0.30) the note's reference case
#   4 - strong AR(2)        phi = (0.80, -0.35) long memory horizon
#
# All settings: skew-t (dist = "t"), K = 1 knot, lambda = 3.
# With K = 1 the Voronoi membership is constant, so w carries no signal
# and the temporal effect lives entirely in tau and z -- this isolates
# the AR(2) signal exactly as Remark 8 of the note recommends.
#
# time-block holdout geometry (stored for run-settings.R / scores.R):
#   block_H      - forecast horizon per block (lead times 1..H)
#   block_seams  - last observed time T_o of each block; block b fits on
#                  the contiguous prefix y[, 1:T_o] and forecasts the
#                  window (T_o, T_o + H]  (Definition 7 of the note).
#########################################################################

# clean out previous variables
rm(list = ls())

# run from the script's own directory so the relative source()/save()
# paths below resolve regardless of the caller's working directory.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1]),
    winslash = "/", mustWork = FALSE
  ))
  if (dir.exists(script_dir)) setwd(script_dir)
}

# libraries
library(fields)
library(SpatialTools)
library(Rcpp)

# necessary functions (one directory deeper than code/analysis/simstudy,
# hence ../../../ rather than ../../)
source("../../../R/ar2/update_params_cpp.R")
source("../../../R/ar2/auxfunctions.R")
source("../../../R/ar2/auxfunctions_ar2.R")

# ---- data settings ----------------------------------------------------
beta.t <- c(10, 0, 0)
nu.t <- 0.5
gamma.t <- 0.9
rho.t <- 1
lambda.t <- 3
tau.alpha.t <- 6
tau.beta.t <- 16
nknots.t <- 1
dist.t <- "t"

# AR(2) coefficient pair per setting (shared by tau*, z*, w*).
phi.settings <- list(
  c(0.00, 0.00),   # 1 - i.i.d. control
  c(0.12, -0.05),  # 2 - weak
  c(0.60, -0.30),  # 3 - moderate (reference case in the note)
  c(0.80, -0.35)   # 4 - strong
)
setting.label <- c("iid", "weak", "moderate", "strong")
nsettings <- length(phi.settings)

# fail loudly if any non-control pair leaves the stationarity triangle
# (Definition 3 of the note); a non-stationary recursion would explode.
for (setting in 2:nsettings) {
  pair <- phi.settings[[setting]]
  if (!check_ar2_stability(pair[1], pair[2])) {
    stop(sprintf(
      "Setting %d AR(2) pair phi=(%.3f, %.3f) is non-stationary.",
      setting, pair[1], pair[2]
    ))
  }
}

# ---- covariate data (identical spatial layout to the Morris study) ----
set.seed(20)
s <- cbind(runif(144, 0, 10), runif(144, 0, 10))
ns <- nrow(s)
nt <- 200
nsets <- 50

x <- array(1, c(ns, nt, 3))
for (t in 1:nt) {
  x[, t, 2] <- s[, 1]
  x[, t, 3] <- s[, 2]
}

# ---- time-block holdout geometry -------------------------------------
# Five expanding-window blocks: block b fits on the contiguous prefix
# y[, 1:block_seams[b]] and forecasts the next block_H steps. Seams are
# spread across the record so the five seam states are near-independent.
block_H <- 15L
block_seams <- c(50L, 80L, 110L, 140L, 170L)
stopifnot(max(block_seams) + block_H <= nt)

# ---- storage ----------------------------------------------------------
y <- array(NA, dim = c(ns, nt, nsets, nsettings))
tau.t <- z.t <- knots.t <- vector("list", length = nsettings)
phi.path <- vector("list", length = nsettings)

# ---- generation loop --------------------------------------------------
for (setting in 1:nsettings) {
  pair <- phi.settings[[setting]]
  is.iid <- all(pair == 0)
  phi.arg <- if (is.iid) 0 else pair

  tau.t.setting <- array(NA, dim = c(nknots.t, nt, nsets))
  z.t.setting <- array(NA, dim = c(nknots.t, nt, nsets))
  knots.t.setting <- array(NA, dim = c(nknots.t, 2, nt, nsets))

  for (set in 1:nsets) {
    set.seed(setting * 100 + set)
    data <- rpotspatTS_arp(
      nt = nt, x = x, s = s, beta = beta.t, gamma = gamma.t, nu = nu.t,
      rho = rho.t, lambda = lambda.t, tau.alpha = tau.alpha.t,
      tau.beta = tau.beta.t, nknots = nknots.t, dist = dist.t,
      phi.z = phi.arg, phi.w = phi.arg, phi.tau = phi.arg
    )

    y[, , set, setting] <- data$y
    tau.t.setting[, , set] <- data$tau
    z.t.setting[, , set] <- data$z
    knots.t.setting[, , , set] <- data$knots
  }

  tau.t[[setting]] <- tau.t.setting
  z.t[[setting]] <- z.t.setting
  knots.t[[setting]] <- knots.t.setting
  phi.path[[setting]] <- list(
    type = if (is.iid) "iid" else "fixed",
    label = setting.label[setting],
    phi1 = rep(pair[1], nt),
    phi2 = rep(pair[2], nt),
    phi_pair = pair
  )
  cat(sprintf(
    "setting %d (%s): phi = (%.2f, %.2f) -- %d datasets generated\n",
    setting, setting.label[setting], pair[1], pair[2], nsets
  ))
}

# ---- save -------------------------------------------------------------
save(y, tau.t, z.t, knots.t, phi.path,
  ns, nt, s, nsets, nsettings,
  block_H, block_seams, setting.label,
  x, # covariate data, identical for all datasets
  file = "simdata.RData"
)

cat(sprintf(
  "\nWrote simdata.RData: y is [%d sites x %d times x %d sets x %d settings]\n",
  ns, nt, nsets, nsettings
))
cat(sprintf(
  "Time-block holdout: H = %d, seams at t = %s\n",
  block_H, paste(block_seams, collapse = ", ")
))
