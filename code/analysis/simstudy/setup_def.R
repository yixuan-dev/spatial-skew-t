#########################################################################
# Simulation study with deformed-covariance settings.
# All settings: skew-t-1 (lambda = 3, dist = "t", nknots = 1).
#
# Settings 1–3: geometric anisotropy (linear deformation via CorFxDef)
#   1 - theta = 0,    ratio = 1.00   (isotropic baseline)
#   2 - theta = pi/4, ratio = 0.50   (moderate anisotropy, 45 deg)
#   3 - theta = pi/6, ratio = 0.25   (strong anisotropy, 30 deg)
#
# Settings 4–6: nonlinear deformation (Sampson & Guttorp 1992, via CorFxNL)
#   4 - nl_axial    D(s1,s2) = (s1, exp(0.15*s1)*s2)
#   5 - nl_sine     D(s1,s2) = (s1 + 1.5*sin(0.4*s2), s2)
#   6 - nl_composed D = axial(alpha=0.10) o rotation(pi/6)
#########################################################################

rm(list = ls())

library(fields)
library(SpatialTools)

source("../../R/ar2/auxfunctions.R")

# Pure-R override of mem() to avoid the Rcpp build of g.Rcpp.
mem <- function(s, knots) {
  d <- fields::rdist(s, knots)
  apply(d, 1, which.min)
}

# ── Common parameters ────────────────────────────────────────────────────────
beta.t      <- c(10, 0, 0)
nu.t        <- 0.5
gamma.t     <- 0.9
rho.t       <- 1
tau.alpha.t <- 3
tau.beta.t  <- 8

# ── Per-setting parameters (6 settings) ──────────────────────────────────────
nsettings   <- 6
deform_type <- c("linear", "linear", "linear",
                 "nl_axial", "nl_sine", "nl_composed")
dist.def    <- rep("t", nsettings)
nknots.def  <- rep(1L, nsettings)
lambda.def  <- rep(3, nsettings)

# Linear-deformation parameters (NA for NL settings)
theta.def <- c(0, pi / 4, pi / 6, NA, NA, NA)
ratio.def <- c(1.00, 0.50, 0.25,  NA, NA, NA)

# ── Sites & covariates ────────────────────────────────────────────────────────
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

# ── Nonlinear deformation functions (settings 4–6) ───────────────────────────
D_fns <- vector("list", nsettings)

D_fns[[4]] <- function(s) {
  cbind(s[, 1], exp(0.15 * s[, 1]) * s[, 2])
}

D_fns[[5]] <- function(s) {
  cbind(s[, 1] + 1.5 * sin(0.4 * s[, 2]), s[, 2])
}

D_rot <- function(s) {
  th <- pi / 6
  cbind(cos(th) * s[, 1] - sin(th) * s[, 2],
        sin(th) * s[, 1] + cos(th) * s[, 2])
}
D_ax10 <- function(s) cbind(s[, 1], exp(0.10 * s[, 1]) * s[, 2])
D_fns[[6]] <- function(s) D_ax10(D_rot(s))

# ── Jacobian check (settings 4–6 must be bijective on domain) ────────────────
jac_det <- function(D, s_grid, h = 1e-5) {
  e1 <- cbind(rep(h, nrow(s_grid)), 0)
  e2 <- cbind(0, rep(h, nrow(s_grid)))
  g1 <- (D(s_grid + e1) - D(s_grid - e1)) / (2 * h)
  g2 <- (D(s_grid + e2) - D(s_grid - e2)) / (2 * h)
  g1[, 1] * g2[, 2] - g1[, 2] * g2[, 1]
}
grid_chk <- as.matrix(expand.grid(s1 = seq(0, 10, 0.5),
                                  s2 = seq(0, 10, 0.5)))
for (setting in 4:nsettings) {
  jd <- jac_det(D_fns[[setting]], grid_chk)
  stopifnot("D has negative Jacobian (folding)" = all(jd > 0))
}

# ── Storage ───────────────────────────────────────────────────────────────────
y            <- array(NA, dim = c(ns, nt, nsets, nsettings))
tau.t        <- vector("list", length = nsettings)
z.t          <- vector("list", length = nsettings)
knots.t      <- vector("list", length = nsettings)
f_coords_list <- vector("list", length = nsettings)

# ── Main loop ─────────────────────────────────────────────────────────────────
for (setting in seq_len(nsettings)) {
  nknots <- nknots.def[setting]
  tau.t.setting   <- array(NA, dim = c(nknots, nt, nsets))
  z.t.setting     <- array(NA, dim = c(nknots, nt, nsets))
  knots.t.setting <- array(NA, dim = c(nknots, nt, 2, nsets))

  nl       <- deform_type[setting] != "linear"
  cov_type <- if (nl) "nl_deformed" else "deformed"
  f_coords <- if (nl) D_fns[[setting]](s) else NULL
  theta_i  <- if (!nl) theta.def[setting] else 0
  ratio_i  <- if (!nl) ratio.def[setting] else 1

  f_coords_list[[setting]] <- f_coords

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
      cov.type  = cov_type,
      theta     = theta_i,
      ratio     = ratio_i,
      f_coords  = f_coords
    )

    y[, , set, setting]        <- data$y
    tau.t.setting[, , set]     <- data$tau
    z.t.setting[, , set]       <- data$z
    knots.t.setting[, , , set] <- data$knots
  }

  cat(sprintf("finished setting %d (%s)\n", setting, deform_type[setting]))

  tau.t[[setting]]   <- tau.t.setting
  z.t[[setting]]     <- z.t.setting
  knots.t[[setting]] <- knots.t.setting
}

# ── Metadata data frame ───────────────────────────────────────────────────────
settings.def <- data.frame(
  setting     = seq_len(nsettings),
  dist        = dist.def,
  nknots      = nknots.def,
  lambda      = lambda.def,
  theta       = theta.def,
  ratio       = ratio.def,
  deform_type = deform_type
)

# ── Save ──────────────────────────────────────────────────────────────────────
save(
  y, tau.t, z.t, knots.t, settings.def, f_coords_list,
  ns, nt, s, nsets, ntest, x,
  file = "simdata_def.RData"
)

cat(sprintf(
  "saved simdata_def.RData: %d settings (%d linear, %d nonlinear), dim(y) = [%s]\n",
  nsettings,
  sum(deform_type == "linear"),
  sum(deform_type != "linear"),
  paste(dim(y), collapse = ", ")
))
