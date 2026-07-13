#########################################################################
# setup_lm_append.R  --  idempotent extender that ADDS the long-memory
# settings 17-19 to an existing simdata.RData WITHOUT regenerating settings
# 1-16 (avoids re-running the slow max-stable / Brown-Resnick generators and
# their SpatialExtremes dependency).
#
# Settings 17-19 are byte-identical to what a full setup.R run would produce:
# same rpotspatTS_arp generator, same seed = setting*100 + set, same scalar
# constants. Consistency is ASSERTED by replaying setting 9 and requiring it to
# match the stored y[,,,9] exactly -- which also confirms the make*ARP -> spec
# refactor left the numeric AR path bit-for-bit unchanged.
#
# setup.R remains the source of truth for the design.
#
#   $R = "C:\Program Files\R\R-4.5.1\bin\Rscript.exe"
#   & $R setup_lm_append.R
#########################################################################

source("ar2_load.R")   # loads simdata.RData + backend (rpotspatTS_arp, arfima_yw_projection)
options(warn = 1)      # generation should not warn; do not treat warnings as errors here

## ---- scalar constants: MUST match setup.R (asserted below via setting 9) ----
beta.t <- c(10, 0, 0)
nu.t <- 0.5
gamma.t <- 0.9
tau.alpha.t <- 6
tau.beta.t <- 16
rho.gen <- 1      # rho.t[4]
lambda.gen <- 3   # lambda.t[4]
dist.gen <- "t"   # dist.t[4]

ns <- dim(y)[1]
stopifnot(dim(y)[2] == nt, dim(y)[3] == nsets)

# ARFIMA d per setting, mirroring build_phi_triple() in setup.R.
lm_d <- c("17" = 0.20, "18" = 0.35, "19" = 0.45)
lm_spec <- function(d) list(type = "arfima", d = d)

# Generate one (setting) worth of datasets, mirroring setup.R's setting>=9 block
# exactly (array shapes, assignment order, seed). `triple` gives the per-channel
# phi/arfima spec.
gen_setting <- function(setting, triple, nknots) {
  y.a <- array(NA, dim = c(ns, nt, nsets))
  tau.a <- array(NA, dim = c(nknots, nt, nsets))
  z.a <- array(NA, dim = c(nknots, nt, nsets))
  knots.a <- array(NA, dim = c(nknots, nt, 2, nsets))
  for (set in 1:nsets) {
    set.seed(setting * 100 + set)
    data <- rpotspatTS_arp(
      nt = nt, x = x, s = s, beta = beta.t, gamma = gamma.t, nu = nu.t,
      rho = rho.gen, lambda = lambda.gen, tau.alpha = tau.alpha.t,
      tau.beta = tau.beta.t, nknots = nknots, dist = dist.gen,
      phi.z = triple$phi.z, phi.w = triple$phi.w, phi.tau = triple$phi.tau
    )
    y.a[, , set] <- data$y
    tau.a[, , set] <- data$tau
    z.a[, , set] <- data$z
    knots.a[, , , set] <- data$knots
  }
  list(y = y.a, tau = tau.a, z = z.a, knots = knots.a)
}

## ---- consistency gate: replay setting 9 and require an exact match ----------
cat("Consistency check: replaying setting 9 (AR(2) phi = (0.80, -0.35))...\n")
phi_strong <- c(0.80, -0.35)
chk <- gen_setting(9, list(phi.z = phi_strong, phi.tau = phi_strong, phi.w = phi_strong), 1)
max_diff <- max(abs(chk$y - y[, , , 9]))
if (!is.finite(max_diff) || max_diff > 1e-8) {
  stop(sprintf(
    "Setting-9 replay does NOT match stored y[,,,9] (max abs diff = %.3e); constants or generator drifted. Aborting.",
    max_diff))
}
cat(sprintf("  OK: setting-9 replay matches (max abs diff = %.3e) -> constants match, AR path preserved.\n\n", max_diff))

## ---- extend the storage objects to hold settings 17-19 ---------------------
target_ns <- 19L
if (dim(y)[4] < target_ns) {
  new_y <- array(NA, dim = c(ns, nt, nsets, target_ns))
  new_y[, , , seq_len(dim(y)[4])] <- y
  y <- new_y
}

## ---- generate long-memory settings 17-19 -----------------------------------
for (setting in c(17L, 18L, 19L)) {
  d <- lm_d[[as.character(setting)]]
  spec <- lm_spec(d)
  triple <- list(phi.z = spec, phi.tau = spec, phi.w = spec)
  out <- gen_setting(setting, triple, 1)

  y[, , , setting] <- out$y
  tau.t[[setting]] <- out$tau
  z.t[[setting]] <- out$z
  knots.t[[setting]] <- out$knots

  yw <- arfima_yw_projection(d)
  phi.path[[setting]] <- list(
    type = "arfima", d = d, hurst = yw$hurst,
    phi.z = spec, phi.tau = spec, phi.w = spec,
    yw_ar1 = yw$ar1, yw_ar2 = yw$ar2
  )
  cat(sprintf("setting %d: ARFIMA d = %.2f (Hurst %.2f), YW AR(2) phi2 = %.3f -- %d datasets\n",
              setting, d, yw$hurst, yw$ar2["phi2"], nsets))
}

## ---- save (mirrors setup.R's save signature) -------------------------------
save(y, tau.t, z.t, knots.t, phi.path, ns, nt, s, nsets, ntest, x,
  file = "simdata.RData")
cat(sprintf("\nsimdata.RData updated: dim(y) = [%s]\n", paste(dim(y), collapse = ", ")))
