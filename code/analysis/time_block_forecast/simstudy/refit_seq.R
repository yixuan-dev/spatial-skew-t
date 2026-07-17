#########################################################################
# refit_seq.R -- FAITHFULLY reproduce the original run's RNG stream.
#
# Why: a standalone refit of (setting 4, dataset 1, block 3) came out
# perfectly healthy -- recursion marginal sd -> 1, sigma ~ 1.4-1.8, no clamp,
# no shape collapse -- yet the ORIGINAL run's block 3 exploded (predictive SD
# 104 vs a data marginal SD of 4.38). The only difference: run-settings.R
# calls set.seed() ONCE and then walks blocks 1..5 in order, so by block 3 the
# RNG is deep into a different stream. This script replicates that exactly
# (same seed, same update=1000, same mcmc -> forecast_block order per block)
# but SAVES every fit.
#
# In the original run: block 1 was CLEAN, block 3 EXPLODED. Reproducing that
# contrast inside one chain gives a controlled comparison -- same data, same
# chain, one good block and one bad one -- so we can see exactly which
# quantity goes wrong, instead of guessing.
#
#   & $R refit_seq.R
#########################################################################

source("../../simstudy/ar2_load.R", chdir = TRUE)
source("./time_block_helpers.R")
load("./simdata.RData", envir = .GlobalEnv)
options(warn = 1)

setting <- 4L
dataset_id <- 1L
method_id <- 2L          # AR(2)
nblocks_run <- 3L        # block 1 (clean) .. block 3 (exploded)
iters <- 20000
burn <- 10000
update <- 1000           # MATCH the original exactly
K <- 1L

blocks <- tbf_blocks(block_seams, block_H, nt)
y.d <- y[, , dataset_id, setting]
marg_sd <- sd(y[, , , setting])
cat(sprintf("data marginal sd = %.2f\n", marg_sd))
cat(sprintf("original run: block1 predictive SD ~5 (clean), block3 ~104 (EXPLODED)\n\n"))

# ---- replicate run_method()'s RNG stream exactly -----------------------
set.seed(get_tbf_seed(setting, method_id, dataset_id))

fits <- vector("list", nblocks_run)
for (b in seq_len(nblocks_run)) {
  blk <- blocks[[b]]
  y.train <- y.d[, blk$train_times, drop = FALSE]
  x.train <- x[, blk$train_times, , drop = FALSE]
  x.block <- x[, blk$test_times, , drop = FALSE]

  fit <- mcmc(
    y = y.train, x = x.train, s = s,
    method = "t", skew = TRUE,
    thresh.all = 0, thresh.quant = TRUE,
    nknots = 1,
    iterplot = FALSE, iters = iters, burn = burn, update = update,
    min.s = c(0, 0), max.s = c(10, 10),
    temporalw = TRUE, temporaltau = TRUE, temporalz = TRUE,
    ar2_w = TRUE, ar2_tau = TRUE, ar2_z = TRUE,
    z.init = 0.01,
    rho.upper = 15, nu.upper = 10
  )

  # forecast_block consumes RNG too -- must run it to stay aligned
  yhat <- forecast_block(
    fit = fit, seam = blk$seam, H = block_H,
    x_block = x.block, s = s, ar2 = TRUE
  )

  fits[[b]] <- list(fit = fit, blk = blk, yhat_sd = sapply(1:block_H,
                    function(h) sd(as.numeric(yhat[, , h]))))
  cat(sprintf("block %d (seam %3d) done | predictive SD: lead1 %7.2f  lead15 %8.2f  %s\n",
              b, blk$seam, fits[[b]]$yhat_sd[1], fits[[b]]$yhat_sd[block_H],
              if (fits[[b]]$yhat_sd[block_H] > 3 * marg_sd) "<<< EXPLODED" else "(clean)"))
  rm(yhat); gc()
}
save(fits, file = "refit_seq_fits.RData")
cat("\nfits saved -> refit_seq_fits.RData\n")

# ---- side-by-side diagnosis: what differs between the blocks? ----------
cat("\n================ SIDE-BY-SIDE DIAGNOSIS ================\n")
cat(sprintf("%-22s %12s %12s %12s\n", "quantity", "block 1", "block 2", "block 3"))
row <- function(nm, f) {
  v <- sapply(fits, f)
  cat(sprintf("%-22s %12s %12s %12s\n", nm,
              formatC(v[1], format = "g", digits = 4),
              formatC(v[2], format = "g", digits = 4),
              formatC(v[3], format = "g", digits = 4)))
}
row("pred SD lead 1",   function(o) o$yhat_sd[1])
row("pred SD lead 15",  function(o) o$yhat_sd[block_H])
row("a (shape) mean",   function(o) mean(o$fit$tau.alpha) / 2)
row("a min",            function(o) min(o$fit$tau.alpha) / 2)
row("b (rate) mean",    function(o) mean(o$fit$tau.beta) / 2)
row("lambda mean",      function(o) mean(o$fit$lambda))
row("tau min",          function(o) min(o$fit$tau))
row("tau at 1e-6 clamp",function(o) sum(o$fit$tau <= 1e-6 + 1e-12))
row("z min",            function(o) min(o$fit$z))

# seam state on the Gaussian scale + recursion behaviour, per block
seam_stats <- function(o) {
  fit <- o$fit; blk <- o$blk
  a <- fit$tau.alpha / 2; b <- fit$tau.beta / 2
  sc <- c(blk$seam - 1L, blk$seam)
  ts <- zs <- numeric(nrow(fit$tau))
  TSsd <- numeric(block_H)
  TS <- matrix(NA_real_, nrow(fit$tau), block_H)
  for (m in seq_len(nrow(fit$tau))) {
    tobs <- matrix(fit$tau[m, sc], K, 2L)
    zobs <- matrix(fit$z[m, sc], K, 2L)
    tsm <- gamma.cop(tobs, a[m], b[m])
    zsm <- hn.cop(zobs, sig = 1 / sqrt(tobs))
    ts[m] <- tsm[1, 2]; zs[m] <- zsm[1, 2]
    TS[m, ] <- ar2_recurse(tsm, fit$phi.tau[m, ], block_H, K)
  }
  tau.fc <- qgamma(pnorm(TS), rep(a, block_H), rep(b, block_H))
  sig <- 1 / sqrt(tau.fc)
  list(ts_mean = mean(ts[is.finite(ts)]), ts_sd = sd(ts[is.finite(ts)]),
       ts_min = min(ts[is.finite(ts)]), ts_nonfin = sum(!is.finite(ts)),
       zs_nonfin = sum(!is.finite(zs)),
       rec_sd15 = sd(TS[, block_H], na.rm = TRUE),
       sig_med = median(sig[, block_H], na.rm = TRUE),
       sig_max = max(sig[, block_H], na.rm = TRUE))
}
S <- lapply(fits, seam_stats)
cat("\n-- seam state on the Gaussian copula scale (should look like N(0,1)) --\n")
for (nm in c("ts_mean", "ts_sd", "ts_min", "ts_nonfin", "zs_nonfin",
             "rec_sd15", "sig_med", "sig_max")) {
  cat(sprintf("%-22s %12s %12s %12s\n", nm,
              formatC(S[[1]][[nm]], format = "g", digits = 4),
              formatC(S[[2]][[nm]], format = "g", digits = 4),
              formatC(S[[3]][[nm]], format = "g", digits = 4)))
}
cat("\nrec_sd15 = marginal sd of ar2_recurse output at lead 15 (MUST be ~1)\n")
cat("sig_*    = 1/sqrt(tau.fc) at lead 15; data marginal sd is", round(marg_sd, 2), "\n")
