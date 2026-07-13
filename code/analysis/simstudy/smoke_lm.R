#########################################################################
# smoke_lm.R  --  fast (small-iters) smoke test of the fitting pipeline on
# the long-memory setting 17. Confirms that the i.i.d. (method 1), AR(2)
# (method 7) and AR(1) (method 9) analysis models all FIT without error on
# ARFIMA data, and that phi has the expected shape (AR(2) -> 2 cols, AR(1) ->
# 1 col). Writes real result files to results/17-<method>-1.RData so the
# post-fit pipeline (scores/tables/posterior/plots) can be exercised.
#
# NOTE: iters here are deliberately tiny -- these are NOT publication fits.
# Re-run run-settings.R with full iters for the real study.
#
#   & $R smoke_lm.R
#########################################################################

source("./ar2_load.R", chdir = TRUE)
source("./helpers.R", chdir = TRUE)
options(warn = 1)

setting <- 17L
dataset_id <- 1L
methods <- c(1L, 7L, 9L)
results_dir <- "results"
if (!dir.exists(results_dir)) dir.create(results_dir)

iters <- 600
burn <- 200
update <- 200

catalog <- get_simstudy_method_catalog(include_maxstable = TRUE)

obs <- c(rep(TRUE, nrow(y) - ntest), rep(FALSE, ntest))
y.d <- y[, , dataset_id, setting]
y.o <- y.d[obs, ]
x.o <- x[obs, , , drop = FALSE]
s.o <- s[obs, , drop = FALSE]
x.p <- x[!obs, , , drop = FALSE]
s.p <- s[!obs, , drop = FALSE]

for (m in methods) {
  spec <- catalog[catalog$method_id == m, , drop = FALSE]
  set.seed(get_simstudy_seed(setting, m, dataset_id))
  cat(sprintf("\n=== method %d (%s) ===\n", m, spec$label[1]))

  fit <- mcmc(
    y = y.o, x = x.o, s = s.o, s.pred = s.p, x.pred = x.p,
    method = spec$method[1], skew = isTRUE(spec$skew[1]),
    thresh.all = spec$thresh_all[1], thresh.quant = isTRUE(spec$thresh_quant[1]),
    nknots = spec$nknots[1], iterplot = FALSE,
    iters = iters, burn = burn, update = update,
    min.s = c(0, 0), max.s = c(10, 10),
    temporalw = isTRUE(spec$temporalw[1]),
    temporaltau = isTRUE(spec$temporaltau[1]),
    temporalz = isTRUE(spec$temporalz[1]),
    ar2_w = isTRUE(spec$ar2_w[1]),
    ar2_tau = isTRUE(spec$ar2_tau[1]),
    ar2_z = isTRUE(spec$ar2_z[1]),
    rho.upper = 15, nu.upper = 10
  )

  # shape / sanity report
  ptau <- fit$phi.tau
  if (is.null(ptau)) {
    cat("  phi.tau: NULL (i.i.d. in time)\n")
  } else {
    cat(sprintf("  phi.tau: %d col(s); posterior mean = %s\n",
                ncol(ptau), paste(sprintf("%.3f", colMeans(ptau)), collapse = ", ")))
    if (ncol(ptau) == 2) {
      cat(sprintf("  -> AR(2) phi2 posterior mean = %.3f (YW target for d=0.20: %.3f)\n",
                  mean(ptau[, 2]), arfima_yw_projection(0.20)$ar2["phi2"]))
    }
  }

  outfile <- build_simstudy_result_file(results_dir, setting, m, dataset_id)
  fit.1 <- fit                       # scores.R / posterior.R expect object 'fit.1'
  save(fit.1, file = outfile)
  cat(sprintf("  saved -> %s\n", outfile))
}

cat("\nSmoke test complete.\n")
