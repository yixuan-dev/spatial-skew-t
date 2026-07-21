# convergence/04_multichain_run.R
# ---------------------------------------------------------------------------
# Run ONE additional MCMC chain on the fold-1 training data of a headline
# setting, with a fresh seed and overdispersed initial values, for the
# multi-chain Gelman-Rubin check in 05_multichain_diagnostics.R.
#
# Chain 1 is the production fold-1 chain (seed = setting*100 + 1, default
# inits) and is never rerun here; valid chain ids are 2-4. The data slicing
# and call construction below mirror us-all-run.R lines 378-491 at commit
# 51fa383; the stopifnot() block pins the invariants so silent drift in
# either file fails loudly instead of quietly targeting a different
# posterior.
#
# Usage (from the US-all directory, one process per chain):
#   Rscript convergence/04_multichain_run.R <setting> <chain> [dev|prod]
# or via env vars US_ALL_CONV_SETTING / US_ALL_CONV_CHAIN / US_ALL_RUN_MODE.
# US_ALL_CONV_PRED=1 additionally draws held-out predictions as production
# did; default off, which leaves the parameter chains' target unchanged
# (both samplers draw predictions only when iter > burn and never feed them
# back into parameter updates).
#
# Output: output/us-all/results/convergence/multichain-<setting>-fold1-chain<c>.rds
#         (slim: scalar chains + metadata only, a few MB)
# ---------------------------------------------------------------------------

rm(list = ls())

library(compiler)
enableJIT(3)

.this <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), winslash = "/", mustWork = FALSE)) else "."
})
setwd(file.path(.this, ".."))
source("convergence/00_conv_lib.R")

# --- arguments -------------------------------------------------------------
cli <- commandArgs(trailingOnly = TRUE)
setting <- suppressWarnings(as.integer(
  if (length(cli) >= 1) cli[1] else Sys.getenv("US_ALL_CONV_SETTING", unset = "")
))
chain <- suppressWarnings(as.integer(
  if (length(cli) >= 2) cli[2] else Sys.getenv("US_ALL_CONV_CHAIN", unset = "")
))
RUN_MODE <- tolower(
  if (length(cli) >= 3) cli[3] else Sys.getenv("US_ALL_RUN_MODE", unset = "prod")
)
with_pred <- Sys.getenv("US_ALL_CONV_PRED", unset = "0") %in% c("1", "true", "yes")

if (is.na(setting) || is.na(chain)) {
  stop("usage: Rscript convergence/04_multichain_run.R <setting> <chain 2-4> [dev|prod]",
       call. = FALSE)
}
if (!setting %in% conv_multichain_settings) {
  stop(sprintf("setting %d has no multichain spec (supported: %s); MRTS setting 204 is reserved but not implemented, see README",
               setting, paste(conv_multichain_settings, collapse = ", ")), call. = FALSE)
}
if (!chain %in% 2:4) stop("chain must be 2, 3, or 4 (chain 1 = production fit)", call. = FALSE)
if (!RUN_MODE %in% c("dev", "prod")) stop("run mode must be dev or prod", call. = FALSE)

mcmc_ctrl <- switch(RUN_MODE,
  dev = list(iters = 2000, burn = 1000, update = 200),
  prod = list(iters = 30000, burn = 25000, update = 500)
)

row <- conv_settings[conv_settings$setting == setting, , drop = FALSE]
backend <- row$backend                      # 55 -> legacy, 111 -> ar2
inits <- conv_chain_inits[[as.character(chain)]]
seed <- conv_chain_seed(setting, chain)

cat("==============================\n")
cat(sprintf("multichain rerun | setting %d | chain %d | backend %s | %s\n",
            setting, chain, backend, RUN_MODE))
cat(sprintf("seed %d | iters %d | burn %d | predictions %s\n",
            seed, mcmc_ctrl$iters, mcmc_ctrl$burn, if (with_pred) "on" else "off"))
cat("inits:", paste(sprintf("%s=%s", names(inits),
    vapply(inits, function(x) paste(x, collapse = ","), "")), collapse = " "), "\n")

# --- backend + data (mirrors us-all-run.R:42-63) ---------------------------
backend_env <- new.env(parent = globalenv())
if (backend == "ar2") {
  sys.source("./ar2_load.R", envir = backend_env, chdir = TRUE)
} else {
  sys.source("./package_load.R", envir = backend_env, chdir = TRUE)
}
Y_data <- get("Y", envir = backend_env)
X_data <- get("X", envir = backend_env)
S_data <- get("S", envir = backend_env)
cv_folds <- get("cv.lst", envir = backend_env)
beta_init_default <- get("beta.init", envir = backend_env)
# Both loaders define the sampler as `mcmc` in their env (ar2_load.R sources
# mcmc_ar2.R last, overriding any earlier `mcmc`); the formals distinguish
# the AR2 sampler (has ar2_tau) from the legacy one.
run_mcmc <- get("mcmc", envir = backend_env)
if (backend == "ar2") {
  stopifnot("ar2_tau" %in% names(formals(run_mcmc)))
} else {
  stopifnot(!"ar2_tau" %in% names(formals(run_mcmc)),
            "temporaltau" %in% names(formals(run_mcmc)))
}

# --- fold-1 slicing (mirrors us-all-run.R:378-398; all 9 settings use CMAQ,
# none of the multichain settings use MRTS) --------------------------------
set.seed(seed)
val.idx <- cv_folds[[1]]
y.o <- Y_data[-val.idx, ]
X.o <- X_data[-val.idx, , , drop = FALSE]
S.o <- S_data[-val.idx, ]
X.p <- X_data[val.idx, , , drop = FALSE]
S.p <- S_data[val.idx, ]

# --- call construction (mirrors us-all-run.R:452-491) ----------------------
call_common <- list(
  y = y.o,
  s = S.o,
  x = X.o,
  method = row$method,
  skew = row$skew,
  keep.knots = FALSE,
  thresh.all = row$thresh,
  thresh.quant = FALSE,
  nknots = row$nknots,
  iters = mcmc_ctrl$iters,
  burn = mcmc_ctrl$burn,
  update = mcmc_ctrl$update,
  iterplot = FALSE,
  beta.init = beta_init_default,
  tau.init = inits$tau.init,          # production: pick_tau_init("t") = 0.05
  tau.alpha.init = inits$tau.alpha.init,
  tau.beta.init = inits$tau.beta.init,
  gamma.init = inits$gamma.init,
  rho.init = inits$rho.init,
  rho.upper = 5,
  nu.init = inits$nu.init,
  nu.upper = 10,
  z.init = inits$z.init,
  min.s = c(-2.25, -1.55),
  max.s = c(2.35, 1.30)
)
if (row$skew) call_common$lambda.init <- inits$lambda.init
if (with_pred) {
  call_common$x.pred <- X.p
  call_common$s.pred <- S.p
}
if (row$temporal) {
  call_common$temporaltau <- TRUE
  call_common$temporalw <- TRUE
  call_common$temporalz <- TRUE
}
if (backend == "ar2") {
  call_common$ar2_tau <- isTRUE(row$temporal)
  call_common$ar2_w <- isTRUE(row$temporal)
  call_common$ar2_z <- isTRUE(row$temporal)
  call_common$phi.tau.init <- inits$phi
  call_common$phi.w.init <- inits$phi
  call_common$phi.z.init <- inits$phi
}

# Pin the invariants of the production run this must replicate.
stopifnot(
  identical(dim(y.o), c(400L, 31L)),
  identical(dim(X.o), c(400L, 31L, 2L)),
  identical(dim(S.o), c(400L, 2L)),
  call_common$thresh.all == row$thresh,
  call_common$nknots == row$nknots,
  RUN_MODE != "prod" || (call_common$iters == 30000 && call_common$burn == 25000),
  # AR2 phi inits must sit in the stationarity triangle
  backend != "ar2" || all(abs(inits$phi[2]) < 1,
                          inits$phi[1] + inits$phi[2] < 1,
                          inits$phi[2] - inits$phi[1] < 1)
)

# --- run -------------------------------------------------------------------
t0 <- proc.time()[3]
fit <- do.call(run_mcmc, call_common)
elapsed <- proc.time()[3] - t0
cat(sprintf("chain finished in %.1f min\n", elapsed / 60))

# --- slim + save -----------------------------------------------------------
chains <- extract_scalar_chains(fit)
backend_detected <- detect_backend(fit)
rm(fit); invisible(gc(verbose = FALSE))

paths <- conv_paths()
out <- multichain_file(setting, chain, paths)
saveRDS(list(
  setting = setting,
  fold = 1L,
  chain = chain,
  seed = seed,
  run_mode = RUN_MODE,
  backend = backend,
  backend_detected = backend_detected,
  inits = inits,
  with_pred = with_pred,
  mcmc_control = mcmc_ctrl,
  chains = chains,
  elapsed_sec = unname(elapsed),
  r_version = R.version.string,
  finished_at = Sys.time()
), out)
cat(sprintf("wrote %s (%.2f MB, %d draws x %d params)\n",
            out, file.size(out) / 2^20, nrow(chains), ncol(chains)))
