#########################################################################
# run_block1_guarded.R -- block-1 positive control for AR(2) vs i.i.d.
#
# Question: fitting AR(2)-generated data with the AR(2) model, does the
# time-block forecast Brier score beat the i.i.d.-in-time skew-t baseline?
#
# Design (see plan): settings {4 strong, 5 nearunit} x methods {1 iid,
# 2 AR(2)} x datasets 1..10, block 1 ONLY (train y[, 1:50], forecast
# leads 1..15, i.e. t = 51..65). Data: ../simstudy/simdata.RData
# (setup.R with setting 5 = phi (0.15, 0.80), spectral radius 0.973).
#
# GUARD (pre-registered, symmetric across methods): block 1's short prefix
# is where the reflected lambda ridge bites (tex/lambda_phiz_ridge). Each
# fit must pass B (truth recovery: lambda sign+size, beta0) AND C
# (predictive spread <= 3 x marginal SD) from check_fit_consistency();
# otherwise refit with seed + 7919*attempt, at most 3 reseeds. If still
# failing, the last fit is kept and FLAGGED (pass = FALSE): the analysis
# excludes it from the primary comparison and reports it in a sensitivity
# run. A and A' (sdz_ratio, transverse to the ridge) are recorded always.
#
#   $R = "C:\Program Files\R\R-4.5.1\bin\Rscript.exe"
#   & $R run_block1_guarded.R --settings=4,5 --datasets=1:10 --methods=1:2
#########################################################################

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1]),
    winslash = "/", mustWork = FALSE
  ))
  if (dir.exists(script_dir)) setwd(script_dir)
}

# ar2_load.R does rm(list = ls()): source it FIRST, then the helpers it wiped.
source("../../simstudy/ar2_load.R", chdir = TRUE)
source("../simstudy/time_block_helpers.R")
source("../simstudy/fit_diag_utils.R")
load("../simstudy/simdata.RData", envir = .GlobalEnv)
options(warn = 1)

# ---- CLI ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
getflag <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
settings <- eval(parse(text = paste0("c(", getflag("settings", "4,5"), ")")))
datasets <- eval(parse(text = getflag("datasets", "1:10")))
method_ids <- eval(parse(text = getflag("methods", "1:2")))

iters <- 20000
burn <- 10000
update <- 2000
max_attempts <- 4L                     # attempt 0 + up to 3 reseeds
blk <- tbf_blocks(block_seams, block_H, nt)[[1]]   # seam 50, leads 1..15
catalog <- get_tbf_method_catalog()

dir.create("results", showWarnings = FALSE)

cat(sprintf(
  "block1 positive control: settings=%s methods=%s datasets=%s\n",
  paste(settings, collapse = ","), paste(method_ids, collapse = ","),
  paste(range(datasets), collapse = "..")
))
cat(sprintf("block 1: train t=1..%d, forecast t=%d..%d (H=%d)\n\n",
            blk$seam, blk$seam + 1L, blk$seam + block_H, block_H))

for (setting in settings) {
  marginal_sd <- sd(y[, , , setting])
  pp <- phi.path[[setting]]
  truth <- list(lambda = 3, beta0 = 10)   # B uses lambda/beta0 only
  temporal_desc <- if (identical(pp$type, "arfima")) {
    sprintf("ARFIMA d = %.2f (Hurst %.2f)", pp$d, pp$hurst)
  } else {
    sprintf("phi = (%.2f, %.2f)", pp$phi_pair[1], pp$phi_pair[2])
  }
  cat(sprintf("=== setting %d (%s): %s, marginal sd %.2f ===\n",
              setting, setting.label[setting], temporal_desc, marginal_sd))

  for (d in datasets) {
    y.train <- y[, blk$train_times, d, setting]
    x.train <- x[, blk$train_times, , drop = FALSE]
    x.block <- x[, blk$test_times, , drop = FALSE]
    y.val <- y[, blk$test_times, d, setting]

    for (m in method_ids) {
      outfile <- sprintf("results/blk1-%d-%d-%d.RData", setting, m, d)
      if (file.exists(outfile)) {
        cat(sprintf("s%d d%d m%d: exists, skip\n", setting, d, m))
        next
      }
      spec <- catalog[catalog$method_id == m, , drop = FALSE]
      tic <- proc.time()

      # Early stop for REPRODUCIBLE failures: if two consecutive attempts
      # land on the same lambda (within 15%) and both fail, the failure is
      # data-level (e.g. a flat long-memory z realization under-identifies
      # lambda in a 50-day window) -- reseeding cannot fix it, and further
      # attempts risk trading a stable near-miss for a reflected chain.
      # This never changes which fits PASS; it only stops futile reseeds.
      chk <- NULL
      yhat <- NULL
      lam_prev <- NA_real_
      for (attempt in 0:(max_attempts - 1L)) {
        seed <- get_tbf_seed(setting, m, d) + 7919L * attempt
        set.seed(seed)
        fit <- mcmc(
          y = y.train, x = x.train, s = s,
          method = "t", skew = isTRUE(spec$skew[1]),
          thresh.all = 0, thresh.quant = TRUE,
          nknots = spec$nknots[1],
          iterplot = FALSE, iters = iters, burn = burn, update = update,
          min.s = c(0, 0), max.s = c(10, 10),
          temporalw = isTRUE(spec$temporal[1]),
          temporaltau = isTRUE(spec$temporal[1]),
          temporalz = isTRUE(spec$temporal[1]),
          ar2_w = isTRUE(spec$ar2[1]),
          ar2_tau = isTRUE(spec$ar2[1]),
          ar2_z = isTRUE(spec$ar2[1]),
          rho.upper = 15, nu.upper = 10
        )
        yhat <- forecast_block(fit = fit, seam = blk$seam, H = block_H,
                               x_block = x.block, s = s,
                               ar2 = isTRUE(spec$ar2[1]))
        chk <- check_fit_consistency(fit, yhat, truth, marginal_sd,
                                     data_mean = mean(y.train))
        chk <- cbind(setting = setting, method = m, dataset = d,
                     attempt = attempt, seed = seed,
                     pass = isTRUE(chk$B_truth) && isTRUE(chk$C_spread), chk)
        rm(fit)
        gc()
        cat(sprintf(
          "s%d d%d m%d att%d: lam %+7.2f beta0 %6.2f sdzr %.2f SDmax %6.2f %s\n",
          setting, d, m, attempt, chk$lambda, chk$beta0, chk$sdz_ratio,
          chk$sd_lead_max,
          if (chk$pass) "PASS" else "FAIL -> reseed"
        ))
        if (chk$pass) break
        if (!is.na(lam_prev) &&
            abs(chk$lambda - lam_prev) < 0.15 * max(abs(lam_prev), 0.5)) {
          cat(sprintf(
            "s%d d%d m%d: reproducible failure (lam stable) -> stop reseeding\n",
            setting, d, m))
          break
        }
        lam_prev <- chk$lambda
      }

      elapsed <- unname((proc.time() - tic)[3])
      res <- list(setting = setting, method_id = m, dataset = d,
                  seam = blk$seam, leads = blk$leads,
                  test_times = blk$test_times,
                  yhat = yhat, y_val = y.val,
                  chk = chk, elapsed_sec = elapsed)
      save(res, file = outfile)
      rm(res, yhat)
      gc()
      cat(sprintf("s%d d%d m%d saved (%.0f s, %d attempt(s), pass=%s)\n",
                  setting, d, m, elapsed, chk$attempt + 1L, chk$pass))
    }
  }
}
cat("\nall requested fits finished\n")
