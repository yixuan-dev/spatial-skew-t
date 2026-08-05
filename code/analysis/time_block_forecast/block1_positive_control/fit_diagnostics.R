#########################################################################
# fit_diagnostics.R -- refit time-block cases WITH the fit retained, and run
# self-consistency assertions that turn a silent fit failure into a loud one.
#
# Moved here from ../simstudy on 2026-08-05. It reads the simstudy data
# and refits simstudy cases, but it is an assertion harness, and the
# assertions belong with the study that gates on them: the simstudy arm
# records descriptive diagnostics only (see ../simstudy/fit_diag_utils.R)
# while run_block1_guarded.R gates on B and C. Outputs now land under this
# directory's output/diag/.
#
# Motivation: run-settings.R saves only `forecast`, never `fit`, so lambda / z /
# beta0 / phi are unrecoverable after a run -- and every diagnosis of the
# time-block divergence has been blocked by that. This standalone harness
# (it does NOT modify run-settings.R) refits chosen cases, keeps per-fit
# summaries, and checks four assertions:
#
#   A  z self-consistency : mean(z) ~ mean(1/sqrt(tau)) * sqrt(2/pi)
#                           (the model IS z ~ HalfNormal(0, 1/sqrt(tau)))
#   B  truth recovery     : lambda ~ 3, beta0 ~ 10, phi.z ~ truth  (sim only)
#   C  predictive spread  : sd(yhat[,,h]) <= k * marginal_sd for all h
#                           (a stationary AR forecast converges to climatology)
#   D  ridge check        : report beta0 + lambda*mean(z) vs the data mean
#
# The immediate use is to diagnose the residual explosions concentrated in
# block 4 (seam 140) that survive the z.init fix. Default mode is SOFT (collect
# & report every failing block); pass --strict to hard-stop on the first fail
# (for later reuse as a guard inside a fitting driver).
#
#   $R = "C:\Program Files\R\R-4.5.1\bin\Rscript.exe"
#   & $R fit_diagnostics.R --setting=4 --datasets=1:5
#########################################################################

# NOTE: ar2_load.R starts with rm(list = ls()), which wipes everything defined
# before it. So we source FIRST, then define functions and parse args. (This is
# the same trap run-settings.R documents.)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1]),
    winslash = "/", mustWork = FALSE))
  if (dir.exists(script_dir)) setwd(script_dir)
}
source("../../simstudy/ar2_load.R", chdir = TRUE)
source("../simstudy/time_block_helpers.R")
load("../simstudy/simdata.RData", envir = .GlobalEnv)
options(warn = 1)

# ---- the assertions ----------------------------------------------------
# check_fit_consistency() lives in fit_assertions.R alongside this
# harness. Both moved here from ../simstudy on 2026-08-05: the simstudy
# arm records descriptive diagnostics only and does not gate on A/B/C, so
# the assertions belong with the study that does.
source("./fit_assertions.R")

# ---- CLI ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
getflag <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
setting <- as.integer(getflag("setting", "4"))
method_id <- as.integer(getflag("method", "2"))       # 2 = AR(2)
strict <- any(args == "--strict")

datasets <- eval(parse(text = getflag("datasets", "1:5")))
iters <- 20000; burn <- 10000; update <- 2000
K <- 1L
blocks <- tbf_blocks(block_seams, block_H, nt)
nblk <- length(blocks)
marginal_sd <- sd(y[, , , setting])

# ground truth (simulation) -- phi shared across channels; see setup.R
phi_by_setting <- list(`1` = c(0, 0), `2` = c(0.12, -0.05),
                       `3` = c(0.60, -0.30), `4` = c(0.80, -0.35))
truth <- list(lambda = 3, beta0 = 10,
              phi = phi_by_setting[[as.character(setting)]])

cat(sprintf("fit_diagnostics: setting %d (phi truth = (%.2f, %.2f)), method %d\n",
            setting, truth$phi[1], truth$phi[2], method_id))
cat(sprintf("  datasets = %s | marginal sd = %.2f | spread threshold = %.1f\n\n",
            paste(datasets, collapse = ","), marginal_sd, 3 * marginal_sd))

# ---- refit (faithful RNG: seed once, walk blocks in order) -------------
rows <- list()
for (d in datasets) {
  y.d <- y[, , d, setting]
  set.seed(get_tbf_seed(setting, method_id, d))     # SAME seed as the pipeline
  for (b in seq_len(nblk)) {
    blk <- blocks[[b]]
    y.train <- y.d[, blk$train_times, drop = FALSE]
    x.train <- x[, blk$train_times, , drop = FALSE]
    x.block <- x[, blk$test_times, , drop = FALSE]

    fit <- mcmc(
      y = y.train, x = x.train, s = s,
      method = "t", skew = TRUE, thresh.all = 0, thresh.quant = TRUE,
      nknots = 1, iterplot = FALSE, iters = iters, burn = burn, update = update,
      min.s = c(0, 0), max.s = c(10, 10),
      temporalw = TRUE, temporaltau = TRUE, temporalz = TRUE,
      ar2_w = TRUE, ar2_tau = TRUE, ar2_z = TRUE,
      rho.upper = 15, nu.upper = 10              # NO z.init: use backend default
    )
    yhat <- forecast_block(fit = fit, seam = blk$seam, H = block_H,
                           x_block = x.block, s = s, ar2 = TRUE)

    chk <- check_fit_consistency(fit, yhat, truth, marginal_sd,
                                 data_mean = mean(y.train))
    chk <- cbind(dataset = d, block = b, seam = blk$seam, chk)
    rows[[length(rows) + 1]] <- chk

    fails <- c(if (isFALSE(chk$A_zconsist)) "A", if (isFALSE(chk$B_truth)) "B",
               if (isFALSE(chk$C_spread)) "C")
    tag <- if (length(fails)) paste0("  <<< FAIL ", paste(fails, collapse = ""),
                                     if (!chk$C_spread) sprintf(" (SD %.1f)", chk$sd_lead_max) else "") else "  ok"
    cat(sprintf("d%d blk%d (seam %3d): lambda %+8.2f  zbar %.3f (want %.3f)  phi1.z %+.2f  SDmax %6.2f%s\n",
                d, b, blk$seam, chk$lambda, chk$zbar, chk$z_pred, chk$phi1.z,
                chk$sd_lead_max, tag))
    if (strict && length(fails)) stop(sprintf("strict: d%d blk%d failed %s", d, b, paste(fails, collapse = "")))
    rm(fit, yhat); gc()
  }
}

diag <- do.call(rbind, rows)
dir.create("output/diag", recursive = TRUE, showWarnings = FALSE)
save(diag, truth, marginal_sd, file = sprintf("output/diag/fit_diag_%d.RData", setting))

# ---- summary table by block --------------------------------------------
cat("\n================= PER-BLOCK SUMMARY =================\n")
cat(sprintf("%-6s %-6s | %8s %8s %8s | %8s %8s %8s\n",
            "block", "seam", "lambda", "zbar/want", "phi1.z", "SDmax", "%failC", "%failA"))
for (b in seq_len(nblk)) {
  sub <- diag[diag$block == b, ]
  if (!nrow(sub)) next
  cat(sprintf("%-6d %-6d | %+8.2f %8s %+8.2f | %8.2f %7.0f%% %7.0f%%\n",
              b, sub$seam[1], mean(sub$lambda),
              sprintf("%.2f/%.2f", mean(sub$zbar), mean(sub$z_pred)),
              mean(sub$phi1.z), mean(sub$sd_lead_max),
              100 * mean(!sub$C_spread), 100 * mean(!sub$A_zconsist)))
}
cat(sprintf("\nsaved -> output/diag/fit_diag_%d.RData\n", setting))
cat("Assertion A=z-self-consistency  B=truth-recovery  C=predictive-spread\n")
