# Fit-health gate for us-all-full-204.R, run BEFORE predict-cincy-204.R.
#
# The three assertions come from the two post-mortems of the z.init bug:
#   1. Self-consistency (tex/z_init_bug "What should have caught it"): the
#      model asserts z ~ HN(0, tau^{-1/2}), so E[z] = tau^{-1/2} sqrt(2/pi).
#      A fitted mean(z) far below that means z is still frozen.
#   2. No reflected ridge (tex/lambda_phiz_ridge): the (beta0, lambda, zbar)
#      ridge has a reflected branch where lambda flips sign with inflated
#      magnitude and beta0 compensates; it corrupts every proper score and
#      appears stochastically in ~1/4 of chains, so it must be checked per
#      fit, not assumed away.
#   3. Basic chain sanity: no NaN in the kept draws, phi's inside (-1, 1).
#
# Usage:  Rscript check-full-204.R [dev]
#   "dev" checks results/us-all-full-204-dev.RData, else the prod file.
# Exits non-zero on failure so it can gate a pipeline.
rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)
suffix <- if (length(args) > 0 && tolower(args[1]) == "dev") "-dev" else ""
infile <- paste0("results/us-all-full-204", suffix, ".RData")
cat("checking:", infile, "\n")
load(infile)

fails <- character(0)
note  <- function(ok, msg) {
  cat(if (ok) "  PASS  " else "  FAIL  ", msg, "\n", sep = "")
  if (!ok) fails <<- c(fails, msg)
}

## ---- 1. z self-consistency -------------------------------------------------
# nknots = 1: fit$z and fit$tau are (kept iters) x nt matrices.
z   <- fit$z
tau <- fit$tau
stopifnot(!is.null(z), !is.null(tau))
sig  <- 1 / sqrt(tau)                     # per-draw, per-day HN scale
ez   <- mean(sig) * sqrt(2 / pi)          # model-implied E[z]
zbar <- mean(z)
ratio <- zbar / ez
cat(sprintf("  z: mean=%.3f  model E[z]=%.3f  ratio=%.3f  (frozen fit ~ 0.1)\n",
            zbar, ez, ratio))
note(!all(z == 0), "z is not identically zero")
note(ratio > 0.5 && ratio < 2.0,
     sprintf("mean(z)/E[z] = %.3f within [0.5, 2] (z unfrozen, self-consistent)", ratio))

## ---- 2. lambda not on the reflected branch ---------------------------------
lam <- fit$lambda
p.neg <- mean(lam < 0)
lam.m <- mean(lam)
cat(sprintf("  lambda: mean=%.3f  sd=%.3f  P(<0)=%.2f  range=[%.2f, %.2f]\n",
            lam.m, sd(lam), p.neg, min(lam), max(lam)))
# Reflection signature: |lambda| inflated (|.| >> healthy 0-5 scale) with a
# committed sign. A posterior straddling 0 is fine (ozone residuals are
# near-symmetric); a posterior locked at lambda < -4 is the ridge.
note(!(abs(lam.m) > 4 && (p.neg > 0.95 || p.neg < 0.05) && lam.m < 0),
     "lambda not committed to the reflected (negative, inflated) branch")
note(sd(lam) < 10,
     sprintf("lambda sd = %.2f < 10 (a prior-random-walk sd is ~20)", sd(lam)))

## ---- 3. chain sanity -------------------------------------------------------
note(!anyNA(fit$beta) && !anyNA(z) && !anyNA(tau), "no NaN/NA in kept draws")
for (nm in c("phi.z", "phi.w", "phi.tau")) {
  ph <- fit[[nm]]
  if (!is.null(ph)) {
    note(all(is.finite(ph)) && all(abs(ph) < 1),
         sprintf("%s finite and inside (-1, 1); mean=%.3f sd=%.3f", nm, mean(ph), sd(ph)))
    if (nm == "phi.w" && runtime_info$nknots == 1L) {
      # K = 1: the partition update is guarded by nknots > 1, so phi.w never
      # moves -- structurally trivial, not a pathology (see convergence notes).
      cat("  SKIP  phi.w mixing check (nknots = 1: partition update inactive)\n")
    } else {
      note(sd(ph) > 0, sprintf("%s is mixing (sd > 0)", nm))
    }
  }
}
b0 <- mean(fit$beta[, 1])
cat(sprintf("  beta0 mean=%.2f (a ridge-compensating beta0 sits far above the data mean)\n", b0))

## ---- verdict ---------------------------------------------------------------
cat("\n")
if (length(fails) > 0) {
  cat("HEALTH CHECK FAILED (", length(fails), " assertion(s) )\n", sep = "")
  quit(save = "no", status = 1)
}
cat("HEALTH CHECK PASSED -- OK to run predict-cincy-204.R\n")
