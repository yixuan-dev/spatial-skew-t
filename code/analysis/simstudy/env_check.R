#########################################################################
# env_check.R -- one-shot environment check + auto-install for the
#                simstudy fit/score pipeline (run-settings.R + scores.R).
#
# Run from anywhere:
#   Rscript env_check.R          # check, auto-install missing, smoke test
#   Rscript env_check.R --fast   # skip the full ar2 pipeline load (~1 min)
#
# Exit 0 = environment ready.  Non-zero = something is missing; the
# message says exactly what and how to fix it.
#
# What the pipeline actually needs (grep'd from the sources):
#   ar2_load.R              -> fields, SpatialTools, Rcpp
#   update_params_cpp.R     -> inline + RcppArmadillo (runtime C++ compile!)
#   helpers.R (MRTS basis)  -> autoFRK
#   scores.R                -> scoringRules
# plus a working C++ toolchain (Rtools on Windows) because the sampler
# compiles its C++ kernels at load time via inline::cxxfunction.
#########################################################################

fast <- "--fast" %in% commandArgs(trailingOnly = TRUE)
repo <- "https://cloud.r-project.org"

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0L) {
  sd <- dirname(normalizePath(sub("^--file=", "", script_arg[1]),
                              winslash = "/", mustWork = FALSE))
  if (dir.exists(sd)) setwd(sd)          # -> code/analysis/simstudy
}

fail <- function(...) { cat("\n[ENV CHECK FAILED] ", sprintf(...), "\n", sep = ""); quit(status = 1) }
step <- function(...) cat("==", sprintf(...), "\n")

## 1. R version --------------------------------------------------------
step("R %s on %s", getRversion(), R.version$platform)
if (getRversion() < "4.2.0") fail("R >= 4.2 required (found %s).", getRversion())

## 2. packages: check + auto-install ----------------------------------
required <- c("fields", "SpatialTools", "Rcpp", "RcppArmadillo",
              "inline", "autoFRK", "scoringRules")
missing <- required[!vapply(required, requireNamespace, TRUE, quietly = TRUE)]
if (length(missing)) {
  step("installing missing packages: %s", paste(missing, collapse = ", "))
  install.packages(missing, repos = repo)
  still <- missing[!vapply(missing, requireNamespace, TRUE, quietly = TRUE)]
  if (length(still)) fail("could not install: %s", paste(still, collapse = ", "))
} else {
  step("all %d required packages present", length(required))
}

## 3. C++ toolchain (Rtools on Windows) --------------------------------
if (.Platform$OS.type == "windows" && !nzchar(Sys.which("make"))) {
  fail(paste0("'make' not found -- install Rtools for R %s.%s from\n",
              "  https://cran.r-project.org/bin/windows/Rtools/\n",
              "then open a NEW terminal and re-run this check."),
       R.version$major, strsplit(R.version$minor, "[.]")[[1]][1])
}
step("build tool: make = %s", Sys.which("make"))

## 4. compile smoke test (same plugin the sampler uses) ----------------
step("compiling a test C++ kernel via inline/RcppArmadillo ...")
suppressMessages(library(inline))
test_fn <- tryCatch(
  cxxfunction(signature(x = "numeric"),
              body = "arma::vec v = Rcpp::as<arma::vec>(x); return Rcpp::wrap(arma::accu(v));",
              plugin = "RcppArmadillo"),
  error = function(e) fail("C++ compile failed: %s", conditionMessage(e)))
if (abs(test_fn(c(1, 2, 3)) - 6) > 1e-12) fail("compiled kernel returned a wrong result.")
step("compile smoke test passed")

## 5. full pipeline load (compiles all sampler kernels) ----------------
## NOTE: ar2_load.R starts with rm(list = ls()), which wipes this script's
## helper functions from the global env -- so past this point use plain
## cat()/quit() only, never step()/fail().
if (!fast) {
  step("loading the full ar2 pipeline (ar2_load.R; compiles all kernels) ...")
  ok <- tryCatch({ source("./ar2_load.R", chdir = TRUE); TRUE },
                 error = function(e) { cat("  error: ", conditionMessage(e), "\n"); FALSE })
  if (!ok || !exists("mcmc")) {
    cat("\n[ENV CHECK FAILED] ar2_load.R did not load cleanly (see error above).\n")
    quit(status = 1)
  }
  cat("== ar2 pipeline loads: mcmc() available\n")
}

cat("\n[ENV OK] ready for run-settings.R / scores.R\n")
