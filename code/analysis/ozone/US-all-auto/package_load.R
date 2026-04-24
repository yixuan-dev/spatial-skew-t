# package_load.R
#
# Loads packages and MCMC source files for the US-all-auto pipeline.
# Must be sourced AFTER .us_all_auto_root and backend are defined.
#
# backend == "ar2"    -> AR2 model functions (mcmc_ar2)
# backend == "legacy" -> legacy model functions (mcmc / mcmc_cont_lambda)

if (!exists("backend")) {
  stop("'backend' must be set before sourcing package_load.R", call. = FALSE)
}
if (!exists(".us_all_auto_root")) {
  stop("'.us_all_auto_root' must be set before sourcing package_load.R", call. = FALSE)
}

library(fields)
library(SpatialTools)
library(mvtnorm)
library(Rcpp)
library(compiler)
library(autoFRK)
enableJIT(3)

.r_root <- normalizePath(file.path(.us_all_auto_root, "../../../R"), winslash = "/")
.ozone_dir <- normalizePath(file.path(.us_all_auto_root, ".."), winslash = "/")

if (backend == "ar2") {
  source(file.path(.r_root, "ar2", "auxfunctions_ar2.R"))
  source(file.path(.r_root, "ar2", "auxfunctions.R"))
  source(file.path(.r_root, "ar2", "update_params_ar2.R"))
  source(file.path(.r_root, "ar2", "update_params.R"))
  source(file.path(.r_root, "ar2", "update_params_cpp.R"))
  source(file.path(.r_root, "ar2", "mcmc_ar2.R"))
} else {
  source(file.path(.r_root, "mcmc_cont_lambda.R"), chdir = TRUE)
  source(file.path(.r_root, "auxfunctions.R"))
  source(file.path(.ozone_dir, "max-stab", "MCMC4MaxStable.R"), chdir = TRUE)
}

set_blas_threads <- function(n = 1) {
  if (exists("openblas.set.num.threads", mode = "function")) {
    get("openblas.set.num.threads", mode = "function")(n)
  }
}
set_blas_threads(1)
