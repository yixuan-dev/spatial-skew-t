rm(list = ls())

library(fields)
library(SpatialTools)
library(mvtnorm)
library(Rcpp)
library(compiler)
enableJIT(3)

load("us-all-setup.RData")
# Optional private utility file used on the original author's machine.
# Keep experiments portable by sourcing only when available.
useful_r_candidates <- c(
  "~/repos-git/usefulR/usefulfunctions.R",
  "../usefulR/usefulfunctions.R",
  "./usefulfunctions.R"
)
useful_r_found <- useful_r_candidates[file.exists(path.expand(useful_r_candidates))]
if (length(useful_r_found) > 0) {
  source(useful_r_found[1], chdir = TRUE)
}

source_required <- function(path, chdir = FALSE) {
  if (!file.exists(path)) {
    stop(sprintf("Required source file not found: %s", path), call. = FALSE)
  }
  source(path, chdir = chdir)
}

source_required("../../../R/ar2/auxfunctions_ar2.R")
source_required("../../../R/ar2/auxfunctions.R")
source_required("../../../R/ar2/update_params_ar2.R")
source_required("../../../R/ar2/update_params.R")
source_required("../../../R/ar2/update_params_cpp.R")
source_required("../../../R/ar2/mcmc_ar2.R")

set_blas_threads <- function(n = 1) {
  if (exists("openblas.set.num.threads", mode = "function")) {
    get("openblas.set.num.threads", mode = "function")(n)
  }
}

if (Sys.info()["nodename"] == "cwl-mth-sam-001") {
  # setMKLthreads(1)
  set_blas_threads(1)
  do.upload <- TRUE
} else if (Sys.info()["sysname"] == "Darwin") {
  do.upload <- TRUE
} else {
  do.upload <- FALSE
  # set number of threads to use
  set_blas_threads(1)
}
