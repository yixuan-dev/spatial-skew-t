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
source("../../../R/mcmc_cont_lambda.R", chdir = T)
source("../../../R/auxfunctions.R")
source("../max-stab/MCMC4MaxStable.R", chdir = T)

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
