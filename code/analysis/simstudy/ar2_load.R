rm(list = ls())

library(fields)
library(SpatialTools)
library(Rcpp)

#### Load simdata
load(file='./simdata.RData')

source_required <- function(path, chdir = FALSE) {
  if (!file.exists(path)) {
    stop(sprintf("Required source file not found: %s", path), call. = FALSE)
  }
  source(path, chdir = chdir)
}

source_required("../../R/ar2/auxfunctions_ar2.R")
source_required("../../R/ar2/auxfunctions.R")
source_required("../../R/ar2/update_params_ar2.R")
source_required("../../R/ar2/update_params.R")
source_required("../../R/ar2/update_params_cpp.R")
source_required("../../R/ar2/mcmc_ar2.R")
source_required('./max-stab/Bayes_GEV.R')
source_required('./max-stab/MCMC4MaxStable.R', chdir = TRUE)

options(warn=2)
