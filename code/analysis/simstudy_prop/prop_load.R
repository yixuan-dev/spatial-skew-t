library(fields)
library(emulator)

source_required <- function(path, chdir = FALSE) {
  if (!file.exists(path)) {
    stop(sprintf("Required source file not found: %s", path), call. = FALSE)
  }
  source(path, chdir = chdir)
}

resolve_simstudy_prop_data <- function() {
  candidates <- c(
    "./simdata.RData",
    "../simstudy/simdata.RData"
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0L) {
    stop(
      "Unable to locate simdata.RData. Checked ./simdata.RData and ../simstudy/simdata.RData",
      call. = FALSE
    )
  }
  normalizePath(existing[1], winslash = "/", mustWork = TRUE)
}

load(file = resolve_simstudy_prop_data())

source_required("../../R/prop/auxfunctions.R")
source_required("../../R/prop/update_params_cpp.R")
source_required("../../R/prop/update_params.R")
source_required("../../R/prop/prop_utils.R")
source_required("../../R/prop/prop_basis.R")
source_required("../../R/prop/prop_covariance.R")
source_required("../../R/prop/prop_imputation.R")
source_required("../../R/prop/mcmc_prop.R")

options(warn = 2)
