library(fields)
library(emulator)

source_required <- function(path, chdir = FALSE) {
  if (!file.exists(path)) {
    stop(sprintf("Required source file not found: %s", path), call. = FALSE)
  }
  source(path, chdir = chdir)
}

load(file = "ozone-prop-setup.RData")

source_required("../../R/prop/auxfunctions.R")
source_required("../../R/prop/update_params_cpp.R")
source_required("../../R/prop/update_params.R")
source_required("../../R/prop/prop_utils.R")
source_required("../../R/prop/prop_basis.R")
source_required("../../R/prop/prop_covariance.R")
source_required("../../R/prop/prop_imputation.R")
source_required("../../R/prop/prop_modules.R")
source_required("../../R/prop/mcmc_prop.R")

set_blas_threads <- function(n = 1) {
  if (exists("openblas.set.num.threads", mode = "function")) {
    get("openblas.set.num.threads", mode = "function")(n)
  }
}
set_blas_threads(1)
