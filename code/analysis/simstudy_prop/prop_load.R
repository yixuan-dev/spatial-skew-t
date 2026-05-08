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

# Resolve an arbitrary dataset path. If data_arg is NULL/empty, use the
# default simdata.RData lookup. Otherwise try the given path as-is and
# fall back to ../simstudy/<basename> (so --data=simdata_def.RData works
# when the file lives only in the sibling simstudy/ folder).
resolve_simstudy_prop_data_path <- function(data_arg = NULL) {
  if (is.null(data_arg) || !nzchar(data_arg)) {
    return(resolve_simstudy_prop_data())
  }
  if (file.exists(data_arg)) {
    return(normalizePath(data_arg, winslash = "/", mustWork = TRUE))
  }
  sibling <- file.path("..", "simstudy", basename(data_arg))
  if (file.exists(sibling)) {
    return(normalizePath(sibling, winslash = "/", mustWork = TRUE))
  }
  stop(
    sprintf(
      "Data file not found: %s (also tried %s)",
      data_arg,
      sibling
    ),
    call. = FALSE
  )
}

load(file = resolve_simstudy_prop_data())

source_required("../../R/prop/auxfunctions.R")
source_required("../../R/prop/update_params_cpp.R")
source_required("../../R/prop/update_params.R")
source_required("../../R/prop/prop_utils.R")
source_required("../../R/prop/prop_basis.R")
source_required("../../R/prop/prop_covariance.R")
source_required("../../R/prop/prop_imputation.R")
source_required("../../R/prop/prop_modules.R")
source_required("../../R/prop/mcmc_prop.R")

options(warn = 2)
