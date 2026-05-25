# Loads the canonical AR(2) infrastructure from ../ar2/.
# Edit AR(2) code in code/R/ar2/ only; this loader picks it up.

local({
  this_file <- sys.frame(1)$ofile
  if (is.null(this_file)) {
    this_file <- normalizePath("load_ar2.R", mustWork = FALSE)
  }
  ar2_dir <- normalizePath(file.path(dirname(this_file), "..", "ar2"),
                           mustWork = TRUE)

  # Note: mcmc_ar2.R is intentionally NOT sourced here — it only defines
  # mcmc(), which mcmc_prop.R overwrites with its own definition.
  source(file.path(ar2_dir, "auxfunctions.R"))
  source(file.path(ar2_dir, "auxfunctions_ar2.R"))
  source(file.path(ar2_dir, "update_params.R"))
  source(file.path(ar2_dir, "update_params_ar2.R"))
  source(file.path(ar2_dir, "update_params_cpp.R"))
})
