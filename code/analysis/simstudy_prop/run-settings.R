#########################################################################
# simstudy_prop naming-parallel wrapper
# Mirrors code/analysis/simstudy/run-settings.R naming, but delegates to
# the prop engine in run-prop.R.
#########################################################################

script_dir_override <- trimws(Sys.getenv("SIMSTUDY_PROP_SCRIPT_DIR", unset = ""))
if (nzchar(script_dir_override) && dir.exists(script_dir_override)) {
  setwd(normalizePath(script_dir_override, winslash = "/", mustWork = TRUE))
} else {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_arg) > 0) {
    script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = FALSE)
    script_dir <- dirname(script_path)
    if (dir.exists(script_dir)) {
      setwd(script_dir)
    }
  }
}

source("./run-prop.R", chdir = TRUE)
