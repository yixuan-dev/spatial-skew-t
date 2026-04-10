# Shared bootstrap utilities for US-all ozone result pipelines.

set_working_dir_to_script <- function() {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_arg) == 0) {
    return(invisible(NULL))
  }

  script_path <- normalizePath(
    sub("^--file=", "", script_arg[1]),
    winslash = "/",
    mustWork = FALSE
  )
  script_dir <- dirname(script_path)
  if (dir.exists(script_dir)) {
    setwd(script_dir)
  }

  invisible(script_dir)
}

safe_as_integer <- function(x) {
  suppressWarnings(as.integer(x))
}

default_probs <- function(include_999 = FALSE) {
  probs <- c(0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995)
  if (isTRUE(include_999)) {
    probs <- c(probs, 0.999)
  }
  probs
}

load_us_all_context <- function(
    setup_file = "us-all-setup.RData",
    settings_file = "settings.csv",
    aux_file = "../../../R/auxfunctions.R"
) {
  if (!file.exists(setup_file)) {
    stop("Setup file not found: ", setup_file)
  }
  if (!file.exists(settings_file)) {
    stop("Settings file not found: ", settings_file)
  }
  if (!file.exists(aux_file)) {
    stop("Auxiliary functions file not found: ", aux_file)
  }

  load_env <- new.env(parent = globalenv())
  load(setup_file, envir = load_env)

  required_objects <- c("Y", "cv.lst")
  missing_objects <- required_objects[!vapply(
    required_objects,
    function(nm) exists(nm, envir = load_env, inherits = FALSE),
    logical(1)
  )]
  if (length(missing_objects) > 0) {
    stop(
      "Missing required objects in ",
      setup_file,
      ": ",
      paste(missing_objects, collapse = ", ")
    )
  }

  source(aux_file)
  settings <- read.csv(settings_file, stringsAsFactors = FALSE)
  if ("setting" %in% names(settings)) {
    settings$setting_num <- safe_as_integer(settings$setting)
  }

  list(
    Y = get("Y", envir = load_env),
    cv_lst = get("cv.lst", envir = load_env),
    settings = settings,
    setup_env = load_env
  )
}
