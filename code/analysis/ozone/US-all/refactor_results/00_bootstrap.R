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
  target_dir <- script_dir
  if (identical(basename(script_dir), "refactor_results")) {
    target_dir <- dirname(script_dir)
  }

  if (dir.exists(target_dir)) {
    setwd(target_dir)
  }

  invisible(target_dir)
}

safe_as_integer <- function(x) {
  suppressWarnings(as.integer(x))
}

default_probs <- function(include_999 = FALSE) {
  probs <- default_threshold_probs()
  if (isTRUE(include_999)) {
    probs <- c(probs, 0.999)
  }
  probs
}

default_threshold_probs <- function() {
  c(seq(0, 0.9, by = 0.1), 0.95, 0.98, 0.99, 0.995)
}

us_all_output_root <- function() {
  cwd_name <- basename(normalizePath(getwd(), winslash = "/", mustWork = FALSE))
  base_dir <- if (identical(cwd_name, "refactor_results")) ".." else "."
  file.path(base_dir, "output", "us-all")
}

ensure_us_all_output_dirs <- function() {
  root <- us_all_output_root()
  dirs <- c(
    root,
    file.path(root, "results"),
    file.path(root, "tables"),
    file.path(root, "plots"),
    file.path(root, "logs")
  )
  for (d in dirs) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(dirs)
}

us_all_output_path <- function(..., subdir = c("results", "tables", "plots", "logs")) {
  subdir <- match.arg(subdir)
  file.path(us_all_output_root(), subdir, ...)
}

load_us_all_context <- function(
    setup_file = "us-all-setup.RData",
    settings_file = "settings.csv",
    aux_file = "../../../R/auxfunctions.R") {
  resolve_existing_path <- function(candidates, label) {
    candidates <- unique(candidates[nzchar(candidates)])
    hit <- candidates[file.exists(candidates)]
    if (length(hit) == 0) {
      stop(label, " not found. Checked: ", paste(candidates, collapse = ", "))
    }
    hit[1]
  }

  setup_path <- resolve_existing_path(
    c(setup_file, file.path("..", setup_file)),
    "Setup file"
  )

  settings_path <- resolve_existing_path(
    c(settings_file, file.path("..", settings_file)),
    "Settings file"
  )

  aux_path <- resolve_existing_path(
    c(
      aux_file,
      file.path("..", aux_file),
      "../../../R/auxfunctions.R",
      "../../../../R/auxfunctions.R"
    ),
    "Auxiliary functions file"
  )

  load_env <- new.env(parent = globalenv())
  load(setup_path, envir = load_env)

  required_objects <- c("Y", "cv.lst")
  missing_objects <- required_objects[!vapply(
    required_objects,
    function(nm) exists(nm, envir = load_env, inherits = FALSE),
    logical(1)
  )]
  if (length(missing_objects) > 0) {
    stop(
      "Missing required objects in ",
      setup_path,
      ": ",
      paste(missing_objects, collapse = ", ")
    )
  }

  source(aux_path)
  settings <- read.csv(settings_path, stringsAsFactors = FALSE)
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
