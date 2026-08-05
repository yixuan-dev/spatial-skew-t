#########################################################################
# inventory_tbf.R -- fit inventory + missing-cell call generator for the
# time-block Experiment A arm (5 expanding-window blocks per cell).
#
# The time-block analogue of code/analysis/simstudy/DESKTOP-61SBCCI/
# inventory_gen.R: no MRTS K axis, and the validity test is about the
# 5-block structure and the lambda prior rather than about fit.1.
#
# Usage (from code/analysis/time_block_forecast/simstudy):
#   Rscript DESKTOP-61SBCCI/inventory_tbf.R <setting> <datasets> <methods> <workers> <calls_out>
#   e.g. Rscript DESKTOP-61SBCCI/inventory_tbf.R 5 1:5 "c(1,2,4)" 12 DESKTOP-61SBCCI/calls_s5_d1-5.txt
#
# A fit counts as PRESENT only if it is plausibly complete AND was produced
# under the HN prior. The prior check is the important one: a fit left over
# from an earlier non-HN run is otherwise indistinguishable on disk, and
# reusing it would silently mix two priors inside one score cache.
#
# stdout ends with two machine-readable lines:
#   MISSING <n>
#   CALLS <n>
#########################################################################

args <- commandArgs(trailingOnly = TRUE)

# Which lambda-prior arm to inventory. Selects results_<arm>/ and decides
# what counts as a valid fit; defaults to hn, the arm every campaign so
# far has used.
prior_arm <- "hn"
prior_flag <- grep("^--prior=", args, value = TRUE)
if (length(prior_flag)) {
  prior_arm <- tolower(sub("^--prior=", "", prior_flag[1]))
  if (!prior_arm %in% c("hn", "n")) {
    stop(sprintf("--prior must be hn or n, got '%s'", prior_arm), call. = FALSE)
  }
  args <- args[!grepl("^--prior=", args)]
}
hn_want <- identical(prior_arm, "hn")

if (length(args) != 5) {
  stop(paste("usage: inventory_tbf.R [--prior=hn|n] <setting>",
             "<datasets e.g. 1:5> <methods e.g. c(1,2,4)> <workers>",
             "<calls_out>"), call. = FALSE)
}
setting <- as.integer(args[1])
datasets <- eval(parse(text = args[2]))
methods <- eval(parse(text = args[3]))
workers <- as.integer(args[4])
calls_out <- args[5]
stopifnot(setting >= 1L, all(datasets >= 1), all(methods >= 1), workers >= 1)

# run from code/analysis/time_block_forecast/simstudy (the driver guarantees
# this; be robust anyway)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0L) {
  sd_ <- dirname(normalizePath(sub("^--file=", "", script_arg[1]),
                               winslash = "/", mustWork = FALSE))
  if (dir.exists(sd_)) setwd(file.path(sd_, ".."))
}

N_BLOCKS <- 5L
MIN_BYTES <- 700e6   # a complete 5-block cell is ~810 MB

fname <- function(m, d) {
  sprintf("results_%s/%d-%d-%d.RData", prior_arm, setting, m, d)
}

grid <- expand.grid(m = methods, d = datasets)
grid$file <- mapply(fname, grid$m, grid$d)
grid$exists <- file.exists(grid$file)

# A cell is written once, at the end of all 5 blocks, so a short file means
# the process was killed mid-save. Size is a cheap screen for every file;
# the newest few get a full structural load as well.
short <- grid$exists & file.size(grid$file) < MIN_BYTES
if (any(short)) {
  cat("TRUNCATED (treated as missing, will be re-run):\n")
  print(grid$file[short])
  grid$exists[short] <- FALSE
}

ex <- grid[grid$exists, ]
if (nrow(ex) > 0) {
  mt <- file.mtime(ex$file)
  newest <- ex$file[order(mt, decreasing = TRUE)][seq_len(min(3, nrow(ex)))]
  bad <- character(0)
  for (f in newest) {
    ok <- tryCatch({
      e <- new.env()
      load(f, envir = e)
      fc <- e$forecast
      !is.null(fc) &&
        length(fc$blocks) == N_BLOCKS &&
        all(vapply(fc$blocks, function(b) {
          !is.null(b$yhat) && length(dim(b$yhat)) == 3L && !anyNA(b$yhat[1, , ])
        }, TRUE)) &&
        !is.null(e$fit_diag) && nrow(e$fit_diag) == N_BLOCKS &&
        identical(isTRUE(e$runtime_info$control$lambda_positive), hn_want)
    }, error = function(err) FALSE)
    if (!ok) bad <- c(bad, f)
  }
  if (length(bad)) {
    cat(sprintf("INVALID (bad structure or not a %s fit; will be re-run):\n",
                if (hn_want) "--hn" else "--prior=n"))
    print(bad)
    grid$exists[grid$file %in% bad] <- FALSE
  } else {
    cat(sprintf("newest %d files all load cleanly under the %s prior\n",
                length(newest), if (hn_want) "HN" else "N"))
  }
}

miss <- grid[!grid$exists, c("m", "d")]
cat(sprintf("setting %d, datasets %s, methods %s: existing valid %d / %d, missing %d\n",
            setting, args[2], args[3], sum(grid$exists), nrow(grid), nrow(miss)))

spec <- function(v) {                       # compact R index spec for the CLI
  v <- sort(unique(v))
  if (length(v) == 1) return(as.character(v))
  if (all(diff(v) == 1)) return(sprintf("%d:%d", v[1], v[length(v)]))
  sprintf("c(%s)", paste(v, collapse = ","))
}

# Group methods whose missing dataset-sets are identical, so a fresh chunk
# collapses to a single call and therefore a single PSOCK cluster.
lines <- character(0)
if (nrow(miss) > 0) {
  key <- vapply(split(miss$d, miss$m), function(x) paste(sort(x), collapse = ","), "")
  for (dsig in unique(key)) {
    ms <- as.integer(names(key)[key == dsig])
    ds <- as.integer(strsplit(dsig, ",")[[1]])
    lines <- c(lines, sprintf(
      'Rscript run-settings.R %s --setting=%d "%s" %d "%s"',
      if (hn_want) "--hn" else "--prior=n",
      setting, spec(ds), workers, spec(ms)))
  }
}

writeLines(lines, calls_out)
cat(sprintf("MISSING %d\n", nrow(miss)))
cat(sprintf("CALLS %d\n", length(lines)))
