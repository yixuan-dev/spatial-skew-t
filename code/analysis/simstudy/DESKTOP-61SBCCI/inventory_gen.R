#########################################################################
# inventory_gen.R -- generalized fit inventory + missing-cell call
# generator for the nonsta n=30 repeats run (settings 16 / 20).
#
# Generalizes output/s20_resume/{s20_inventory.R, s20_resume_gen.R}:
# parameterized setting, dataset range, and worker count.
#
# Usage (from code/analysis/simstudy):
#   Rscript DESKTOP-61SBCCI/inventory_gen.R <setting> <datasets> <workers> <calls_out>
#   e.g. Rscript DESKTOP-61SBCCI/inventory_gen.R 20 11:15 16 DESKTOP-61SBCCI/calls_s20_d11-15.txt
#
# Scans results_nonsta/ for the datasets x methods(1:5) x K grid, validates
# the newest 12 existing files (mid-write corruption after an interrupt),
# and writes one run-settings.R call per line to <calls_out>.
# K=0 baseline cells get their own call WITHOUT the 5th CLI arg, because
# build_run_plan drops the baseline whenever a K list is passed.
#
# stdout ends with two machine-readable lines:
#   MISSING <n>
#   CALLS <n>
#########################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("usage: inventory_gen.R <setting> <datasets e.g. 11:15> <workers> <calls_out>",
       call. = FALSE)
}
setting  <- as.integer(args[1])
datasets <- eval(parse(text = args[2]))
workers  <- as.integer(args[3])
calls_out <- args[4]
stopifnot(setting %in% c(16L, 20L) || setting %in% 1:21,
          all(datasets >= 1), workers >= 1)

# run from code/analysis/simstudy (driver guarantees this; be robust anyway)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0L) {
  sd <- dirname(normalizePath(sub("^--file=", "", script_arg[1]),
                              winslash = "/", mustWork = FALSE))
  if (dir.exists(sd)) setwd(file.path(sd, ".."))
}

methods <- 1:5
Ks <- c(0, 3, 5, 8, 10, 12, 15, 20, 25, 30, 40, 50)

fname <- function(m, d, K) {
  if (K == 0) sprintf("results_nonsta/%d-%d-%d.RData", setting, m, d)
  else        sprintf("results_nonsta/%d-%d-%d-K%d.RData", setting, m, d, K)
}

grid <- expand.grid(m = methods, d = datasets, K = Ks)
grid$file <- mapply(fname, grid$m, grid$d, grid$K)
grid$exists <- file.exists(grid$file)

# validate the newest existing files (workers may have been killed mid-save)
ex <- grid[grid$exists, ]
if (nrow(ex) > 0) {
  mt <- file.mtime(ex$file)
  newest <- ex$file[order(mt, decreasing = TRUE)][seq_len(min(12, nrow(ex)))]
  bad <- character(0)
  for (f in newest) {
    ok <- tryCatch({ e <- new.env(); load(f, envir = e); !is.null(e$fit.1) },
                   error = function(err) FALSE)
    if (!ok) bad <- c(bad, f)
  }
  if (length(bad)) {
    cat("CORRUPT (treated as missing, will be re-run):\n")
    print(bad)
    grid$exists[grid$file %in% bad] <- FALSE
  } else {
    cat(sprintf("newest %d files all load cleanly\n", length(newest)))
  }
}

miss <- grid[!grid$exists, c("m", "d", "K")]
cat(sprintf("setting %d, datasets %s: existing valid %d / %d, missing %d\n",
            setting, args[2], sum(grid$exists), nrow(grid), nrow(miss)))

spec <- function(v) {                       # compact R index spec for the CLI
  v <- sort(unique(v))
  if (length(v) == 1) return(as.character(v))
  if (all(diff(v) == 1)) return(sprintf("%d:%d", v[1], v[length(v)]))
  sprintf("c(%s)", paste(v, collapse = ","))
}

lines <- character(0)

emit_group <- function(sub, k_spec) {
  # group methods by identical missing dataset-sets; one call per group
  out <- character(0)
  key <- vapply(split(sub$d, sub$m), function(x) paste(sort(x), collapse = ","), "")
  for (dsig in unique(key)) {
    ms <- as.integer(names(key)[key == dsig])
    ds <- as.integer(strsplit(dsig, ",")[[1]])
    kpart <- if (is.null(k_spec)) "" else sprintf(' "%s"', k_spec)
    out <- c(out, sprintf(
      'Rscript run-settings.R --data=simdata_nonsta.RData --setting=%d "%s" %d 2 "%s"%s',
      setting, spec(ds), workers, spec(ms), kpart))
  }
  out
}

# K = 0 baseline: its own call(s), no 5th arg
sub0 <- miss[miss$K == 0, ]
if (nrow(sub0) > 0) lines <- c(lines, emit_group(sub0, NULL))

# K > 0: merge K values whose missing (m,d) sets are identical, so a fresh
# chunk becomes a single call covering all K at once (one PSOCK cluster).
missK <- miss[miss$K > 0, ]
if (nrow(missK) > 0) {
  sig <- vapply(split(paste(missK$m, missK$d), missK$K),
                function(x) paste(sort(x), collapse = ";"), "")
  for (s in unique(sig)) {
    ks <- as.integer(names(sig)[sig == s])
    sub <- missK[missK$K == ks[1], ]
    lines <- c(lines, emit_group(sub, spec(ks)))
  }
}

writeLines(lines, calls_out)
cat(sprintf("MISSING %d\n", nrow(miss)))
cat(sprintf("CALLS %d\n", length(lines)))
