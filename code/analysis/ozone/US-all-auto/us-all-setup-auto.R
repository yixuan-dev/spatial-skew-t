# us-all-setup-auto.R
#
# Loads the shared ozone dataset and produces us-all-setup-auto.RData with a
# reproducible train / validation / test site partition.
#
# Saved objects:
#   Y, X, S      - full dataset (row indices = site IDs)
#   beta.init, tau.init - MCMC initialisation values
#   split.lst    - list(train = int[300], val = int[100], test = int[400])
#                  row indices into Y / X / S
#
# How to run (PowerShell):
#   cd D:\Github\spatial-skew-t\code\analysis\ozone\US-all-auto
#   Rscript us-all-setup-auto.R
#
# Optional environment variables:
#   US_ALL_AUTOSELECT_N_TRAIN   default 300
#   US_ALL_AUTOSELECT_N_VAL     default 100
#   US_ALL_AUTOSELECT_N_TEST    default 400
#   US_ALL_AUTOSELECT_SEED      default 2024
#   US_ALL_AUTOSELECT_OUT       output filename (default: us-all-setup-auto.RData)

rm(list = ls())

# ---- Locate this script ----

.this_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg) > 0)
    dirname(normalizePath(sub("^--file=", "", arg[1]), winslash = "/", mustWork = FALSE))
  else
    normalizePath(".", winslash = "/", mustWork = FALSE)
})
setwd(.this_dir)

.us_all_root <- normalizePath(file.path(.this_dir, "../US-all"),
                               winslash = "/", mustWork = FALSE)

source(file.path(.this_dir, "utils.R"))

# ---- Helper ----

parse_env_int <- function(name, default) {
  val <- suppressWarnings(as.integer(trimws(Sys.getenv(name, unset = ""))))
  if (length(val) == 1L && is.finite(val) && val > 0L) val else as.integer(default)
}

# ---- Split configuration ----

n_train    <- parse_env_int("US_ALL_AUTOSELECT_N_TRAIN", 300L)
n_val      <- parse_env_int("US_ALL_AUTOSELECT_N_VAL",   100L)
n_test     <- parse_env_int("US_ALL_AUTOSELECT_N_TEST",  400L)
split_seed <- parse_env_int("US_ALL_AUTOSELECT_SEED",   2024L)

cat("=== us-all-setup-auto.R ===\n")
cat(sprintf("Split config : %d train | %d val | %d test  (seed %d)\n\n",
    n_train, n_val, n_test, split_seed))

# ---- Locate and load shared dataset ----
# Prefer a local copy in US-all-auto/; fall back to the sibling US-all/.

setup_path <- local({
  local_copy <- file.path(.this_dir,   "us-all-setup.RData")
  sibling    <- file.path(.us_all_root, "us-all-setup.RData")
  if (file.exists(local_copy)) local_copy else sibling
})

if (!file.exists(setup_path))
  stop("us-all-setup.RData not found in US-all-auto/ or US-all/.\n",
       "Copy the file here or run US-all/us-all-setup.R first.")

cat("Loading:", setup_path, "\n")
src <- new.env(parent = emptyenv())
load(setup_path, envir = src)

for (obj in c("Y", "X", "S", "cv.lst", "beta.init", "tau.init")) {
  if (!exists(obj, envir = src, inherits = FALSE))
    stop("Required object missing from us-all-setup.RData: ", obj)
}

Y        <- get("Y",        envir = src)
X        <- get("X",        envir = src)
S        <- get("S",        envir = src)
cv.lst   <- get("cv.lst",   envir = src)
beta.init <- get("beta.init", envir = src)
tau.init  <- get("tau.init",  envir = src)

cat(sprintf("Dataset loaded : %d sites x %d time points\n", nrow(Y), ncol(Y)))

# ---- Build split.lst ----
# Pool all site indices covered by cv.lst, then partition randomly.

all_sites <- sort(unique(unlist(cv.lst)))
n_sites   <- length(all_sites)
cat(sprintf("CV-covered sites: %d\n\n", n_sites))

if (n_train + n_val + n_test != n_sites) {
  n_test <- n_sites - n_train - n_val
  message(sprintf("Adjusting n_test to %d so that train+val+test = %d.", n_test, n_sites))
}
if (n_val <= 0L || n_test <= 0L)
  stop("n_val and n_test must both be positive after size adjustment.")

set.seed(split_seed)
perm <- sample(all_sites)

split.lst <- list(
  train = sort(perm[seq_len(n_train)]),
  val   = sort(perm[seq.int(n_train + 1L, n_train + n_val)]),
  test  = sort(perm[seq.int(n_train + n_val + 1L, n_sites)])
)

stopifnot(length(intersect(split.lst$train, split.lst$val))  == 0L)
stopifnot(length(intersect(split.lst$train, split.lst$test)) == 0L)
stopifnot(length(intersect(split.lst$val,   split.lst$test)) == 0L)
stopifnot(length(split.lst$train) + length(split.lst$val) + length(split.lst$test) == n_sites)

cat(sprintf("split.lst : %d train | %d val | %d test\n",
    length(split.lst$train), length(split.lst$val), length(split.lst$test)))
cat("Sanity checks passed (no overlap, correct total).\n\n")

# ---- Exceedance counts by split ----

print_exceedance_summary(Y, split.lst)

# ---- Save ----

.out_file <- trimws(Sys.getenv("US_ALL_AUTOSELECT_OUT", unset = ""))
if (!nzchar(.out_file))
  .out_file <- sprintf("us-all-setup-auto-%d-%d-%d.RData",
                       length(split.lst$train),
                       length(split.lst$val),
                       length(split.lst$test))
out_path <- file.path(.this_dir, .out_file)
if (file.exists(out_path))
  stop("Output file already exists: ", out_path,
       "\nSet US_ALL_AUTOSELECT_OUT to a different filename.", call. = FALSE)
save(Y, X, S, beta.init, tau.init, split.lst, file = out_path)

cat("Saved:", out_path, "\n")
cat("Objects: Y, X, S, beta.init, tau.init, split.lst\n")
cat("\nRun next:  Rscript run_settings_val.R\n")
