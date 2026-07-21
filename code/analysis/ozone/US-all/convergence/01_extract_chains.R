# convergence/01_extract_chains.R
# ---------------------------------------------------------------------------
# Extract the scalar posterior chains from the big results/us-all-<N>.RData
# fits (0.9-1.8 GB each) into small per-setting caches so the diagnostics
# and plotting scripts never pay the GB-scale loads again.
#
# Scope: the nine thesis headline settings (see 00_conv_lib.R), overridable
# via US_ALL_CONV_SETTINGS as a comma list, e.g. "204,111".
# Re-extraction: existing caches are skipped unless US_ALL_CONV_FORCE=1.
#
# Output: output/us-all/results/convergence/chains-us-all-<N>.rds
# ---------------------------------------------------------------------------

rm(list = ls())
.this <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), winslash = "/", mustWork = FALSE)) else "."
})
setwd(file.path(.this, ".."))  # run from the US-all directory
source("convergence/00_conv_lib.R")

paths <- conv_paths()

env_settings <- trimws(Sys.getenv("US_ALL_CONV_SETTINGS", unset = ""))
targets <- if (nzchar(env_settings)) {
  as.integer(trimws(unlist(strsplit(env_settings, ",", fixed = TRUE))))
} else {
  conv_settings$setting
}
stopifnot(all(targets %in% conv_settings$setting))
force <- Sys.getenv("US_ALL_CONV_FORCE", unset = "0") %in% c("1", "true", "yes")

t0 <- proc.time()[3]
for (s in targets) {
  cf <- cache_file(s, paths)
  if (file.exists(cf) && !force) {
    cat(sprintf("setting %d: cache exists, skipping (%s)\n", s, cf))
    next
  }
  slim <- extract_setting_from_rdata(s, paths)
  saveRDS(slim, cf)
  cat(sprintf("  wrote %s (%.1f MB)\n", cf, file.size(cf) / 2^20))
}
cat(sprintf("done: %d setting(s) in %.1f min\n", length(targets), (proc.time()[3] - t0) / 60))
