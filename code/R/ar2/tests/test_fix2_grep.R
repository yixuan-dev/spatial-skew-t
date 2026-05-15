# ============================================================================
# Structural smoke test for Fix 2.
#
# Cheap insurance against the "someone deleted one of the t+2 blocks" failure
# mode. The patched code adds four `(t + 2) <= nt` conditionals in
# update_params_ar2.R (one per place y_t appears as a lag-2 lag in another
# day's conditional density):
#   1. updateTauTS_AR2,  nknots == 1 branch
#   2. updateTauTS_AR2,  nknots >  1 branch
#   3. updateZTS_AR2
#   4. updateKnotsTS_AR2
#
# This test asserts exactly four matches. It does NOT verify the contents of
# each block -- typos inside a block are not caught here; the math test
# (test_fix2_lag2_mh.R) and the integration test (test_fix2_production.R)
# cover that for one call site. The grep is the broadest, cheapest layer of
# regression protection across all four call sites.
# ============================================================================

script_path <- (function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
        return(normalizePath(sub("^--file=", "", file_arg[1])))
    }
    return(normalizePath(sys.frames()[[1]]$ofile))
})()
ar2_dir <- normalizePath(file.path(dirname(script_path), ".."))
src_path <- file.path(ar2_dir, "update_params_ar2.R")

src <- readLines(src_path)

fix2.pattern  <- "\\(t \\+ 2\\) <= nt"
expected      <- 4

fix2.hits <- which(grepl(fix2.pattern, src))
n.hits    <- length(fix2.hits)

cat(sprintf("Source: %s\n", src_path))
cat(sprintf("Pattern: %s\n", fix2.pattern))
cat(sprintf("Hits:    %d  (expected %d)\n", n.hits, expected))
if (n.hits > 0) {
    cat("Lines:\n")
    for (ln in fix2.hits) {
        cat(sprintf("  L%-4d  %s\n", ln, trimws(src[ln])))
    }
}

if (n.hits == expected) {
    cat("\nPASS: all four t+2 correction blocks present in update_params_ar2.R\n")
    quit(status = 0)
} else if (n.hits < expected) {
    cat(sprintf("\nFAIL: %d t+2 correction blocks missing -- Fix 2 has been\n",
                expected - n.hits))
    cat("      partially regressed. Re-apply the lag-2 MH term in the\n")
    cat("      offending update_*TS_AR2 function(s).\n")
    quit(status = 1)
} else {
    cat(sprintf("\nFAIL: expected %d matches but found %d. Either a new call\n",
                expected, n.hits))
    cat("      site was added (update this test to bump 'expected') or a\n")
    cat("      pattern collision occurred (refine 'fix2.pattern').\n")
    quit(status = 1)
}
