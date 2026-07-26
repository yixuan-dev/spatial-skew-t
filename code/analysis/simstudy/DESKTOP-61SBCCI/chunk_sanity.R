#########################################################################
# chunk_sanity.R -- fast validity check of one chunk score cache; the
# n30 driver only deletes a chunk's fit files after this passes.
#
# Usage (from code/analysis/simstudy):
#   Rscript DESKTOP-61SBCCI/chunk_sanity.R <chunk.RData> <setting> <datasets>
#   e.g. Rscript DESKTOP-61SBCCI/chunk_sanity.R \
#          output/results/scores20_nonsta_d11-15.RData 20 11:15
#
# Checks (rules derived from the verified n=10 production caches):
#   - setting / datasets / methods(1:5) / mrts_ks (full 12-K grid) match
#   - score arrays, elapsed_sec and parameter-interval arrays contain no
#     NA -- except `lambda`, which is legitimately NA for the non-skew
#     methods 1/3/5 and must be non-NA only for the skew-t methods 2/4
#   - elapsed_sec strictly positive
#########################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("usage: chunk_sanity.R <chunk.RData> <setting> <datasets e.g. 11:15>",
       call. = FALSE)
}
cache_file <- args[1]
setting_expect <- as.integer(args[2])
datasets_expect <- sort(as.integer(eval(parse(text = args[3]))))

fail <- function(...) { cat("[SANITY FAILED] ", sprintf(...), "\n", sep = ""); quit(status = 1) }

if (!file.exists(cache_file)) fail("cache not found: %s", cache_file)
e <- new.env()
load(cache_file, envir = e)

if (!identical(as.integer(e$setting), setting_expect)) {
  fail("setting is %s, expected %d", e$setting, setting_expect)
}
if (!identical(sort(as.integer(e$datasets)), datasets_expect)) {
  fail("datasets are [%s], expected [%s]",
       paste(e$datasets, collapse = ","), paste(datasets_expect, collapse = ","))
}
if (!identical(as.integer(e$methods), 1:5)) {
  fail("methods are [%s], expected 1:5", paste(e$methods, collapse = ","))
}
K_full <- c(0L, 3L, 5L, 8L, 10L, 12L, 15L, 20L, 25L, 30L, 40L, 50L)
if (!identical(as.integer(e$mrts_ks), K_full)) {
  fail("mrts_ks are [%s], expected full 12-K grid", paste(e$mrts_ks, collapse = ","))
}

no_na <- c("brier.score", "quant.score", "energy.score", "vario.score",
           "crps.score", "pred.rmse", "recovery.rmse", "elapsed_sec",
           "beta.0", "gamma", "nu", "rho",
           "tau.alpha", "tau.beta")
for (nm in no_na) {
  if (!exists(nm, envir = e)) fail("object `%s` missing from cache", nm)
  o <- get(nm, envir = e)
  if (anyNA(o)) fail("`%s` has %d NA cells", nm, sum(is.na(o)))
}

# beta.1 / beta.2 exist only in the K=0 baseline design [1, s1, s2]; the
# MRTS fits use mrts(S, k) as the full design (p_base = 0), so K>0 cells
# are legitimately NA (same pattern as the verified n=10 caches).
for (nm in c("beta.1", "beta.2")) {
  if (!exists(nm, envir = e)) fail("object `%s` missing from cache", nm)
  o <- get(nm, envir = e)
  k0 <- o[, , , "0", drop = FALSE]
  if (anyNA(k0)) fail("`%s` has %d NA cells at K=0", nm, sum(is.na(k0)))
}

# lambda: NA allowed only for non-skew methods (1/3/5)
lam <- e$lambda
skew_methods <- intersect(c("2", "4"), dimnames(lam)[["method"]])
lam_skew <- lam[, , skew_methods, , drop = FALSE]
if (anyNA(lam_skew)) fail("`lambda` has NA for skew-t methods 2/4")

if (any(e$elapsed_sec <= 0)) fail("elapsed_sec has non-positive entries")

cat(sprintf("[SANITY OK] %s: setting %d, datasets %s, %d methods x %d K\n",
            cache_file, setting_expect,
            paste(range(datasets_expect), collapse = ".."),
            length(e$methods), length(e$mrts_ks)))
