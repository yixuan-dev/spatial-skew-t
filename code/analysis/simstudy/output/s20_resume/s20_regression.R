# Regression check for the setting-20 n=10 rescore: the new cache's
# datasets 1:3 x K{0,30,40,50} cells must reproduce the old n=3 cache
# exactly (same seed -> same fit -> same score).  Nonzero exit on failure.
setwd("d:/Github/spatial-skew-t/code/analysis/simstudy")
en <- new.env(); load("output/results/scores20_nonsta.RData",        envir = en)
eo <- new.env(); load("output/results/scores20_nonsta.n3.bak.RData", envir = eo)

ok <- TRUE
for (nm in c("brier.score", "quant.score", "energy.score", "vario.score",
             "pred.rmse", "recovery.rmse")) {
  a <- get(nm, envir = en); b <- get(nm, envir = eo)
  nda <- length(dim(a))
  # dataset & K positions in the NEW cache matching the OLD one
  di <- match(eo$datasets, en$datasets)
  ki <- match(dimnames(b)[[nda]], dimnames(a)[[nda]])
  sub <- if (nda == 4) a[, di, , ki, drop = FALSE] else a[di, , ki, drop = FALSE]
  same <- isTRUE(all.equal(as.vector(sub), as.vector(b), tolerance = 1e-10))
  cat(sprintf("%-14s : %s\n", nm, if (same) "MATCH" else "** MISMATCH **"))
  if (!same) ok <- FALSE
}
if (!ok) { cat("REGRESSION CHECK FAILED\n"); quit(status = 1) }
cat("REGRESSION CHECK PASSED\n")
