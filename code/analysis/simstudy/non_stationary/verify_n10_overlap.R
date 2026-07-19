#########################################################################
# verify_n10_overlap.R -- regression guard before deleting fits.
#
# After re-scoring a setting to n=10 / full-K, confirm that the cells that
# ALSO existed in the old n=3 cache (datasets 1:3, K in {0,30,40,50}) are
# bit-identical.  Same seed -> same fits -> same scores, so any mismatch
# means the scoring code drifted or a fit file was corrupted -- in which
# case we must NOT delete the fits.
#
# Usage:
#   Rscript verify_n10_overlap.R <setting> <old_backup.RData> <new_cache.RData>
# Exit code 0 = match (safe to delete fits); 1 = mismatch (keep fits).
#########################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("need: <setting> <old_backup> <new_cache>")
setting <- as.integer(args[1])
old_f   <- args[2]
new_f   <- args[3]
tol     <- 1e-6

to3 <- function(a) if (length(dim(a)) == 4) apply(a, c(2, 3, 4), mean, na.rm = TRUE) else a

eo <- new.env(); load(old_f, envir = eo)
en <- new.env(); load(new_f, envir = en)

ds_o <- eo$datasets; ds_n <- en$datasets
if (!all(ds_o %in% ds_n)) stop("old datasets not a subset of new datasets")

fields <- c("brier.score", "energy.score", "vario.score", "pred.rmse", "recovery.rmse")
worst <- 0
ok <- TRUE
for (nm in fields) {
  ao <- to3(get(nm, envir = eo))   # [ds, methods, ks]
  an <- to3(get(nm, envir = en))
  ks_o <- as.numeric(dimnames(ao)[[3]]); ks_n <- as.numeric(dimnames(an)[[3]])
  ks <- intersect(ks_o, ks_n)
  di <- match(ds_o, ds_n)
  ki_o <- match(ks, ks_o); ki_n <- match(ks, ks_n)
  sub_o <- ao[, , ki_o, drop = FALSE]
  sub_n <- an[di, , ki_n, drop = FALSE]
  d <- abs(sub_o - sub_n)
  mx <- max(d[is.finite(d)], 0)
  na_mismatch <- sum(is.na(sub_o) != is.na(sub_n))
  worst <- max(worst, mx)
  flag <- if (mx > tol || na_mismatch > 0) { ok <- FALSE; "FAIL" } else "ok"
  cat(sprintf("  %-14s overlap ds=%s ks={%s}  max|diff|=%.3g  NA-mismatch=%d  [%s]\n",
              nm, paste(ds_o, collapse = ","), paste(ks, collapse = ","),
              mx, na_mismatch, flag))
}
cat(sprintf("setting %d: worst |diff| = %.3g over overlap cells\n", setting, worst))
if (!ok) { cat("MISMATCH -- do NOT delete fits\n"); quit(status = 1) }
cat("OVERLAP VERIFIED -- safe to delete fits\n"); quit(status = 0)
