# Scan scores5_0mrts.RData (baseline, no MRTS) and scores5_15mrts.RData (MRTS K=15)
# DGP for setting 5: Skew-t, K=5, lambda=3 -- expected matching method: 4

inspect_env <- function(path) {
  cat("\n========================================\n")
  cat(sprintf(" %s\n", path))
  cat("========================================\n")
  e <- new.env()
  load(path, envir = e)
  cat("objects:", paste(sort(ls(e)), collapse = ", "), "\n")
  for (nm in c("probs", "mrts_k", "setting", "methods", "datasets")) {
    if (exists(nm, envir = e, inherits = FALSE)) {
      v <- get(nm, envir = e, inherits = FALSE)
      cat(sprintf("  %-10s = %s\n", nm, paste(v, collapse = ",")))
    }
  }
  cat(sprintf("  dim(brier.score) = %s\n",
              paste(dim(e$brier.score), collapse = " x ")))
  e
}

report <- function(e, tag) {
  cat("\n========================================\n")
  cat(sprintf(" %s\n", tag))
  cat("========================================\n")
  bs <- e$brier.score
  probs <- e$probs
  methods <- e$methods
  # [n_probs, n_datasets, n_methods]
  mean_bs <- apply(bs, c(1, 3), mean, na.rm = TRUE)
  rownames(mean_bs) <- sprintf("q=%.3f", probs)
  colnames(mean_bs) <- sprintf("m%d", methods)
  cat("Mean absolute Brier:\n")
  print(round(mean_bs, 5))

  cat("\nRelative to Gaussian:\n")
  ref_col <- which(methods == 1L)
  if (length(ref_col) == 1L) {
    rel <- sweep(mean_bs, 1, mean_bs[, ref_col], "/")
    print(round(rel, 4))
  } else {
    cat("(no method 1 in this cache)\n")
  }

  cat("\nBest method per quantile:",
      methods[apply(mean_bs, 1, which.min)], "\n")

  invisible(mean_bs)
}

e0  <- inspect_env("scores5_0mrts.RData")
e15 <- inspect_env("scores5_15mrts.RData")

mb0  <- report(e0,  "scores5_0mrts  (setting 5, mrts_k = 0)")
mb15 <- report(e15, "scores5_15mrts (setting 5, mrts_k = 15)")

# ---- Side-by-side comparison for the matching method (m4) ----------
cat("\n========================================\n")
cat(" Side-by-side: relative Brier of m2 vs m4\n")
cat(" (DGP = Skew-t K=5; matching method should be m4)\n")
cat("========================================\n")
mk_rel <- function(mb) {
  sweep(mb, 1, mb[, "m1"], "/")
}
rel0  <- mk_rel(mb0)
rel15 <- mk_rel(mb15)
side <- cbind(
  m2_k0  = rel0[,  "m2"], m4_k0  = rel0[,  "m4"],
  m2_k15 = rel15[, "m2"], m4_k15 = rel15[, "m4"]
)
rownames(side) <- rownames(rel0)
print(round(side, 4))

# ---- Bonus: per-dataset count of "m4 best" vs "m2 best" -----------
cat("\n========================================\n")
cat(" Per-dataset best-method count (across 50 datasets, q=0.99)\n")
cat("========================================\n")
count_best <- function(e, tag) {
  bs <- e$brier.score  # [probs, datasets, methods]
  probs <- e$probs
  methods <- e$methods
  qi <- which(abs(probs - 0.99) < 1e-9)
  mat <- bs[qi, , ]   # [datasets, methods]
  best <- methods[apply(mat, 1, function(r) {
    if (all(is.na(r))) NA_integer_ else which.min(r)
  })]
  tbl <- table(factor(best, levels = methods))
  cat(sprintf("  %s : ", tag))
  print(tbl)
  invisible(tbl)
}
count_best(e0,  "k=0 ")
count_best(e15, "k=15")
