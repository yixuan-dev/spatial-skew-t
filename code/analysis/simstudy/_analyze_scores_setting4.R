# Compare scores4_0mrts.RData (new schema, mrts_k=0 baseline only)
#   vs scores4.RData (legacy 4D: probs x datasets x methods(6) x settings(8))
# DGP for setting 4: Skew-t, K=1, lambda=3
# Expected matching method: 2 (Skew-t, K=1)

inspect_env <- function(path) {
  cat("\n========================================\n")
  cat(sprintf(" %s\n", path))
  cat("========================================\n")
  e <- new.env()
  load(path, envir = e)
  cat("objects:", paste(sort(ls(e)), collapse = ", "), "\n")
  e
}

e0 <- inspect_env("scores4_0mrts.RData")
e1 <- inspect_env("scores4.RData")

# === scores4_0mrts: [probs, datasets, methods] =================
cat("\n========================================\n")
cat(" scores4_0mrts (baseline, methods 1..5)\n")
cat("========================================\n")
bs0 <- e0$brier.score
probs0 <- e0$probs
methods0 <- e0$methods
mean_bs0 <- apply(bs0, c(1, 3), mean, na.rm = TRUE)
rownames(mean_bs0) <- sprintf("q=%.3f", probs0)
colnames(mean_bs0) <- sprintf("m%d", methods0)
cat("Mean absolute Brier:\n")
print(round(mean_bs0, 5))
cat("\nRelative to Gaussian:\n")
print(round(sweep(mean_bs0, 1, mean_bs0[, 1], "/"), 4))
cat("\nBest method per quantile:",
    methods0[apply(mean_bs0, 1, which.min)], "\n")

# === scores4: [probs(11), datasets(50), methods(6), settings(8)] ===
cat("\n========================================\n")
cat(" scores4 legacy — setting=4 slice (DGP Skew-t K=1)\n")
cat("========================================\n")
bs1 <- e1$brier.score
probs1 <- e1$probs
cat("dim(brier.score) =", paste(dim(bs1), collapse = " x "), "\n")
cat("axes: probs(11) x datasets(50) x methods(6) x settings(8)\n\n")

setting_target <- 4L
n_settings1 <- dim(bs1)[4]
n_methods1  <- dim(bs1)[3]

# slice the target setting -> [probs, datasets, methods]
bs1_set4 <- bs1[, , , setting_target, drop = FALSE]
dim(bs1_set4) <- dim(bs1_set4)[1:3]
mean_bs1 <- apply(bs1_set4, c(1, 3), mean, na.rm = TRUE)
rownames(mean_bs1) <- sprintf("q=%.3f", probs1)
colnames(mean_bs1) <- sprintf("m%d", seq_len(n_methods1))
cat("Mean absolute Brier (setting 4 slice):\n")
print(round(mean_bs1, 5))
cat("\nRelative to Gaussian:\n")
print(round(sweep(mean_bs1, 1, mean_bs1[, 1], "/"), 4))
cat("\nBest method per quantile:",
    apply(mean_bs1, 1, function(r) {
      if (all(is.na(r))) NA_integer_ else which.min(r)
    }), "\n")

# === Reconciliation ============================================
cat("\n========================================\n")
cat(" Reconciliation: scores4_0mrts vs scores4[,,, setting=4]\n")
cat("========================================\n")
common_m <- intersect(methods0, seq_len(n_methods1))
mean0_common <- mean_bs0[, sprintf("m%d", common_m)]
mean1_common <- mean_bs1[, sprintf("m%d", common_m)]
d <- mean0_common - mean1_common
cat(sprintf("methods compared: %s\n", paste(common_m, collapse = ", ")))
cat(sprintf("max |diff|  = %.6g\n", max(abs(d), na.rm = TRUE)))
cat(sprintf("mean |diff| = %.6g\n", mean(abs(d), na.rm = TRUE)))

# === Winner across all settings in legacy cache ================
cat("\n========================================\n")
cat(" Winning method per setting (legacy cache)\n")
cat("========================================\n")
dgp_label <- c("Gaussian","t K=1","t K=5","Skew-t K=1","Skew-t K=5",
               "Max-stable","Skew-t (sub-thresh exp)","Brown-Resnick")
for (s in seq_len(n_settings1)) {
  arr <- bs1[, , , s, drop = FALSE]
  dim(arr) <- dim(arr)[1:3]
  m_mat <- apply(arr, c(1, 3), mean, na.rm = TRUE)
  if (all(is.na(m_mat))) {
    cat(sprintf("  setting %d (%s): all NA\n", s, dgp_label[s]))
    next
  }
  best <- apply(m_mat, 1, function(r) {
    if (all(is.na(r))) NA_integer_ else which.min(r)
  })
  cat(sprintf("  setting %d (%-25s): best per q = %s\n",
              s, dgp_label[s], paste(best, collapse = " ")))
}
