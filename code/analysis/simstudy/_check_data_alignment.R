# Diagnose: does fit$y (frozen when fit was created) match current simdata.RData?
load("simdata.RData")
y_new <- y; ntest_new <- ntest

load("simdata_old.RData")
y_old <- y; ntest_old <- ntest

cat("simdata.RData      dim(y) =", paste(dim(y_new), collapse = "x"),
    " ntest =", ntest_new, "\n")
cat("simdata_old.RData  dim(y) =", paste(dim(y_old), collapse = "x"),
    " ntest =", ntest_old, "\n")

# Compare y for setting 5, dataset 1
diff_new_old <- max(abs(y_new[, , 1, 5] - y_old[, , 1, 5]), na.rm = TRUE)
cat("\nmax |y_new - y_old| at (setting=5, dataset=1) =", diff_new_old, "\n")

source("../../R/auxfunctions.R")
probs <- c(0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995)

# Compute Brier for setting=5 dataset=1 across all 5 methods, twice:
# once against simdata_OLD (what fits were trained on) and once against NEW.
# Compare with scores5_0mrts.RData (computed May 7 against OLD data).
run_one <- function(y_arr, ntest_v, label) {
  obs <- c(rep(TRUE, nrow(y_arr) - ntest_v), rep(FALSE, ntest_v))
  thr <- quantile(y_arr[, , 1, 5], probs = probs, na.rm = TRUE, names = FALSE)
  val <- y_arr[!obs, , 1, 5]
  res <- matrix(NA_real_, nrow = 5, ncol = length(probs),
                dimnames = list(sprintf("m%d", 1:5), sprintf("q%.3f", probs)))
  for (m in 1:5) {
    f <- sprintf("results/5-%d-1.RData", m)
    if (!file.exists(f)) next
    e <- new.env(); load(f, envir = e)
    fit <- e$fit.1
    if (is.null(fit$yp)) next
    res[m, ] <- as.numeric(BrierScore(fit$yp, thr, val))
  }
  cat(sprintf("\n[%s] Brier @ q=0.90/0.95/0.99/0.995 (dataset 1, all 5 methods):\n",
              label))
  print(round(res[, c(1, 6, 10, 11)], 5))
  invisible(res)
}
bs_old <- run_one(y_old, ntest_old, "simdata_OLD")
bs_new <- run_one(y_new, ntest_new, "simdata_NEW")

# Compare to scores5_0mrts (was computed against the data version that
# existed on May 7 -- should match OLD if simdata_old is that snapshot).
e_s <- new.env(); load("scores5_0mrts.RData", envir = e_s)
bs_cache <- t(e_s$brier.score[, 1, ])  # [methods, probs] for dataset 1
rownames(bs_cache) <- sprintf("m%d", e_s$methods)
colnames(bs_cache) <- sprintf("q%.3f", e_s$probs)
cat("\n[scores5_0mrts cached] Brier dataset 1, all 5 methods:\n")
print(round(bs_cache[, c(1, 6, 10, 11)], 5))

cat("\nmax |OLD - cache| =", max(abs(bs_old - bs_cache), na.rm = TRUE), "\n")
cat("max |NEW - cache| =", max(abs(bs_new - bs_cache), na.rm = TRUE), "\n")
