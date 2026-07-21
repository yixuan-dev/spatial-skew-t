# brier_bootstrap_400.R
# ---------------------------------------------------------------------------
# 400/400 (two-fold CV) companion to US-all-auto/brier_bootstrap.R, for the
# thesis's MAIN ozone comparison. Each results/us-all-<setting>.RData holds
# fit[[1]] and fit[[2]], the fits trained on one fold and predicting the other
# (cv.lst[[d]] = validation sites of fold d; the folds partition all 800
# sites), so pooling the two folds scores every site exactly once as a
# validation site.
#
# Same estimand and machinery as the auto-select splits: per-site Brier SSE
# via BrierScoreSite, paired difference per site, pooled ratio sum(SSE)/
# sum(days), site-clustered bootstrap. Resampling is STRATIFIED BY FOLD
# (sites resampled within their fold), since all sites in a fold share one
# training set and posterior.
#
# Output (kept beside the other splits' CSVs so the table generator sees all
# three): US-all-auto/output/us-all-auto/tables/brier_bootstrap_pairs_400-400.csv
# ---------------------------------------------------------------------------

rm(list = ls())
.this <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), winslash = "/", mustWork = FALSE)) else "."
})
setwd(.this)
.r_root <- normalizePath(file.path(.this, "../../../R"), winslash = "/", mustWork = FALSE)
source(local({ ar2 <- file.path(.r_root, "ar2", "auxfunctions.R"); if (file.exists(ar2)) ar2 else file.path(.r_root, "auxfunctions.R") }))

target_probs <- c(0.90, 0.95, 0.98, 0.99, 0.995)
B <- 2000L; seed <- 2024L
pairs <- list(
  c(111, 61,  "AR(2) vs AR(1)  | skew-t, K=7, T=50"),
  c(204, 51,  "MRTS(10) vs none | skew-t, K=1, T=0, TS")
)

se <- new.env(); load("us-all-setup.RData", envir = se)
Y <- get("Y", se); cv.lst <- get("cv.lst", se)
stopifnot(length(cv.lst) == 2L, length(intersect(cv.lst[[1]], cv.lst[[2]])) == 0L)
thresholds <- quantile(Y, probs = target_probs, na.rm = TRUE)

cat(sprintf("Two-fold CV: %d + %d validation sites (all %d)\n",
            length(cv.lst[[1]]), length(cv.lst[[2]]), nrow(Y)))
for (k in seq_along(target_probs))
  cat(sprintf("  p=%.3f -> %.2f ppb  (n_exceed=%d)\n", target_probs[k], thresholds[k],
              sum(Y > thresholds[k], na.rm = TRUE)))

# per-site SSE across both folds: site order = c(cv.lst[[1]], cv.lst[[2]])
site_sse <- function(sid) {
  fe <- new.env(); load(sprintf("results/us-all-%d.RData", sid), envir = fe)
  fit <- get("fit", fe)
  stopifnot(length(fit) >= 2L)
  out <- vector("list", 2L)
  for (d in 1:2) {
    val <- cv.lst[[d]]
    Yv <- Y[val, , drop = FALSE]
    nday_d <- rowSums(!is.na(Yv))
    per_site_mean <- BrierScoreSite(preds = fit[[d]]$yp, thresholds = thresholds, validate = Yv)
    out[[d]] <- per_site_mean * nday_d              # SSE = mean * ndays
  }
  cat(sprintf("  scored setting %-4d (folds 1+2, draws=%d)\n", sid, dim(fit[[1]]$yp)[1]))
  rbind(out[[1]], out[[2]])                          # 800 x nthresh, fold-1 sites first
}

need <- unique(unlist(lapply(pairs, function(p) as.integer(p[1:2]))))
have <- vapply(need, function(s) file.exists(sprintf("results/us-all-%d.RData", s)), logical(1))
if (any(!have)) cat("  [skip missing]", paste(need[!have], collapse = ", "), "\n")
SSE <- setNames(lapply(need[have], site_sse), as.character(need[have]))
pairs <- Filter(function(p) all(as.character(as.integer(p[1:2])) %in% names(SSE)), pairs)

sites_all <- c(cv.lst[[1]], cv.lst[[2]])
nday <- rowSums(!is.na(Y[sites_all, , drop = FALSE]))
fold_id <- rep(1:2, times = c(length(cv.lst[[1]]), length(cv.lst[[2]])))
idx_f1 <- which(fold_id == 1); idx_f2 <- which(fold_id == 2)

boot_pair <- function(sidA, sidB, k) {
  d <- SSE[[as.character(sidA)]][, k] - SSE[[as.character(sidB)]][, k]
  w <- nday
  ratio <- function(idx) sum(d[idx]) / sum(w[idx])
  point <- ratio(seq_along(d))
  set.seed(seed)
  bootv <- replicate(B, {
    idx <- c(sample(idx_f1, replace = TRUE), sample(idx_f2, replace = TRUE))  # stratified
    ratio(idx)
  })
  ci <- quantile(bootv, c(.025, .975), names = FALSE)
  p <- min(1, 2 * min(mean(bootv <= 0), mean(bootv >= 0)))
  data.frame(level = target_probs[k],
             BS_A = sum(SSE[[as.character(sidA)]][, k]) / sum(w),
             BS_B = sum(SSE[[as.character(sidB)]][, k]) / sum(w),
             Delta = point, CI_lo = ci[1], CI_hi = ci[2], boot_p = p,
             n_exceed = sum(Y > thresholds[k], na.rm = TRUE))
}

out <- do.call(rbind, lapply(pairs, function(p) {
  A <- as.integer(p[1]); Bx <- as.integer(p[2])
  rows <- do.call(rbind, lapply(seq_along(target_probs), function(k) boot_pair(A, Bx, k)))
  cbind(pair = p[3], A = A, B = Bx, rows)
}))

fmt <- out
for (c0 in c("BS_A", "BS_B", "Delta", "CI_lo", "CI_hi")) fmt[[c0]] <- sprintf("%+.5f", out[[c0]])
fmt$boot_p <- sprintf("%.3f", out$boot_p)
fmt$sig <- ifelse(out$CI_lo > 0 | out$CI_hi < 0, "*", "")
cat("\n=== 400/400 paired Brier differences (Delta = A - B; negative => A better) ===\n")
cat("=== 95% fold-stratified site bootstrap CI; '*' = CI excludes 0 ===\n\n")
print(fmt[, c("pair", "level", "BS_A", "BS_B", "Delta", "CI_lo", "CI_hi", "boot_p", "n_exceed", "sig")],
      row.names = FALSE)

out_dir <- normalizePath(file.path(.this, "../US-all-auto/output/us-all-auto/tables"),
                         winslash = "/", mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
outfile <- file.path(out_dir, "brier_bootstrap_pairs_400-400.csv")
write.csv(out, outfile, row.names = FALSE)
cat(sprintf("\nWritten: %s\n", outfile))
