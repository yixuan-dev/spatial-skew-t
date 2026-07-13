#########################################################################
# THE HEADLINE FIGURE: nonsta setting 1 vs setting 6.
#
# Both settings are driven by the SAME cosine field f1(s).  Setting 1 puts it in
# the MEAN (g(s) = 5 f1 + 3 f2); setting 6 puts it in the correlation RANGE
# (rho(s) = rho0 * exp(0.7 f1(s))).  f1 is a smooth function of location, so the
# MRTS basis spans it -- setting 1 proves exactly that -- and the model has the
# coordinates either way.  So if MRTS recovers setting 1 and fails on setting 6,
# the failure cannot be "the model lacked the information".  It is the MOMENT
# that is wrong: mean-only augmentation cannot reach the second moment.
#
# Everything is plotted RELATIVE TO THE SAME METHOD'S OWN K=0 FIT, per dataset.
# Two reasons:
#
#   1. It is the question actually being asked (see non_stationary/FUCK.md):
#      "does adding MRTS covariates help THIS method?", not "does this method
#      beat Gaussian?".  A ratio against the same method's no-MRTS baseline
#      answers the first; a ratio against Gaussian answers neither.
#
#   2. recovery.rmse has a per-method CONSTANT OFFSET that a ratio cancels.  The
#      true mean is 10 + structure, but the fitted intercept also absorbs
#      lambda*E[z_t] ~ 5.4 (the skew term's mean).  So even where the truth is
#      dead flat -- setting 6, where the true mean IS exactly 10 everywhere --
#      the absolute recovery.rmse sits around 4, not 0.  That offset is constant
#      in K, so score(K)/score(0) removes it.
#
# The K comparison is PAIRED: get_simstudy_seed() deliberately excludes mrts_k,
# so every K for a given (method, dataset) shares an RNG stream.  We therefore
# take the ratio WITHIN a dataset and average the ratios, not the other way
# round -- that keeps the pairing and cuts the variance.
#
# Output: output/plots/headline_nonsta_1v6.pdf  +  a printed summary table.
#########################################################################

if (basename(getwd()) != "simstudy") setwd("d:/Github/spatial-skew-t/code/analysis/simstudy")

dir.create("output/plots", recursive = TRUE, showWarnings = FALSE)
out_pdf <- "output/plots/headline_nonsta_1v6.pdf"

load_scores <- function(setting) {
  fp <- sprintf("output/results/scores%d_nonsta.RData", setting)
  if (!file.exists(fp)) stop(sprintf("missing %s -- run scores.R --setting=%d first", fp, setting))
  e <- new.env(parent = emptyenv())
  load(fp, envir = e)
  e
}

s1 <- load_scores(1)
s6 <- load_scores(6)

method_lab <- c("1 Gaussian", "2 Skew-t K=1", "3 Sym-t K=1 T80",
                "4 Skew-t K=5", "5 Sym-t K=5 T80")
mcol <- c("#1b6ca8", "#d1495b", "#66a182", "#e8a33d", "#8367c7")

# ---- which datasets are actually fitted in BOTH settings ---------------
# setting 1 was run on datasets 1-10; the setting-6 pilot on 1-3. Only the
# intersection is a fair paired comparison.
fitted_ds <- function(e) {
  a <- e$recovery.rmse                       # [dataset, method, K]
  which(apply(a, 1, function(v) any(!is.na(v))))
}
ds <- intersect(fitted_ds(s1), fitted_ds(s6))
cat(sprintf("datasets fitted in setting 1 : %s\n", paste(fitted_ds(s1), collapse = ",")))
cat(sprintf("datasets fitted in setting 6 : %s\n", paste(fitted_ds(s6), collapse = ",")))
cat(sprintf("=> paired comparison on datasets: %s  (n = %d)\n\n",
            paste(ds, collapse = ","), length(ds)))
if (length(ds) == 0) stop("no datasets fitted in both settings")

Ks <- as.numeric(dimnames(s1$recovery.rmse)[[3]])
stopifnot("K grids differ between settings 1 and 6" =
            identical(dimnames(s1$recovery.rmse)[[3]], dimnames(s6$recovery.rmse)[[3]]))
k0 <- which(Ks == 0)

# ---- relative curves ---------------------------------------------------
# rel[m, k] = mean over datasets of  score[d, m, k] / score[d, m, K=0]
rel_recovery <- function(e) {
  a <- e$recovery.rmse[ds, , , drop = FALSE]          # [d, m, K]
  r <- sweep(a, c(1, 2), a[, , k0], "/")
  apply(r, c(2, 3), mean, na.rm = TRUE)               # [m, K]
}

# brier.score is [prob, dataset, method, K]; average over the 11 exceedance
# probabilities first, then take the paired ratio.
rel_brier <- function(e, probs_use = NULL) {
  b <- e$brier.score
  pi_use <- if (is.null(probs_use)) seq_len(dim(b)[1]) else which(dimnames(b)[[1]] %in% probs_use)
  a <- apply(b[pi_use, ds, , , drop = FALSE], c(2, 3, 4), mean, na.rm = TRUE)  # [d, m, K]
  r <- sweep(a, c(1, 2), a[, , k0], "/")
  apply(r, c(2, 3), mean, na.rm = TRUE)
}

R1 <- rel_recovery(s1); R6 <- rel_recovery(s6)
B1 <- rel_brier(s1);    B6 <- rel_brier(s6)
# the far tail, where MRTS is most likely to be oversold
B1t <- rel_brier(s1, c("0.99", "0.995"))
B6t <- rel_brier(s6, c("0.99", "0.995"))

# ---- plot --------------------------------------------------------------
panel <- function(M, ttl, ylab, ylim = NULL) {
  if (is.null(ylim)) ylim <- range(c(1, M), na.rm = TRUE) * c(0.98, 1.02)
  matplot(Ks, t(M), type = "b", pch = 16, lty = 1, lwd = 2, col = mcol,
          xlab = expression(K[MRTS]), ylab = ylab, ylim = ylim,
          main = ttl, cex.main = 1.0)
  abline(h = 1, col = "grey40", lty = 2)   # = no better than the K=0 baseline
}

pdf(out_pdf, width = 12, height = 8.5)
par(mfrow = c(2, 2), mar = c(4.5, 4.8, 3.5, 1), oma = c(3.2, 0, 3.2, 0))

ylim_rec <- range(c(1, R1, R6), na.rm = TRUE) * c(0.95, 1.05)
panel(R1, "Setting 1: f1 in the MEAN\n(MRTS should recover it)",
      "relative recovery RMSE", ylim_rec)
panel(R6, "Setting 6: the SAME f1 in the RANGE\n(MRTS structurally cannot)",
      "relative recovery RMSE", ylim_rec)

ylim_bs <- range(c(1, B1, B6), na.rm = TRUE) * c(0.98, 1.02)
panel(B1, "Setting 1: Brier (mean over exceedance probs)",
      "relative Brier score", ylim_bs)
panel(B6, "Setting 6: Brier (mean over exceedance probs)",
      "relative Brier score", ylim_bs)

legend("bottomright", method_lab, col = mcol, lty = 1, lwd = 2, pch = 16,
       bty = "n", cex = 0.85)
mtext(sprintf("Non-stationary sim: the SAME cosine field f1(s), placed in a different moment  (paired, datasets %s)",
              paste(range(ds), collapse = "-")),
      outer = TRUE, line = 0.8, cex = 1.05, font = 2)
mtext("Everything relative to the same method's own K=0 fit. Below the dashed line = MRTS helps.",
      outer = TRUE, side = 1, line = 1.0, cex = 0.9)

# tail-only Brier, on its own page
par(mfrow = c(1, 2), mar = c(4.5, 4.8, 3.5, 1), oma = c(3.2, 0, 3.2, 0))
ylim_t <- range(c(1, B1t, B6t), na.rm = TRUE) * c(0.98, 1.02)
panel(B1t, "Setting 1: Brier at q = 0.99 / 0.995", "relative Brier score", ylim_t)
panel(B6t, "Setting 6: Brier at q = 0.99 / 0.995", "relative Brier score", ylim_t)
legend("bottomright", method_lab, col = mcol, lty = 1, lwd = 2, pch = 16, bty = "n", cex = 0.85)
mtext("Far-tail Brier -- where a mean-only basis is most likely to be oversold",
      outer = TRUE, line = 0.8, cex = 1.05, font = 2)

dev.off()

# ---- summary -----------------------------------------------------------
fmt <- function(M, ttl) {
  cat("\n===", ttl, "===\n")
  d <- as.data.frame(round(M, 3))
  colnames(d) <- paste0("K", Ks)
  rownames(d) <- method_lab
  print(d)
}
fmt(R1, "setting 1 -- relative recovery RMSE (vs own K=0)")
fmt(R6, "setting 6 -- relative recovery RMSE (vs own K=0)")
fmt(B1, "setting 1 -- relative Brier (vs own K=0)")
fmt(B6, "setting 6 -- relative Brier (vs own K=0)")

cat("\n=== THE CLAIM, in two numbers ===\n")
best <- function(M) apply(M[, -k0, drop = FALSE], 1, min)   # best relative score over K > 0
cat("best relative recovery RMSE achieved by ANY K > 0:\n")
cmp <- data.frame(setting1 = round(best(R1), 3), setting6 = round(best(R6), 3))
rownames(cmp) <- method_lab
print(cmp)
cat("\n  setting 1 << 1  => MRTS recovers the mean structure\n")
cat("  setting 6 ~= 1  => nothing to recover: the same f1 lives in the SECOND moment\n")
cat(sprintf("\nwrote %s\n", out_pdf))
