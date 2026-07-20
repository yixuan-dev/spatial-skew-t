# timeblock_twcrps_export.R
# ---------------------------------------------------------------------------
# Tail-band twCRPS for the time-block forecast arm, to test whether the plain-
# CRPS AR(2) gain is a TAIL effect (goal-relevant: exceedance forecasting) or a
# BULK effect. Mirrors timeblock_crps_export.R's loading + clean-n=10 protocol,
# but the per-cell score is the tail-band integral of Brier
#   twCRPS(F,y) = integral_[q.90,q.99] (F(z) - 1{y<=z})^2 dz
# on a fine z-grid (trapezoid), z-values from the validation block's quantiles.
#
# Writes: output/tables/blk1_paired_ci_twcrps.csv
# ---------------------------------------------------------------------------

rm(list = ls())
.this <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), winslash = "/", mustWork = FALSE)) else "."
})
setwd(.this)

e <- new.env(); load("output/tables/blk1_brier_cache.RData", envir = e)
pass <- e$pass; settings <- e$settings; n_target <- e$n_target
methods <- c(1L, 2L, 4L); mnm <- c("iid", "AR2", "AR1")
setting_desc <- c(`5` = "near-unit-root", `7` = "ARFIMA d=0.45")
H <- 15L; windows <- list(`1-5` = 1:5, `1-15` = 1:15)
band_probs <- c(0.90, 0.99); G <- 31L

contrasts <- list(c("AR1","iid","AR(1) - iid"),
                  c("AR2","iid","AR(2) - iid"),
                  c("AR2","AR1","AR(2) - AR(1)"))

# per-lead tail-band twCRPS (averaged over sites); zgrid shared per dataset
twcrps_by_lead <- function(yhat, y_val, zgrid) {
  vapply(seq_len(dim(yhat)[3]), function(h) {
    yp_h <- yhat[, , h]; y_h <- y_val[, h]
    acc <- numeric(length(y_h)); prevBS <- NULL
    for (k in seq_along(zgrid)) {
      BS <- (colMeans(yp_h > zgrid[k]) - as.numeric(y_h > zgrid[k]))^2
      if (k > 1) acc <- acc + 0.5 * (BS + prevBS) * (zgrid[k] - zgrid[k - 1])
      prevBS <- BS
    }
    mean(acc, na.rm = TRUE)
  }, numeric(1))
}
paired <- function(d) {
  d <- d[is.finite(d)]
  tt <- t.test(d); w1 <- suppressWarnings(wilcox.test(d, alternative = "less", exact = FALSE))
  list(n = length(d), mean = mean(d), lo = tt$conf.int[1], hi = tt$conf.int[2],
       t_p2 = tt$p.value, w_p1 = w1$p.value)
}

rows <- list()
for (si in seq_along(settings)) {
  st <- settings[si]
  clean <- head(which(apply(pass[, , si], 1, all)), n_target)
  tw <- array(NA_real_, c(H, length(clean), length(methods)))
  for (di in seq_along(clean)) {
    zgrid <- NULL
    for (mi in seq_along(methods)) {
      f <- sprintf("results/blk1-%d-%d-%d.RData", st, methods[mi], clean[di])
      if (!file.exists(f)) { cat("missing:", f, "\n"); next }
      le <- new.env(); load(f, envir = le)
      if (is.null(zgrid)) {                          # band from validation block (all 15 leads)
        band <- quantile(le$res$y_val[, 1:H], probs = band_probs, na.rm = TRUE)
        zgrid <- seq(band[1], band[2], length.out = G)
      }
      tw[, di, mi] <- twcrps_by_lead(le$res$yhat, le$res$y_val, zgrid)
      rm(le); gc()
    }
  }
  for (wn in names(windows)) {
    hs <- windows[[wn]]
    wm <- sapply(seq_along(methods), function(mi) apply(tw[hs, , mi, drop = FALSE], 2, mean))
    colnames(wm) <- mnm
    for (ct in contrasts) {
      r <- paired(wm[, ct[1]] - wm[, ct[2]])
      rows[[length(rows) + 1]] <- data.frame(
        setting = st, desc = setting_desc[as.character(st)], leads = wn,
        contrast = ct[3], n = r$n, Delta = r$mean, CI_lo = r$lo, CI_hi = r$hi,
        t_p_2sided = r$t_p2, wilcox_p_1sided = r$w_p1, stringsAsFactors = FALSE)
    }
  }
  cat(sprintf("setting %d done (clean: %s)\n", st, paste(clean, collapse = ",")))
}
out <- do.call(rbind, rows)
write.csv(out, "output/tables/blk1_paired_ci_twcrps.csv", row.names = FALSE)
cat("Written: output/tables/blk1_paired_ci_twcrps.csv\n")
print(out[out$contrast != "AR(1) - iid", c("setting","leads","contrast","Delta","CI_lo","CI_hi","wilcox_p_1sided")], row.names = FALSE)
