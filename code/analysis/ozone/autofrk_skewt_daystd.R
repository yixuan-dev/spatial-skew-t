# Improve the skew-t QQ fit for autoFRK residuals (res2 from
# autofrk_full_qq.R) by diagnosing WHY the global fit's tail was loose.
#
# Finding: 2005-07-09 is a genuine ozone episode day (raw Y up to 153.5 ppb),
# with autoFRK residual sd ~27 vs a typical day's ~7-9 -- i.e. day-to-day
# variance is highly heteroscedastic, not a data error. Pooling all days
# with a single global scale mixes this heteroscedasticity into the
# "heavy tail", which a single skew-t then has to compromise on: mediocre
# fit everywhere instead of a clean bulk fit for typical days.
#
# Fix: standardize each day's residual column by its OWN day sd before
# pooling, isolating the cross-sectional shape (skew/kurtosis within a
# day) from the day-to-day scale differences. Refit skew-t on the
# day-standardized residuals and compare QQ fit to the original global fit.
#
# Run: Rscript code/analysis/ozone/autofrk_skewt_daystd.R

args <- commandArgs(trailingOnly = FALSE)
file.arg <- grep("^--file=", args, value = TRUE)
if (length(file.arg) > 0) {
  setwd(dirname(normalizePath(sub("^--file=", "", file.arg[1]))))
}

library(sn)

load("autofrk_full_fit.RData") # res2: autoFRK residuals (ns x nt), has NA

day.sd <- apply(res2, 2, sd, na.rm = TRUE)
cat("per-day residual sd range:", range(day.sd), "\n")
cat("worst day:", names(which.max(day.sd)), " sd =", max(day.sd), "\n")

res2.daystd <- sweep(res2, 2, day.sd, "/")

fit.skewt <- function(r, label) {
  r <- r[!is.na(r)]
  sel <- selm(r ~ 1, family = "ST")
  dp <- coef(sel, "DP")
  cat(label, ": xi =", dp["xi"], " omega =", dp["omega"],
      " alpha =", dp["alpha"], " nu =", dp["nu"], "\n")
  list(r = sort(r), dp = dp)
}

global <- fit.skewt(as.vector(res2), "global (pooled, no day-scaling)")
daystd <- fit.skewt(as.vector(res2.daystd), "day-standardized")

qq.compare <- function(fitobj, main) {
  r <- fitobj$r
  n <- length(r)
  p <- (1:n) / (n + 1)
  these <- c(1:500, seq(501, n - 500, by = 120), (n - 499):n)
  dp <- fitobj$dp
  xplot <- qst(p[these], xi = dp["xi"], omega = dp["omega"], alpha = dp["alpha"], nu = dp["nu"])
  plot(xplot, r[these],
    xlab = "Theoretical Quantile", ylab = "Observed Quantile",
    main = main, cex = 1.2, cex.lab = 1.3, cex.axis = 1.3, pch = 20
  )
  abline(0, 1)
}

dir.create("plots", showWarnings = FALSE)
pdf("plots/qq-res-autofrk-daystd.pdf", width = 9, height = 4.5)
par(mfrow = c(1, 2), mar = c(5.1, 4.7, 4.1, 2.1))
qq.compare(global, sprintf(
  "pooled (no day-scaling)\nalpha=%.2f, nu=%.1f", global$dp["alpha"], global$dp["nu"]
))
qq.compare(daystd, sprintf(
  "day-standardized\nalpha=%.2f, nu=%.1f", daystd$dp["alpha"], daystd$dp["nu"]
))
dev.off()
cat("saved plots/qq-res-autofrk-daystd.pdf\n")
