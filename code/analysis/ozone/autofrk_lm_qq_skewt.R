# QQ plot of the lm(ozone ~ CMAQ) residuals (pre-autoFRK) against Normal,
# t(10), and an ML-fitted skew-t. Reuses `res` saved by autofrk_full_qq.R
# (autofrk_full_fit.RData) so lm/autoFRK do not need to be refit.
#
# Run: Rscript code/analysis/ozone/autofrk_lm_qq_skewt.R

args <- commandArgs(trailingOnly = FALSE)
file.arg <- grep("^--file=", args, value = TRUE)
if (length(file.arg) > 0) {
  setwd(dirname(normalizePath(sub("^--file=", "", file.arg[1]))))
}

library(sn)

load("autofrk_full_fit.RData") # provides res (lm residuals), res2, fit, S

qq.prep <- function(r) {
  r <- sort(as.vector(r))
  n <- length(r)
  list(
    res = (r - mean(r)) / sd(r),
    p = (1:n) / (n + 1),
    these = c(1:500, seq(501, n - 500, by = 120), (n - 499):n)
  )
}

llike.st <- function(params, y) {
  ll <- dst(y,
    xi = params[1], omega = exp(params[2]),
    alpha = params[3], nu = params[4], log = TRUE
  )
  return(sum(-ll, na.rm = TRUE))
}
fit.st <- function(r) {
  est <- optim(par = c(0, 1, 0, 3), fn = llike.st, y = r)$par
  c(xi = est[1], omega = exp(est[2]), alpha = est[3], nu = est[4])
}

qq.lm <- qq.prep(res)
st.lm <- fit.st(qq.lm$res)
cat("skew-t ML on lm residuals: xi =", st.lm[1], " omega =", st.lm[2],
    " alpha =", st.lm[3], " nu =", st.lm[4], "\n")

qq.panel <- function(xplot, qq, main) {
  plot(xplot, qq$res[qq$these],
    xlab = "Theoretical Quantile", ylab = "Observed Quantile",
    main = main, cex = 1.2, cex.lab = 1.3, cex.axis = 1.3, pch = 20
  )
  abline(0, 1)
}

dir.create("plots", showWarnings = FALSE)
pdf("plots/qq-res-lm-skewt.pdf", width = 12, height = 4.5)
par(mfrow = c(1, 3), mar = c(5.1, 4.7, 4.1, 2.1))
qq.panel(qnorm(qq.lm$p[qq.lm$these]), qq.lm, "lm residuals vs Normal")
qq.panel(qt(qq.lm$p[qq.lm$these], 10), qq.lm, "lm residuals vs t(10)")
qq.panel(
  qst(qq.lm$p[qq.lm$these],
    xi = st.lm[1], omega = st.lm[2], alpha = st.lm[3], nu = st.lm[4]
  ),
  qq.lm,
  sprintf(
    "lm residuals vs skew-t (alpha=%.2f, nu=%.1f)",
    st.lm[3], st.lm[4]
  )
)
dev.off()
cat("saved plots/qq-res-lm-skewt.pdf\n")
