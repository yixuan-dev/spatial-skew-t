# Thesis figure: two-panel QQ plot of the lm(ozone ~ CMAQ) residuals against
# Normal and an ML-fitted skew-t. Reuses `res` saved by autofrk_full_qq.R
# (autofrk_full_fit.RData) so lm/autoFRK do not need to be refit.
# Writes directly to myLatex/pdf/ (same precedent as myLatex/pdf/plots_map.R).
#
# Run: Rscript code/analysis/ozone/qq_lm_2panel_thesis.R

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

out <- "../../../myLatex/pdf/qq-res-lm-2panel.pdf"
pdf(out, width = 9, height = 4.5)
par(mfrow = c(1, 2), mar = c(5.1, 4.7, 4.1, 2.1))
qq.panel(qnorm(qq.lm$p[qq.lm$these]), qq.lm, "Normal")
qq.panel(
  qst(qq.lm$p[qq.lm$these],
    xi = st.lm[1], omega = st.lm[2], alpha = st.lm[3], nu = st.lm[4]
  ),
  qq.lm, "Fitted skew-t"
)
dev.off()
cat("saved", out, "\n")
