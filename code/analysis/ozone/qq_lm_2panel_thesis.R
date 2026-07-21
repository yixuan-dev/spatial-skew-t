# Thesis figure: two-panel QQ plot of the lm(ozone ~ CMAQ) residuals against
# Normal and the skew-t with lambda = 1, a = 10 (thesis/STP notation; sn's
# direct parameters alpha = lambda, nu = a) used for the same diagnostic in
# Morris et al. (2017). Reuses `res` saved by autofrk_full_qq.R
# (autofrk_full_fit.RData) so lm/autoFRK do not need to be refit.
# Writes directly to myLatex/pdf/ (same precedent as myLatex/pdf/plots_map.R).
#
# Both axes are on the unit-variance scale: the residuals are standardized by
# their sample mean/sd, and the skew-t quantiles by the analytic mean/sd of
# st(0, 1, alpha = 1, nu = 10). (The Morris et al. original centered the
# skew-t quantiles but left their scale at omega = 1.)
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

# analytic mean/sd of st(xi = 0, omega = 1, alpha, nu)
st.moments <- function(alpha, nu) {
  delta <- alpha / sqrt(1 + alpha^2)
  b.nu <- sqrt(nu / pi) * exp(lgamma((nu - 1) / 2) - lgamma(nu / 2))
  m <- b.nu * delta
  c(mean = m, sd = sqrt(nu / (nu - 2) - m^2))
}

qq.lm <- qq.prep(res)
st.mom <- st.moments(alpha = 1, nu = 10)

qq.panel <- function(xplot, qq, main) {
  plot(xplot, qq$res[qq$these],
    xlab = "Theoretical Quantile", ylab = "Observed Quantile",
    main = main, cex = 1.2, cex.lab = 1.3, cex.axis = 1.3, pch = 20
  )
  abline(0, 1)
}

out <- "../../../myLatex/pdf/qq-res-lm-2panel.pdf"
pdf(out, width = 9, height = 4.5)
par(mfrow = c(1, 2), mar = c(5.1, 4.7, 4.1, 2.1), font.main = 1)
qq.panel(qnorm(qq.lm$p[qq.lm$these]), qq.lm, "Normal")
qq.panel(
  (qst(qq.lm$p[qq.lm$these], alpha = 1, nu = 10) - st.mom["mean"]) / st.mom["sd"],
  qq.lm,
  expression(paste("Skew-", italic(t), " (", lambda == 1, ", ", italic(a) == 10, ")"))
)
dev.off()
cat("saved", out, "\n")
