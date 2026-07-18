# Fit autoFRK directly to raw ozone Y (no lm/CMAQ detrending first) on the
# full ozone data (all sites with <=50% missing days), and check whether the
# residuals are still heavy-tailed via QQ plots (Normal / t(10) / ML skew-t).
# Companion to autofrk_full_qq.R, which detrends with lm(Y ~ CMAQ) first.
#
# Run: Rscript code/analysis/ozone/autofrk_full_qq_direct.R

args <- commandArgs(trailingOnly = FALSE)
file.arg <- grep("^--file=", args, value = TRUE)
if (length(file.arg) > 0) {
  setwd(dirname(normalizePath(sub("^--file=", "", file.arg[1]))))
}

library(autoFRK)
library(sn)

#### data preprocessing (same as us-all-setup.R but without the 800-site subsample)
load("ozone_data.RData")
S <- cbind(x[s[, 1]], y[s[, 2]])

excl <- which(rowMeans(is.na(Y)) > 0.50)
Y <- Y[-excl, ]
S <- S[-excl, ]
S <- S / 1000

ns <- nrow(Y)
nt <- ncol(Y)
cat("sites:", ns, "days:", nt, "obs:", sum(!is.na(Y)), "\n")

#### autoFRK fit directly on raw Y (mean not removed -- autoFRK's mu = 0 default,
#### so the mean surface itself is absorbed into the MRTS fit)
t0 <- proc.time()
fit <- autoFRK(data = Y, loc = S)
cat("autoFRK elapsed:", round((proc.time() - t0)[3], 1), "sec\n")
cat("selected K:", ncol(fit$G), "\n")

pred <- predict(fit)$pred.value
stopifnot(dim(pred) == c(ns, nt))
res <- Y - pred

save(fit, res, S, file = "autofrk_full_fit_direct.RData")

cat("residual sd (raw Y):", sd(as.vector(Y), na.rm = TRUE),
    " autoFRK residual sd:", sd(res, na.rm = TRUE), "\n")

#### QQ diagnostics (same pattern as autofrk_full_qq.R)
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

qq.fr <- qq.prep(res)
st.fr <- fit.st(qq.fr$res)
cat("skew-t ML on direct-autoFRK residuals: xi =", st.fr[1], " omega =", st.fr[2],
    " alpha =", st.fr[3], " nu =", st.fr[4], "\n")

qq.panel <- function(xplot, qq, main) {
  plot(xplot, qq$res[qq$these],
    xlab = "Theoretical Quantile", ylab = "Observed Quantile",
    main = main, cex = 1.2, cex.lab = 1.3, cex.axis = 1.3, pch = 20
  )
  abline(0, 1)
}

dir.create("plots", showWarnings = FALSE)
pdf("plots/qq-res-autofrk-direct.pdf", width = 12, height = 4.5)
par(mfrow = c(1, 3), mar = c(5.1, 4.7, 4.1, 2.1))
qq.panel(qnorm(qq.fr$p[qq.fr$these]), qq.fr, "direct autoFRK residuals vs Normal")
qq.panel(qt(qq.fr$p[qq.fr$these], 10), qq.fr, "direct autoFRK residuals vs t(10)")
qq.panel(
  qst(qq.fr$p[qq.fr$these],
    xi = st.fr[1], omega = st.fr[2], alpha = st.fr[3], nu = st.fr[4]
  ),
  qq.fr,
  sprintf(
    "direct autoFRK residuals vs skew-t (alpha=%.2f, nu=%.1f)",
    st.fr[3], st.fr[4]
  )
)
dev.off()
cat("saved plots/qq-res-autofrk-direct.pdf\n")
