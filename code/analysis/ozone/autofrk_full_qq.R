# Fit autoFRK to the full ozone data (all sites with <=50% missing days),
# after removing the intercept + CMAQ trend, and check whether the residuals
# are still heavy-tailed via QQ plots (Normal / t(10) / ML-fitted skew-t).
# Companion diagnostic to US-all/us-all-setup.R (qq-res.pdf), which uses the
# 800-site subsample and lm residuals only.
#
# Run: Rscript code/analysis/ozone/autofrk_full_qq.R

# work relative to this script's directory
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

# Exclude if site is missing more than 50% of its days
excl <- which(rowMeans(is.na(Y)) > 0.50)
index <- index[-excl]
Y <- Y[-excl, ]
S <- S[-excl, ]

cmaq <- CMAQ[index, ]
cmaq <- (cmaq - mean(cmaq)) / sd(cmaq)
S <- S / 1000

ns <- nrow(Y)
nt <- ncol(Y)
cat("sites:", ns, "days:", nt, "obs:", sum(!is.na(Y)), "\n")

#### lm detrend: ozone ~ intercept + CMAQ
# stacking by column is equivalent to the day-loop in us-all-setup.R
ozone.lm <- lm(as.vector(Y) ~ as.vector(cmaq))
ozone.beta <- coef(ozone.lm)
cat("lm coefficients (Int, CMAQ):", ozone.beta, "\n")
res <- Y - (ozone.beta[1] + ozone.beta[2] * cmaq)

#### autoFRK fit on the detrended residual matrix
t0 <- proc.time()
fit <- autoFRK(data = res, loc = S)
cat("autoFRK elapsed:", round((proc.time() - t0)[3], 1), "sec\n")
cat("selected K:", ncol(fit$G), "\n")

# in-sample fitted values G %*% w (ns x nt)
pred <- predict(fit)$pred.value
stopifnot(dim(pred) == c(ns, nt))
res2 <- res - pred

save(fit, res, res2, S, file = "autofrk_full_fit.RData")

cat("residual sd  lm:", sd(res, na.rm = TRUE),
    " autoFRK:", sd(res2, na.rm = TRUE), "\n")

#### QQ diagnostics (pattern of us-all-setup.R:320-383)
# sorted standardized residuals + plotting positions, thinned in the middle
qq.prep <- function(r) {
  r <- sort(as.vector(r)) # drops NA
  n <- length(r)
  list(
    res = (r - mean(r)) / sd(r),
    p = (1:n) / (n + 1),
    these = c(1:500, seq(501, n - 500, by = 120), (n - 499):n)
  )
}

# ML skew-t fit (xi location, log-omega scale, alpha slant, nu df)
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
qq.fr <- qq.prep(res2)

st.lm <- fit.st(qq.lm$res)
st.fr <- fit.st(qq.fr$res)
cat("skew-t ML on lm residuals:      xi =", st.lm[1], " omega =", st.lm[2],
    " alpha =", st.lm[3], " nu =", st.lm[4], "\n")
cat("skew-t ML on autoFRK residuals: xi =", st.fr[1], " omega =", st.fr[2],
    " alpha =", st.fr[3], " nu =", st.fr[4], "\n")

qq.panel <- function(xplot, qq, main) {
  plot(xplot, qq$res[qq$these],
    xlab = "Theoretical Quantile", ylab = "Observed Quantile",
    main = main, cex = 1.2, cex.lab = 1.3, cex.axis = 1.3, pch = 20
  )
  abline(0, 1)
}

dir.create("plots", showWarnings = FALSE)
pdf("plots/qq-res-autofrk.pdf", width = 16, height = 4.5)
par(mfrow = c(1, 4), mar = c(5.1, 4.7, 4.1, 2.1))
qq.panel(qnorm(qq.lm$p[qq.lm$these]), qq.lm, "lm residuals vs Normal")
qq.panel(qnorm(qq.fr$p[qq.fr$these]), qq.fr, "autoFRK residuals vs Normal")
qq.panel(qt(qq.fr$p[qq.fr$these], 10), qq.fr, "autoFRK residuals vs t(10)")
qq.panel(
  qst(qq.fr$p[qq.fr$these],
    xi = st.fr[1], omega = st.fr[2], alpha = st.fr[3], nu = st.fr[4]
  ),
  qq.fr,
  sprintf(
    "autoFRK residuals vs skew-t (alpha=%.2f, nu=%.1f)",
    st.fr[3], st.fr[4]
  )
)
dev.off()
cat("saved plots/qq-res-autofrk.pdf\n")
