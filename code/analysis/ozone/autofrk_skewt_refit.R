# Improve the skew-t QQ fit for the autoFRK residuals (res2 from
# autofrk_full_qq.R). The original fit standardized residuals by the sample
# sd first, then ran a single-start Nelder-Mead optim -- both steps are
# unstable for heavy-tailed data (nu ~ 4 means the 4th moment barely exists,
# so sample sd is noisy). This version fits skew-t directly on the raw
# residual scale using sn::selm() (Azzalini's dedicated, better-converging
# MLE), and plots theoretical quantiles from that fit directly (no
# standardization detour).
#
# Run: Rscript code/analysis/ozone/autofrk_skewt_refit.R

args <- commandArgs(trailingOnly = FALSE)
file.arg <- grep("^--file=", args, value = TRUE)
if (length(file.arg) > 0) {
  setwd(dirname(normalizePath(sub("^--file=", "", file.arg[1]))))
}

library(sn)

load("autofrk_full_fit.RData") # res2: autoFRK residuals (ns x nt, has NA)

r <- as.vector(res2)
r <- r[!is.na(r)]
n <- length(r)
cat("n obs:", n, "\n")

#### old method: standardize by sample sd, then single-start optim
r.std <- (r - mean(r)) / sd(r)
llike.st <- function(params, y) {
  ll <- dst(y,
    xi = params[1], omega = exp(params[2]),
    alpha = params[3], nu = params[4], log = TRUE
  )
  return(sum(-ll, na.rm = TRUE))
}
old.est <- optim(par = c(0, 1, 0, 3), fn = llike.st, y = r.std)$par
old.dp <- c(xi = old.est[1], omega = exp(old.est[2]), alpha = old.est[3], nu = old.est[4])
cat("old (standardized, single-start optim): xi =", old.dp[1], " omega =", old.dp[2],
    " alpha =", old.dp[3], " nu =", old.dp[4], "\n")

#### new method: sn::selm on raw scale (dedicated skew-t MLE, multiple internal starts)
sel <- selm(r ~ 1, family = "ST")
dp <- coef(sel, "DP")
cat("new (raw scale, selm MLE): xi =", dp["dp1"], " omega =", dp["dp2"],
    " alpha =", dp["dp3"], " nu =", dp["dp4"], "\n")

#### also try a multi-start refit of the old parameterization for comparison
grid <- expand.grid(
  xi = c(-0.5, 0, 0.5), logomega = log(c(0.5, 1, 1.5)),
  alpha = c(-1, 0, 1), nu = c(3, 6, 10)
)
best <- NULL
best.val <- Inf
for (i in seq_len(nrow(grid))) {
  o <- try(
    optim(par = as.numeric(grid[i, ]), fn = llike.st, y = r.std, method = "BFGS"),
    silent = TRUE
  )
  if (!inherits(o, "try-error") && o$value < best.val) {
    best.val <- o$value
    best <- o$par
  }
}
multi.dp <- c(xi = best[1], omega = exp(best[2]), alpha = best[3], nu = best[4])
cat("multi-start BFGS (standardized): xi =", multi.dp[1], " omega =", multi.dp[2],
    " alpha =", multi.dp[3], " nu =", multi.dp[4], " negloglik =", best.val, "\n")

#### QQ comparison: plot raw sorted residuals against each fit's theoretical
#### quantiles, all evaluated on the SAME raw scale (no separate standardization)
p <- (1:n) / (n + 1)
these <- c(1:500, seq(501, n - 500, by = 120), (n - 499):n)
r.sorted <- sort(r)

qq.panel <- function(xplot, main) {
  plot(xplot[these], r.sorted[these],
    xlab = "Theoretical Quantile", ylab = "Observed Quantile",
    main = main, cex = 1.2, cex.lab = 1.3, cex.axis = 1.3, pch = 20
  )
  abline(0, 1)
}

# old fit's quantiles need to be mapped back to raw scale: raw = mean+sd*std
r.mean <- mean(r); r.sd <- sd(r)
old.q <- r.mean + r.sd * qst(p, xi = old.dp[1], omega = old.dp[2], alpha = old.dp[3], nu = old.dp[4])
new.q <- qst(p, xi = dp["dp1"], omega = dp["dp2"], alpha = dp["dp3"], nu = dp["dp4"])

dir.create("plots", showWarnings = FALSE)
pdf("plots/qq-res-autofrk-skewt-refit.pdf", width = 8, height = 4.5)
par(mfrow = c(1, 2), mar = c(5.1, 4.7, 4.1, 2.1))
qq.panel(old.q, sprintf(
  "old fit (std.+single optim)\nalpha=%.2f, nu=%.1f", old.dp[3], old.dp[4]
))
qq.panel(new.q, sprintf(
  "selm MLE (raw scale)\nalpha=%.2f, nu=%.1f", dp["dp3"], dp["dp4"]
))
dev.off()
cat("saved plots/qq-res-autofrk-skewt-refit.pdf\n")
