setwd("D:/Github/spatial-skew-t/code/analysis/simstudy")
source("../../R/auxfunctions.R")
load("simdata.RData")
# load("./results/5-1-1.RData")
# load("./results/5-2-1.RData")
# load("./results/5-3-1.RData")
# load("./results/5-4-1.RData")
# load("./results/5-5-1.RData")
# load("./results/5-6-1.RData")

set <- 1
setting <- 5
filename <- paste("new_scores", setting, ".RData", sep = "")

obs <- c(rep(T, 100), rep(F, 44))
validate <- y[!obs, , set, setting]
probs <- c(0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995)
thresholds <- quantile(y[, , set, setting], probs = probs, na.rm = T)
brier.score <- array(NA, dim = c(length(probs), 6))

fit <- fit.1
pred <- fit$yp
pred <- fit$yp[10001:20000, , ]
brier.score[, 5] <- BrierScore(pred, thresholds, validate)
brier.score[, 6] <- BrierScore(pred, thresholds, validate, trans = TRUE)
brier.score.mean <- brier.score
bs.mean.ref.gau <- array(NA, dim = c(11, 5))

for (j in 1:5) {
    bs.mean.ref.gau[, j] <- brier.score.mean[, (j + 1)] / brier.score.mean[, 1]
}

ymax <- 3.4
ymin <- 0.3
setting.title <- c(
    "Data: Gaussian", "Data: Symmetric-t (K = 1)",
    "Data: Symmetric-t (K = 5)",
    bquote(paste("Data: Skew-t (K = 1, ", lambda == 3, ")")),
    bquote(paste("Data: Skew-t (K = 5, ", lambda == 3, ")")),
    "Data: Max-stable, Reich and Shaby", "Data: transform below T",
    "Data: Max-stable, Brown-Resnick"
)
methods <- c(
    "Skew-t, K = 1, T = q(0.0)", "Sym-t, K = 1, T = q(0.8)",
    "Skew-t, K = 5, T = q(0.0)", "Sym-t, K = 5, T = q(0.8)",
    "Max-stable, T = q(0.80)"
)
bg <- c("firebrick1", "dodgerblue1", "firebrick1", "dodgerblue1", "gray70")
col <- c("firebrick4", "dodgerblue4", "firebrick4", "dodgerblue4", "gray14")
pch <- c(22, 22, 22, 22, 21)
lty <- c(1, 1, 3, 3, 3)

windows(width = 18, height = 12)
plot(probs, bs.mean.ref.gau[, 1],
    type = "o",
    lty = lty[1], pch = pch[1], col = col[1], bg = bg[1], cex = 1.5,
    ylim = c(ymin, ymax),
    main = as.expression(setting.title[5]),
    ylab = "Relative Brier score", xlab = "Threshold quantile", cex.lab = 2,
    cex.axis = 2, cex.main = 2
)

for (i in 2:5) {
    lines(probs, bs.mean.ref.gau[, i], lty = lty[i], col = col[i])
    points(probs, bs.mean.ref.gau[, i],
        pch = pch[i], col = col[i],
        bg = bg[i], cex = 1.5
    )
    abline(h = 1, lty = 2)
}

legend("center",
    legend = methods, lty = lty, col = col, pch = pch, pt.bg = bg,
    y.intersp = 2, bty = "n", cex = 1.9, lwd = 1.5
)

dev.off()
