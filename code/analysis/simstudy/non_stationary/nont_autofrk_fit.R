#########################################################################
# autoFRK benchmark on simdata_nonsta.RData, ALL settings 1-21.
#
# For each (setting, dataset): fit autoFRK::autoFRK() on the 100 training
# sites (data grand-mean centred; mrts basis, K selected by AIC up to the
# default maxK = 100), predict the 44 held-out sites at every t, and report
#
#   mse      = mean((yhat - y_test)^2) over 44 x 50        (predictive MSE)
#   rel      = mse / mse0, mse0 = per-day spatial-constant predictor
#              (each day predicted by that day's training mean)
#   K        = number of mrts basis functions selected by AIC
#   rec.rmse = RMSE of the time-averaged prediction against the true mean
#              surface at the test sites (truemean.field contract)
#
# autoFRK is the pure-MRTS machine (mean + low-rank random effects, no
# skew, no heavy tails), so this table is the "how far does MRTS alone go"
# reference for every DGP, hours of MCMC not required.
#
# Output: output/nont_autofrk_fit.csv + a per-setting summary table.
#########################################################################

if (basename(getwd()) != "simstudy") setwd("d:/Github/spatial-skew-t/code/analysis/simstudy")

suppressMessages(library(autoFRK))

load("simdata_nonsta.RData")

datasets <- 1:10
obs      <- c(rep(TRUE, ns - ntest), rep(FALSE, ntest))
s.o      <- s[obs, , drop = FALSE]
s.p      <- s[!obs, , drop = FALSE]

res <- vector("list", 0)
tic.all <- proc.time()

for (setting in seq_len(dim(y)[4])) {
  for (set in datasets) {
    y.d  <- y[, , set, setting]
    y.o  <- y.d[obs, ]
    y.p  <- y.d[!obs, ]

    mu0  <- mean(y.o)
    fit  <- autoFRK(data = y.o - mu0, loc = s.o)
    yhat <- predict(fit, newloc = s.p)$pred.value + mu0

    # per-day spatial-constant baseline (each day's training mean)
    y0   <- matrix(colMeans(y.o), ntest, nt, byrow = TRUE)

    mse  <- mean((yhat - y.p)^2)
    mse0 <- mean((y0   - y.p)^2)

    truemean <- truemean.field[[setting]][!obs, set]
    rec.rmse <- sqrt(mean((rowMeans(yhat) - truemean)^2))

    res[[length(res) + 1]] <- data.frame(
      setting = setting, surf_type = settings.nonsta$surf_type[setting],
      set = set, K = ncol(fit$G), mse = mse, rel = mse / mse0,
      rec.rmse = rec.rmse
    )
  }
  cat(sprintf("setting %2d (%-16s) done  [%.1f min elapsed]\n",
              setting, settings.nonsta$surf_type[setting],
              (proc.time() - tic.all)[3] / 60))
}

res <- do.call(rbind, res)
write.csv(res, "output/nont_autofrk_fit.csv", row.names = FALSE)

se <- function(v) sd(v) / sqrt(length(v))
summ <- do.call(rbind, lapply(split(res, res$setting), function(d) data.frame(
  setting   = d$setting[1],
  surf_type = d$surf_type[1],
  K.mean    = mean(d$K),
  mse       = mean(d$mse),
  mse.se    = se(d$mse),
  rel       = mean(d$rel),
  rec.rmse  = mean(d$rec.rmse)
)))

cat(sprintf("\n=== autoFRK benchmark, datasets %s (n = %d per setting) ===\n",
            paste(range(datasets), collapse = "-"), length(datasets)))
print(format(summ, digits = 3), row.names = FALSE)
cat("\nwrote output/nont_autofrk_fit.csv\n")
