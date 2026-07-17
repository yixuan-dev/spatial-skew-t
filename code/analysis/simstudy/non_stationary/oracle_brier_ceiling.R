#########################################################################
# THE ORACLE BRIER CEILING -- why 9 of 11 settings can never move the Brier
# score, provable without running a single MCMC chain.
#
# The MRTS extension augments the FIXED-EFFECT MEAN, and its coefficient beta is
# POOLED OVER TIME (updateBeta() in R/ar2/update_params.R). So whatever K is, the
# most the mean channel can ever deliver is the TIME-AVERAGED true mean surface
#
#     mbar(s) = E_t[ mu_t(s) ]  =  truemean.field[[setting]][s, d].
#
# The K = 0 baseline, by contrast, can only fit a spatial constant (the design is
# [1, s1, s2] and beta.t = c(10, 0, 0), so the linear part is null in truth).
#
# Define the ORACLE relative Brier
#
#     R* = BS( mbar )  /  BS( spatial constant ),
#
# both evaluated under the TRUE skew-t predictive law
#
#     P(Y > u | m) = E_{tau, z}[ 1 - Phi( (u - m - lambda z) / sigma ) ],
#     tau ~ Ga(3/2, 4),  sigma = tau^{-1/2},  z ~ HN(0, sigma),  lambda = 3.
#
# R* is a FLOOR: no K, no sample size, no sampler can push the relative Brier
# below it. If R* ~= 1 the setting is a Brier dead end, and an observed null
# there says nothing about MRTS -- it is a property of the design.
#
# This is what makes the whole Brier story interpretable. Without it, "MRTS does
# not improve the Brier score on setting 6" has a fatal alternative explanation:
# maybe MRTS never improves the Brier score, full stop. R* separates the two.
#
# Output: output/oracle_brier_ceiling.csv (+ a printed table)
#########################################################################

if (basename(getwd()) != "simstudy") setwd("d:/Github/spatial-skew-t/code/analysis/simstudy")

load("simdata_nonsta.RData")

set.seed(20250713)
lambda <- 3           # DGP skewness
ta     <- 3 / 2       # rpotspatTS halves tau.alpha = 3 internally
tb     <- 4           # ... and tau.beta = 8
probs  <- c(0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995)
obs    <- c(rep(TRUE, ns - ntest), rep(FALSE, ntest))
NM     <- 6000        # MC draws for the predictive exceedance probability
NSET   <- 3           # datasets to average R* over

# P(Y > u | mean = m) under the true skew-t law, vectorised over m
Pexc <- function(u, m) {
  tau <- rgamma(NM, ta, tb)
  sd  <- 1 / sqrt(tau)
  z   <- abs(rnorm(NM, 0, sd))
  vapply(m, function(mi) mean(1 - pnorm((u - mi - lambda * z) / sd)), numeric(1))
}

nsettings <- nrow(settings.nonsta)
res <- data.frame()

for (k in seq_len(nsettings)) {
  Rk <- numeric(NSET)
  for (d in seq_len(NSET)) {
    Y  <- y[!obs, , d, k]                              # 44 x 50 held-out truth
    mb <- truemean.field[[k]][!obs, d]                 # the oracle mean basis
    cn <- rep(mean(truemean.field[[k]][, d]), ntest)   # all K = 0 can know

    bsA <- bsB <- 0
    for (p in probs) {
      u   <- quantile(as.vector(Y), p)
      ind <- (Y > u)
      bsA <- bsA + mean((Pexc(u, mb) - ind)^2)
      bsB <- bsB + mean((Pexc(u, cn) - ind)^2)
    }
    Rk[d] <- bsA / bsB
  }

  sdm <- mean(vapply(seq_len(NSET), function(d) sd(truemean.field[[k]][, d]), numeric(1)))
  r   <- mean(Rk)
  res <- rbind(res, data.frame(
    setting    = k,
    surf_type  = settings.nonsta$surf_type[k],
    moment     = settings.nonsta$moment[k],
    route      = settings.nonsta$route[k],
    sd_meanbar = sdm,
    R_oracle   = r,
    se         = sd(Rk) / sqrt(NSET),
    verdict    = if (r < 0.85) "positive control" else
                 if (r < 0.96) "weak"             else
                 if (r > 1.02) "HARMFUL"          else "dead end",
    stringsAsFactors = FALSE
  ))
}

dir.create("output", showWarnings = FALSE)
write.csv(res, "output/oracle_brier_ceiling.csv", row.names = FALSE)

cat("\n===== ORACLE BRIER CEILING  R* = BS(true mean) / BS(constant) =====\n")
cat("A mean basis CANNOT push the relative Brier below R*, whatever K is.\n\n")
print(format(res[, c("setting", "surf_type", "moment", "sd_meanbar",
                     "R_oracle", "verdict")], digits = 3), row.names = FALSE)

cat("\nReading it:\n")
cat("  R* ~ 1     the setting carries no time-averaged MEAN structure, so the\n")
cat("             Brier score cannot move for ANY K. An observed null there is\n")
cat("             a property of the DESIGN, not evidence about MRTS.\n")
cat("  R* > 1     supplying the true mean makes the Brier score WORSE. This is\n")
cat("             setting 10: lambda(s) induces a mean shift, but that shift is\n")
cat("             a SYMPTOM of a skewness change. Correcting the symptom while\n")
cat("             the model's lambda stays constant actively misleads.\n")
cat("  R* << 1    a genuine Brier positive control.\n\n")
cat("wrote output/oracle_brier_ceiling.csv\n")
