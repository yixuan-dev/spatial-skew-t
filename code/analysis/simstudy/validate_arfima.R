#########################################################################
# validate_arfima.R  --  de-risk gate for the long-memory (ARFIMA) DGP.
#
# Pure theory checks, NO MCMC. Run this BEFORE regenerating simdata or
# fitting anything. If any check fails, stop and fix simulate_arfima_standard
# before spending compute on the full study.
#
#   $R = "C:\Program Files\R\R-4.5.1\bin\Rscript.exe"
#   & $R validate_arfima.R
#
# IMPORTANT: the copula precondition is that the process is marginally N(0,1)
# ACROSS THE ENSEMBLE (probability space), i.e. Var(x_t) = 1 at each fixed t.
# A SINGLE long-memory realization has a non-zero sample mean and a sample
# variance below 1 -- a known finite-sample effect of the slowly-decaying
# low-frequency power, NOT a generator bug. So all marginal/ACF checks below
# use ensemble moments over many independent realizations (the columns of one
# simulate_arfima_standard call with n = R rows).
#
# Checks (plan Part E):
#   1. Ensemble marginal is N(0,1): E[x_t] ~ 0, Var(x_t) ~ 1 at every t.
#   2. Ensemble ACF matches theoretical ARFIMA rho(h).
#   3. Ensemble-mean arima(1,0,0)/(2,0,0) recover the Yule-Walker projections.
#   4. Davies-Harte circulant eigenvalues are non-negative across the d grid
#      and both simulation lengths (nt = 50 study, nt = 5000 validation).
#########################################################################

source("../../R/ar2/auxfunctions_ar2.R")

set.seed(2026)
ok <- TRUE
report <- function(pass, msg) {
  cat(sprintf("[%s] %s\n", if (pass) " OK " else "FAIL", msg))
  if (!pass) ok <<- FALSE
}
approx_eq <- function(a, b, tol) all(abs(a - b) <= tol)

d_grid   <- c(0.20, 0.35, 0.45)
nt_study <- 50L      # Morris setup.R length

cat("=== ARFIMA(0,d,0) generator validation ===\n\n")

## ---- pre-registered Yule-Walker projections (the expected-result table) ----
cat("Yule-Walker pseudo-true projections (analysis-model targets):\n")
cat(sprintf("  %4s %6s %9s %18s\n", "d", "Hurst", "AR(1) phi", "AR(2) (phi1, phi2)"))
for (d in c(d_grid, 0.40)) {
  yw <- arfima_yw_projection(d)
  cat(sprintf("  %4.2f %6.2f %9.3f     (%.3f, %.3f)\n",
              yw$d, yw$hurst, yw$ar1, yw$ar2["phi1"], yw$ar2["phi2"]))
}
cat("\n")

## ---- Check 4: eigenvalue non-negativity across grid x lengths ------------
cat("--- Check 4: Davies-Harte embedding is PSD across d grid x length ---\n")
for (d in d_grid) for (nt in c(nt_study, 5000L)) {
  res <- tryCatch({ simulate_arfima_standard(nt = nt, n = 1L, d = d); TRUE },
                  error = function(e) { cat("   ", conditionMessage(e), "\n"); FALSE })
  report(res, sprintf("embeddable: d = %.2f, nt = %d", d, nt))
}
cat("\n")

## ---- Check 1: ensemble marginal N(0,1) at study length -------------------
cat("--- Check 1: ensemble marginal E[x_t]~0, Var(x_t)~1 (d=0.35, R=8000, nt=50) ---\n")
R <- 8000L
X <- simulate_arfima_standard(nt = nt_study, n = R, d = 0.35)   # R x nt
mean_t <- colMeans(X)
var_t  <- apply(X, 2, var)
cat(sprintf("   E[x_t]:   mean over t = % .4f   (max |.| = %.4f)\n", mean(mean_t), max(abs(mean_t))))
cat(sprintf("   Var(x_t): mean over t = % .4f   (range [%.3f, %.3f])\n", mean(var_t), min(var_t), max(var_t)))
report(max(abs(mean_t)) < 0.05, "ensemble mean ~ 0 at every t")
report(approx_eq(var_t, 1, 0.06), "ensemble variance ~ 1 at every t (copula precondition)")

# marginal normality of the pooled ensemble at a fixed t
report(suppressWarnings(shapiro.test(X[1:4000, 25])$p.value) > 0.01,
       "x_t is Gaussian at fixed t (Shapiro on ensemble)")
cat("\n")

## ---- Check 2: ensemble ACF matches theory --------------------------------
cat("--- Check 2: ensemble ACF vs theoretical ARFIMA rho(h), d=0.35 ---\n")
rho_theo <- arfima_acf(0.35, 10)[-1]
# lag-h autocorrelation from ensemble cross-moments (variance ~ 1)
rho_emp <- sapply(1:10, function(h) mean(rowMeans(X[, 1:(nt_study - h)] * X[, (1 + h):nt_study])))
cat("   lag :  ", sprintf("%6d", 1:10), "\n")
cat("   theo:  ", sprintf("%6.3f", rho_theo), "\n")
cat("   emp :  ", sprintf("%6.3f", rho_emp), "\n")
report(max(abs(rho_emp - rho_theo)) < 0.03,
       sprintf("ensemble ACF matches theory (max err %.4f)", max(abs(rho_emp - rho_theo))))
cat("\n")

## ---- Check 3: generator's 2nd-order structure yields the YW projection ---
# Validates the GENERATOR (not an estimator): from ensemble autocovariances
# gamma_hat(0..2) (huge sample ~ population) solve the Yule-Walker equations
# and confirm they match arfima_yw_projection(). arima() is NOT used here: it
# targets the ML/spectral projection, a DIFFERENT estimand than population YW,
# so it would not match even with a perfect generator.
cat("--- Check 3: ensemble Yule-Walker solve matches projection (d=0.40) ---\n")
d <- 0.40
yw <- arfima_yw_projection(d)
Xg <- simulate_arfima_standard(nt = 400L, n = 6000L, d = d)   # ~population moments
g0 <- mean(Xg^2)
g1 <- mean(rowMeans(Xg[, 1:399] * Xg[, 2:400]))
g2 <- mean(rowMeans(Xg[, 1:398] * Xg[, 3:400]))
ar1_yw <- g1 / g0
den <- g0^2 - g1^2
ar2_yw <- c(g1 * (g0 - g2) / den, (g0 * g2 - g1^2) / den)
cat(sprintf("   AR(1): ensemble-YW = %.3f            (target %.3f)\n", ar1_yw, yw$ar1))
cat(sprintf("   AR(2): ensemble-YW = (%.3f, %.3f)   (target (%.3f, %.3f))\n",
            ar2_yw[1], ar2_yw[2], yw$ar2["phi1"], yw$ar2["phi2"]))
report(approx_eq(ar1_yw, yw$ar1, 0.02), "AR(1) ensemble-YW ~ projection")
report(approx_eq(ar2_yw[1], yw$ar2["phi1"], 0.02) && approx_eq(ar2_yw[2], yw$ar2["phi2"], 0.02),
       "AR(2) ensemble-YW ~ projection (positive phi2)")

# Informational (NOT gating): what a finite-N=50 ML fit -- like the backend
# posterior -- actually recovers. Expect phi2 > 0 but biased BELOW the YW
# projection: low-order AR MLE underestimates persistence for long memory at
# short N. The study's phi2-vs-YW plot should expect posteriors under the YW
# line, with the QUALITATIVE claims (phi2 > 0, increasing in d) intact.
Xs <- simulate_arfima_standard(nt = nt_study, n = 300L, d = d)
ml2 <- t(apply(Xs, 1, function(x) tryCatch(
  arima(x, c(2, 0, 0), include.mean = TRUE)$coef[c("ar1", "ar2")],
  error = function(e) c(NA, NA))))
cat(sprintf("   [info] finite N=50 ML AR(2) mean = (%.3f, %.3f); phi2>0 in %.0f%% -- expect posterior near here, below YW\n",
            mean(ml2[, 1], na.rm = TRUE), mean(ml2[, 2], na.rm = TRUE),
            100 * mean(ml2[, 2] > 0, na.rm = TRUE)))
cat("\n")

cat(sprintf("=== RESULT: %s ===\n", if (ok) "ALL CHECKS PASSED" else "FAILURES PRESENT -- do not proceed"))
if (!ok) quit(status = 1)
