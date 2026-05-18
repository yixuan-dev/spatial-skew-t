## AR(2) simulation: manual recursion vs filter, then fit arima(2,0,0)
## and run residual diagnostics on each.

set.seed(123)

n   <- 1000          # series length
phi <- c(0.6, -0.3)  # AR(2) coefficients
sig <- 1             # innovation sd

## --- 1) Manual recursion ----------------------------------------------------
eps_m   <- rnorm(n, sd = sig)
x_man   <- numeric(n)
x_man[1:2] <- rnorm(2)
for (t in 3:n) {
  x_man[t] <- phi[1] * x_man[t - 1] + phi[2] * x_man[t - 2] + eps_m[t]
}

## --- 2) filter() with recursive method --------------------------------------
eps_f <- rnorm(n, sd = sig)
x_flt <- as.numeric(stats::filter(eps_f, filter = phi, method = "recursive"))

## --- Fit arima(2,0,0) to each -----------------------------------------------
fit_man <- arima(x_man, order = c(2, 0, 0))
fit_flt <- arima(x_flt, order = c(2, 0, 0))

cat("\n=== Fit from manual recursion ===\n"); print(fit_man)
cat("\n=== Fit from filter()         ===\n"); print(fit_flt)

## --- Residual analysis ------------------------------------------------------
resid_diag <- function(fit, label) {
  r <- residuals(fit)

  cat("\n--- Residual diagnostics:", label, "---\n")
  cat("mean(resid) =", mean(r), "  sd(resid) =", sd(r), "\n")

  # Ljung-Box test for autocorrelation in residuals
  lb <- Box.test(r, lag = 10, type = "Ljung-Box",
                 fitdf = length(fit$coef) - 1)  # subtract intercept df
  print(lb)

  # Shapiro-Wilk for normality (use a subsample if n > 5000)
  sw <- shapiro.test(if (length(r) > 5000) sample(r, 5000) else r)
  print(sw)

  # Plots
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  plot(r, type = "l", main = paste("Residuals -", label),
       ylab = "resid", xlab = "t"); abline(h = 0, col = "red", lty = 2)
  acf(r,  main = paste("ACF -",  label))
  pacf(r, main = paste("PACF -", label))
  qqnorm(r, main = paste("Normal Q-Q -", label)); qqline(r, col = "red")
  par(op)

  invisible(list(ljung_box = lb, shapiro = sw))
}

pdf("code/ar2_residuals.pdf", width = 8, height = 7)
diag_man <- resid_diag(fit_man, "manual recursion")
diag_flt <- resid_diag(fit_flt, "filter()")
dev.off()

cat("\nDiagnostic plots written to code/ar2_residuals.pdf\n")
