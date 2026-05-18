# AR(2) generation test

Compare two ways to simulate an AR(2) series in R, fit each with `arima(2,0,0)`, and run residual diagnostics.

## General ways to generate AR(2) data in R

1. **`arima.sim`** — standard, handles burn-in, custom innovations, stationarity check.
2. **Manual recursion** — iterate $X_t = \phi_1 X_{t-1} + \phi_2 X_{t-2} + \varepsilon_t$.
3. **`filter()` with `method = "recursive"`** — vectorised recursive filtering of white noise.
4. **`simulate()` on a fitted `arima` / `Arima` object** — for bootstrap / forecasting studies.
5. **Custom innovations** — via `rand.gen` or `innov` in `arima.sim` for heavy-tailed / non-Gaussian errors.

Key considerations: stationarity of AR roots, burn-in to remove start-up bias, innovation distribution, and `set.seed()` for reproducibility.

## Setup

- True parameters: $\phi = (0.6,\ -0.3)$, $\sigma = 1$
- Length: $n = 1000$
- Seed: 123

## Methods compared

### 1. Manual recursion
```r
eps <- rnorm(n, sd = sig)
x   <- numeric(n)
x[1:2] <- rnorm(2)
for (t in 3:n) x[t] <- phi[1]*x[t-1] + phi[2]*x[t-2] + eps[t]
```

### 2. `filter()` (recursive)
```r
eps <- rnorm(n, sd = sig)
x   <- as.numeric(stats::filter(eps, filter = phi, method = "recursive"))
```

Both series then fit with `arima(x, order = c(2,0,0))`.

## Results

### Estimated coefficients

| Method           | $\hat\phi_1$  | $\hat\phi_2$     | $\hat\sigma^2$ |
| ---------------- | ------------- | ---------------- | -------------- |
| Manual recursion | 0.576 (0.030) | $-0.322$ (0.030) | 0.98           |
| `filter()`       | 0.549 (0.030) | $-0.294$ (0.030) | 1.01           |

Both estimates are within ~1 SE of the true $(0.6, -0.3)$.

### Residual diagnostics

| Test                            | Manual recursion     | `filter()`           |
| ------------------------------- | -------------------- | -------------------- |
| Mean of residuals               | $-4.3\times 10^{-5}$ | $-3.0\times 10^{-5}$ |
| SD of residuals                 | 0.991                | 1.007                |
| Ljung-Box (lag 10, fitdf = 2) p | 0.916                | 0.643                |
| Shapiro-Wilk p                  | 0.430                | 0.626                |

- Ljung-Box: no remaining autocorrelation in residuals for either fit.
- Shapiro-Wilk: residuals consistent with Gaussian innovations.
- ACF / PACF / Q-Q plots (see `ar2_residuals.pdf`) show no structure.

## Conclusion

Manual recursion and `filter(..., method = "recursive")` produce statistically equivalent AR(2) series. `arima(2,0,0)` recovers the true coefficients in both cases and the residuals behave as white Gaussian noise. For most simulation work either method is fine; `filter()` is more concise and vectorised, while manual recursion is more transparent and easier to customise (e.g. nonstandard initial conditions or time-varying coefficients).

## Files

- `ar2_sim_fit.R` — simulation + fitting + diagnostics script
- `ar2_residuals.pdf` — residual plots (time series, ACF, PACF, Q-Q) for both fits
