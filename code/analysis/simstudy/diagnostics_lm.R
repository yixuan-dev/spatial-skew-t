#########################################################################
# diagnostics_lm.R  --  long-memory diagnostic figures (plan Part D).
#
# Two figures for the ARFIMA settings 17-19:
#
#  1. IMPLIED-ACF "money plot" (pure theory): for each d, overlay the true
#     ARFIMA(0,d,0) autocorrelation against the ACFs implied by the AR(1) and
#     AR(2) Yule-Walker projections. Shows WHY AR(2) is less wrong than AR(1):
#     AR(1)'s geometric ACF falls below the true hyperbolic decay by lag ~3-5,
#     while AR(2) tracks it out to ~8-10. This is the mechanism plot behind any
#     score difference and needs no fits.
#
#  2. phi2-vs-YW: posterior-mean AR(2) phi2 (method 7) per dataset against the
#     positive Yule-Walker reference, read from output/tables/posterior_
#     summary<setting>.csv. Confirms phi2 > 0 (the signature AR(1) cannot
#     represent). NB at nt=50 the posterior sits BELOW the YW line (finite-
#     sample AR bias, see validate_arfima.R); the qualitative claim
#     (phi2 > 0, increasing in d) is what holds.
#
#   & $R diagnostics_lm.R
#########################################################################

source("../../R/ar2/auxfunctions_ar2.R")   # arfima_acf, arfima_yw_projection

lm_settings <- c("17" = 0.20, "18" = 0.35, "19" = 0.45)
Hmax <- 20L
outdir <- "output/plots/diagnostics"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ACF implied by a stationary AR(2) with coefficients (phi1, phi2).
ar2_implied_acf <- function(phi1, phi2, H) {
  r <- numeric(H + 1L)
  r[1L] <- 1
  r[2L] <- phi1 / (1 - phi2)
  for (h in 2L:H) r[h + 1L] <- phi1 * r[h] + phi2 * r[h - 1L]
  r
}

## ---- Figure 1: implied-ACF money plot ------------------------------------
pdf(file.path(outdir, "implied_acf_money_plot.pdf"), width = 11, height = 4)
op <- par(mfrow = c(1, length(lm_settings)), mar = c(4.2, 4.2, 3, 1), cex.lab = 1.1)
h <- 0:Hmax
for (nm in names(lm_settings)) {
  d <- lm_settings[[nm]]
  yw <- arfima_yw_projection(d)
  rho_true <- arfima_acf(d, Hmax)
  rho_ar1  <- yw$ar1^h
  rho_ar2  <- ar2_implied_acf(yw$ar2["phi1"], yw$ar2["phi2"], Hmax)

  plot(h, rho_true, type = "b", pch = 16, lwd = 2.4, col = "black",
       ylim = c(0, 1), xlab = "lag h", ylab = expression(rho(h)),
       main = sprintf("setting %s: ARFIMA d = %.2f (H = %.2f)", nm, d, yw$hurst))
  lines(h, rho_ar2, type = "b", pch = 17, lwd = 1.8, col = "#1f77b4")
  lines(h, rho_ar1, type = "b", pch = 15, lwd = 1.8, col = "#d62728")
  abline(h = 0.05, lty = 3, col = "grey55")
  legend("topright", bty = "n", lwd = c(2.4, 1.8, 1.8), pch = c(16, 17, 15),
         col = c("black", "#1f77b4", "#d62728"),
         legend = c("true ARFIMA", "AR(2) YW", "AR(1) YW"))
}
par(op)
dev.off()
cat(sprintf("wrote %s\n", file.path(outdir, "implied_acf_money_plot.pdf")))

## ---- Figure 2: phi2-vs-YW (per available posterior summary) ---------------
have <- character(0)
rows <- list()
for (nm in names(lm_settings)) {
  f <- sprintf("output/tables/posterior_summary%s.csv", nm)
  if (!file.exists(f)) next
  csv <- read.csv(f, stringsAsFactors = FALSE)
  sub <- csv[csv$method == 7 & csv$parameter %in% c("phi2.tau", "phi2.z"), ]
  if (nrow(sub) == 0) next
  sub$d <- lm_settings[[nm]]
  sub$setting <- as.integer(nm)
  rows[[nm]] <- sub
  have <- c(have, nm)
}

if (length(rows) > 0) {
  dat <- do.call(rbind, rows)
  pdf(file.path(outdir, "phi2_vs_yw.pdf"), width = 7, height = 5)
  par(mar = c(4.5, 4.5, 3, 1))
  ds <- sort(unique(dat$d))
  yw2 <- sapply(ds, function(dd) arfima_yw_projection(dd)$ar2["phi2"])
  ylim <- range(0, dat$mean_post_mean, yw2, na.rm = TRUE) + c(-0.02, 0.05)
  plot(NA, xlim = range(ds) + c(-0.03, 0.03), ylim = ylim,
       xlab = "ARFIMA d", ylab = expression(hat(phi)[2]),
       main = "AR(2) posterior phi2 vs Yule-Walker projection")
  abline(h = 0, col = "grey70")
  lines(ds, yw2, type = "b", pch = 4, lwd = 2, col = "grey30")
  pts <- dat[dat$parameter == "phi2.tau", ]
  points(jitter(pts$d, amount = 0.004), pts$mean_post_mean, pch = 16, col = "#1f77b4")
  ptz <- dat[dat$parameter == "phi2.z", ]
  points(jitter(ptz$d, amount = 0.004), ptz$mean_post_mean, pch = 1, col = "#ff7f0e")
  legend("topleft", bty = "n", pch = c(4, 16, 1), lwd = c(2, NA, NA),
         col = c("grey30", "#1f77b4", "#ff7f0e"),
         legend = c("YW projection (target)", "posterior phi2.tau", "posterior phi2.z"))
  dev.off()
  cat(sprintf("wrote %s (settings: %s)\n",
              file.path(outdir, "phi2_vs_yw.pdf"), paste(have, collapse = ", ")))
} else {
  cat("phi2-vs-YW: no posterior_summary<setting>.csv found yet; run posterior.R first.\n")
}
