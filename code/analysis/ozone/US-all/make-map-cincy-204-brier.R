# Site-level Brier score map for the Cincinnati window, setting 204.
#
# For each AQS station inside the prediction window, score the model's daily
# exceedance probabilities (taken from the nearest grid cell of the posterior
# predictive sample) against the observed exceedance indicators:
#   Brier_s = mean_t ( P[Y_pred > 75] - 1{Y_obs(s, t) > 75} )^2
# over the station's non-missing days.
#
# Notes on interpretation:
#   - In-sample: these stations were part of the full-data fit (bar the
#     handful excluded for >50% missingness), so this is a calibration
#     diagnostic, not an out-of-sample skill estimate.
#   - Nearest grid cell is at most ~8.5 km from the station (12 km CMAQ grid).
#
# Usage:  Rscript make-map-cincy-204-brier.R [dev]

args <- commandArgs(trailingOnly = TRUE)
suffix <- if (length(args) > 0 && args[1] == "dev") "-dev" else ""

library(fields)

load(paste0("predict-maps-cincy-204", suffix, ".RData"))

# stations: the saved image keeps x, y in raw km (only S.p is /1000-scaled),
# so rescale here; guard in case a future rerun saves them already scaled
if (max(abs(x)) > 100) { x <- x / 1000; y <- y / 1000 }
S.sites <- cbind(x[s[, 1]], y[s[, 2]])

inwin <- S.sites[, 1] >= min(S.p[, 1]) & S.sites[, 1] <= max(S.p[, 1]) &
         S.sites[, 2] >= min(S.p[, 2]) & S.sites[, 2] <= max(S.p[, 2])
S.win <- S.sites[inwin, , drop = FALSE]
Y.win <- Y[inwin, , drop = FALSE]
cat("AQS sites in window:", nrow(S.win), "\n")

# nearest prediction grid cell for each station
g <- apply(rdist(S.win, S.p), 1, which.min)
offset <- rdist(S.win, S.p)[cbind(seq_len(nrow(S.win)), g)]
cat("max station-to-cell offset (km):", round(max(offset) * 1000, 1), "\n")

p.exceed <- 1 - p.below                  # [grid, day]
brier <- rep(NA_real_, nrow(S.win))
ndays <- rep(0L, nrow(S.win))
for (i in seq_len(nrow(S.win))) {
  obs <- Y.win[i, ] > threshold
  ok  <- !is.na(obs)
  ndays[i] <- sum(ok)
  if (ndays[i] > 0) {
    brier[i] <- mean((p.exceed[g[i], ok] - as.numeric(obs[ok]))^2)
  }
}
keep <- !is.na(brier)
cat(sprintf("scored sites: %d  (days per site: %d-%d)\n",
            sum(keep), min(ndays[keep]), max(ndays[keep])))
cat(sprintf("Brier: min=%.3f  median=%.3f  max=%.3f\n",
            min(brier[keep]), median(brier[keep]), max(brier[keep])))

draw <- function() {
  par(mar = c(2.2, 2.2, 3.0, 1.0))
  zr <- c(0, max(brier[keep]))
  cols <- tim.colors(64)
  colidx <- pmin(64L, pmax(1L, ceiling(brier[keep] / zr[2] * 64)))
  plot(S.p[, 1], S.p[, 2], type = "n", asp = 1,
       xlim = range(S.p[, 1]), ylim = range(S.p[, 2]),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "",
       main = paste0("Site-level Brier score, exceedance of 75 ppb, ",
                     "setting 204 (MRTS k=10)", if (nzchar(suffix)) " [dev]"))
  lines(borders.km)
  points(S.win[keep, 1], S.win[keep, 2], pch = 21, cex = 1.7,
         bg = cols[colidx], col = "grey25", lwd = 0.6)
  points(cincy$x, cincy$y, pch = 17, cex = 1.4)
  text(cincy$x, cincy$y, "Cincinnati", pos = 3, offset = 0.5,
       cex = 1.05, font = 2)
  image.plot(legend.only = TRUE, zlim = zr, col = cols,
             legend.width = 1.0, legend.mar = 4.5)
}

dir.create("plots", showWarnings = FALSE)
pdf(sprintf("plots/cincy-204-brier-sites%s.pdf", suffix), width = 8, height = 7.2)
draw()
dev.off()
png(sprintf("plots/cincy-204-brier-sites%s.png", suffix),
    width = 900, height = 810, res = 110)
draw()
dev.off()

# ---------------------------------------------------------------------
# Variant: Brier score as an interpolated field (thin-plate spline over
# the 66 scored sites, evaluated on the prediction grid), with the AQS
# sites overlaid as black circles.  The field between stations is an
# interpolation for visual reading only -- Brier is measured at sites.
# ---------------------------------------------------------------------
tps <- Tps(S.win[keep, ], brier[keep])
brier.field <- predict(tps, S.p)
brier.field <- pmax(as.vector(brier.field), 0)

draw.field <- function() {
  par(mar = c(2.2, 2.2, 3.0, 1.0))
  quilt.plot(x = S.p[, 1], y = S.p[, 2], brier.field, nx = nx, ny = ny,
             zlim = c(0, max(c(brier.field, brier[keep]))),
             xaxt = "n", yaxt = "n",
             main = paste0("Brier score, exceedance of 75 ppb, ",
                           "setting 204 (MRTS k=10)",
                           if (nzchar(suffix)) " [dev]"))
  lines(borders.km)
  points(S.win[keep, 1], S.win[keep, 2], pch = 1, cex = 1.2,
         col = "black", lwd = 1.4)
  points(cincy$x, cincy$y, pch = 17, cex = 1.4)
  text(cincy$x, cincy$y, "Cincinnati", pos = 3, offset = 0.5,
       cex = 1.05, font = 2)
}

pdf(sprintf("plots/cincy-204-brier-field%s.pdf", suffix), width = 8, height = 7.2)
draw.field()
dev.off()
png(sprintf("plots/cincy-204-brier-field%s.png", suffix),
    width = 900, height = 810, res = 110)
draw.field()
dev.off()
cat("field map: plots/cincy-204-brier-field", suffix, ".pdf\n", sep = "")

out <- data.frame(x = S.win[keep, 1], y = S.win[keep, 2],
                  grid.cell = g[keep], ndays = ndays[keep],
                  brier = brier[keep])
write.csv(out, sprintf("output/cincy-204-brier-sites%s.csv", suffix),
          row.names = FALSE)
cat("map: plots/cincy-204-brier-sites", suffix, ".pdf\n", sep = "")
