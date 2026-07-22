# Exceedance map for setting 204 around Cincinnati: P(at least 1 of the 31
# days exceeds 75 ppb) at each grid point, the make-map-71.R statistic
# restricted to the one map we need (no p.1 / p.2 O(nt^2) loops).
#
# Same caveat as the other make-map-*.R: the product over days treats daily
# exceedances as independent given the marginal posterior probabilities
# (marginals are averaged over draws first), which understates the tail.
#
# Usage:  Rscript make-map-cincy-204.R [dev]
rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)
suffix <- if (length(args) > 0 && tolower(args[1]) == "dev") "-dev" else ""

library(fields)
load("../ozone_data.RData")   # borders (and nothing else is used)
borders.km <- borders / 1000
load(paste0("us-all-pred-cincy-204", suffix, ".RData"))

threshold <- 75
np <- dim(y.pred)[2]
nt <- dim(y.pred)[3]
nx <- length(unique(S.p[, 1]))
ny <- length(unique(S.p[, 2]))
stopifnot(nx * ny == np)

# P(at least one day exceeds 75) per grid point
p.below <- matrix(0, np, nt)
for (i in 1:np) { for (t in 1:nt) {
  p.below[i, t] <- mean(y.pred[, i, t] <= threshold)
} }
p.0        <- apply(p.below, 1, prod)
p.atleast1 <- 1 - p.0
cat(sprintf("P(>=1 day > %d): min=%.3f  median=%.3f  max=%.3f\n",
            threshold, min(p.atleast1), median(p.atleast1), max(p.atleast1)))

dir.create("plots", showWarnings = FALSE)
mapfile <- paste0("plots/cincy-204-exceed75", suffix, ".pdf")
pdf(mapfile, width = 8, height = 7.2)
par(mar = c(2.2, 2.2, 3.0, 1.0))
quilt.plot(
  x = S.p[, 1], y = S.p[, 2], p.atleast1, nx = nx, ny = ny,
  zlim = c(0, 1), xaxt = "n", yaxt = "n",
  main = sprintf("P(at least 1 day > %d ppb), setting 204 (MRTS k=10)", threshold)
)
lines(borders.km)
points(cincy$x, cincy$y, pch = 17, cex = 1.3)
text(cincy$x, cincy$y, "Cincinnati", pos = 3, offset = 0.5, cex = 1.05, font = 2)
dev.off()
cat("map:", mapfile, "\n")

rm(y.pred)
save.image(file = paste0("predict-maps-cincy-204", suffix, ".RData"))
