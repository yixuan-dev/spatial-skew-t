# block_check_400.R
# ---------------------------------------------------------------------------
# Spatial block-bootstrap sensitivity for the 400/400 (two-fold) MRTS stars.
# Same per-site SSE differences as brier_bootstrap_400.R, but resampling
# contiguous geographic blocks of sites (grid on the S coordinates) instead of
# single sites, so spatially correlated neighbours stay together. Blocks cut
# across folds, which is fine: the block, with whatever mix of fold-1/fold-2
# sites it contains, is the approximately independent unit.
# ---------------------------------------------------------------------------

rm(list = ls())
.this <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), winslash = "/", mustWork = FALSE)) else "."
})
setwd(.this)
.r_root <- normalizePath(file.path(.this, "../../../R"), winslash = "/", mustWork = FALSE)
source(local({ ar2 <- file.path(.r_root, "ar2", "auxfunctions.R"); if (file.exists(ar2)) ar2 else file.path(.r_root, "auxfunctions.R") }))

target_probs <- c(0.90, 0.95, 0.98, 0.99, 0.995)
B <- 2000L; seed <- 2024L
pairs <- list(c(111, 61, "AR(2) vs AR(1)"), c(204, 51, "MRTS vs none"))
grids <- list(site = NULL, block_coarse = c(4, 3), block_med = c(6, 4), block_fine = c(8, 6))

se <- new.env(); load("us-all-setup.RData", envir = se)
Y <- get("Y", se); cv.lst <- get("cv.lst", se); S <- get("S", se)
thresholds <- quantile(Y, probs = target_probs, na.rm = TRUE)

site_sse <- function(sid) {
  fe <- new.env(); load(sprintf("results/us-all-%d.RData", sid), envir = fe)
  fit <- get("fit", fe)
  out <- vector("list", 2L)
  for (d in 1:2) {
    val <- cv.lst[[d]]; Yv <- Y[val, , drop = FALSE]
    out[[d]] <- BrierScoreSite(fit[[d]]$yp, thresholds, Yv) * rowSums(!is.na(Yv))
  }
  cat(sprintf("  scored %-4d\n", sid))
  rbind(out[[1]], out[[2]])
}
need <- unique(unlist(lapply(pairs, function(p) as.integer(p[1:2]))))
SSE <- setNames(lapply(need, site_sse), as.character(need))

sites_all <- c(cv.lst[[1]], cv.lst[[2]])
coords <- S[sites_all, , drop = FALSE]
nday <- rowSums(!is.na(Y[sites_all, , drop = FALSE]))

assign_blocks <- function(co, nx, ny) {
  xr <- range(co[, 1]); yr <- range(co[, 2])
  xb <- pmin(nx, 1L + floor(nx * (co[, 1] - xr[1]) / (diff(xr) + 1e-9)))
  yb <- pmin(ny, 1L + floor(ny * (co[, 2] - yr[1]) / (diff(yr) + 1e-9)))
  as.integer((yb - 1L) * nx + xb)
}
boot_ratio <- function(d, w, grp) {
  units <- unique(grp); idx <- split(seq_along(d), grp)
  set.seed(seed)
  bootv <- replicate(B, { drawn <- sample(units, replace = TRUE)
    s <- unlist(idx[as.character(drawn)], use.names = FALSE); sum(d[s]) / sum(w[s]) })
  ci <- quantile(bootv, c(.025, .975), names = FALSE)
  list(Delta = sum(d) / sum(w), lo = ci[1], hi = ci[2],
       p = min(1, 2 * min(mean(bootv <= 0), mean(bootv >= 0))), n_units = length(units))
}

rows <- list()
for (p in pairs) {
  A <- as.character(as.integer(p[1])); Bx <- as.character(as.integer(p[2]))
  for (k in seq_along(target_probs)) {
    d <- SSE[[A]][, k] - SSE[[Bx]][, k]
    for (g in names(grids)) {
      grp <- if (is.null(grids[[g]])) seq_along(d) else assign_blocks(coords, grids[[g]][1], grids[[g]][2])
      r <- boot_ratio(d, nday, grp)
      rows[[length(rows) + 1]] <- data.frame(pair = p[3], level = target_probs[k],
        scheme = g, n_units = r$n_units, Delta = r$Delta, CI_lo = r$lo, CI_hi = r$hi,
        boot_p = r$p, sig = if (r$lo > 0 || r$hi < 0) "*" else "", stringsAsFactors = FALSE)
    }
  }
}
out <- do.call(rbind, rows)
fmt <- out
for (c0 in c("Delta", "CI_lo", "CI_hi")) fmt[[c0]] <- formatC(out[[c0]], format = "fg", digits = 2, flag = "+")
fmt$boot_p <- sprintf("%.3f", out$boot_p)
cat("\n=== 400/400 site vs geographic-block bootstrap ===\n\n")
print(fmt[, c("pair", "level", "scheme", "n_units", "Delta", "CI_lo", "CI_hi", "boot_p", "sig")], row.names = FALSE)
write.csv(out, "output/block_check_400.csv", row.names = FALSE)
cat("\nWritten: output/block_check_400.csv\n")
