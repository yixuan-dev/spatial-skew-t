#########################################################################
# recovery_hnprior.R -- parameter recovery under lambda ~ HN(0, 20).
#
# Truth (setup.R): beta = (10, 0, 0), lambda = 3, tau.alpha = 6,
# tau.beta = 16 -- the generator uses rgamma(shape = tau.alpha/2,
# rate = tau.beta/2), and check_fit_consistency() records a = tau.alpha/2,
# b = tau.beta/2, so the truth for the recorded (a, b) is (3, 8).
# phi (shared by tau*, z*, w*): setting 4 = (0.80, -0.35),
# setting 5 = (0.15, 0.80).
#
# Recovered from the on-disk diagnostics of both studies (no refitting):
# beta0, lambda, a, b, phi.z, phi.tau, plus the two z-moment ratios
# (z_ratio and sdz_ratio, both = 1 under the model).
# rho/nu/gamma/beta1/beta2 are not in chk -- see refit_allparams_hn.R.
#
#   & $R recovery_hnprior.R
#########################################################################

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1]),
    winslash = "/", mustWork = FALSE))
  if (dir.exists(script_dir)) setwd(script_dir)
}

TRUTH <- list(beta0 = 10, lambda = 3, a = 3, b = 8,
              phi = list("4" = c(0.80, -0.35), "5" = c(0.15, 0.80)))

pull <- function(pattern, tag) {
  out <- list()
  for (f in Sys.glob(pattern)) {
    e <- new.env(parent = emptyenv())
    load(f, envir = e)
    chk <- e$res$chk
    chk <- chk[nrow(chk), , drop = FALSE]      # kept attempt
    out[[f]] <- data.frame(
      study = tag, setting = chk$setting, method = chk$method,
      dataset = chk$dataset, pass = chk$pass,
      beta0 = chk$beta0, lambda = chk$lambda, a = chk$a, b = chk$b,
      phi1.z = chk$phi1.z, phi2.z = chk$phi2.z,
      phi1.tau = chk$phi1.tau, phi2.tau = chk$phi2.tau,
      z_ratio = chk$z_ratio, sdz_ratio = chk$sdz_ratio,
      stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

tab <- rbind(pull("results/hn-*.RData", "HN"),
             pull("../results/blk1-*.RData", "guard"))
rownames(tab) <- NULL
# the parent results/ dir also holds cells from other studies (settings 6-7,
# method 4); restrict to the 40-cell block-1 positive-control design.
tab <- tab[tab$setting %in% c(4, 5) & tab$method %in% c(1, 2) &
             tab$dataset %in% 1:10, ]
# implied sigma^2 = b / (a - 1); truth 8 / 2 = 4
tab$Esig2 <- with(tab, ifelse(a > 1, b / (a - 1), NA))

fmt <- function(v, truth, digits = 2) {
  sprintf("%.*f (%+.*f)", digits, mean(v), digits, mean(v) - truth)
}
line <- function(v, truth, digits = 2) {
  sprintf("%7.*f %+7.*f %7.*f  %s",
          digits, mean(v), digits, mean(v) - truth, digits, sd(v),
          sprintf("[%.*f, %.*f]", digits, min(v), digits, max(v)))
}

cat("=== parameter recovery, HN(0,20) prior, clean cells only ===\n")
cat("(truth in brackets; bias = mean - truth; range over cells)\n\n")

for (st in c(4, 5)) {
  ph <- TRUTH$phi[[as.character(st)]]
  cat(sprintf("---- setting %d (phi truth %.2f, %+.2f) ----\n", st, ph[1], ph[2]))
  for (tag in c("HN", "guard")) {
    w <- tab[tab$study == tag & tab$setting == st & tab$pass, ]
    cat(sprintf("  %-5s n = %d clean cells\n", tag, nrow(w)))
    cat(sprintf("    %-16s %7s %7s %7s  %s\n",
                "param [truth]", "mean", "bias", "sd", "range"))
    cat(sprintf("    %-16s %s\n", "beta0     [10]", line(w$beta0, 10)))
    cat(sprintf("    %-16s %s\n", "lambda     [3]", line(w$lambda, 3)))
    cat(sprintf("    %-16s %s\n", "a          [3]", line(w$a, 3)))
    cat(sprintf("    %-16s %s\n", "b          [8]", line(w$b, 8)))
    cat(sprintf("    %-16s %s\n", "E[sig2]    [4]", line(w$Esig2, 4)))
    cat(sprintf("    %-16s %s\n", "z_ratio    [1]", line(w$z_ratio, 1)))
    cat(sprintf("    %-16s %s\n", "sdz_ratio  [1]", line(w$sdz_ratio, 1)))
    m2 <- w[w$method == 2 & !is.na(w$phi1.z), ]
    if (nrow(m2)) {
      cat(sprintf("    -- AR(2) cells only (n = %d) --\n", nrow(m2)))
      cat(sprintf("    %-16s %s\n",
                  sprintf("phi1.z  [%+.2f]", ph[1]), line(m2$phi1.z, ph[1])))
      cat(sprintf("    %-16s %s\n",
                  sprintf("phi2.z  [%+.2f]", ph[2]), line(m2$phi2.z, ph[2])))
      cat(sprintf("    %-16s %s\n",
                  sprintf("phi1.tau[%+.2f]", ph[1]), line(m2$phi1.tau, ph[1])))
      cat(sprintf("    %-16s %s\n",
                  sprintf("phi2.tau[%+.2f]", ph[2]), line(m2$phi2.tau, ph[2])))
    }
    cat("\n")
  }
}

cat("=== paired HN vs guard on cells clean in BOTH (did the prior move anything?) ===\n")
h <- tab[tab$study == "HN", ]
g <- tab[tab$study == "guard", ]
key <- c("setting", "method", "dataset")
m <- merge(h, g, by = key, suffixes = c(".hn", ".g"))
m <- m[m$pass.hn & m$pass.g, ]
cat(sprintf("n = %d cells\n", nrow(m)))
cat(sprintf("%-12s %9s %9s %9s %8s\n", "param", "HN", "guard", "meandiff", "p"))
for (p in c("beta0", "lambda", "a", "b", "Esig2", "z_ratio", "sdz_ratio",
            "phi1.z", "phi2.z", "phi1.tau", "phi2.tau")) {
  vh <- m[[paste0(p, ".hn")]]
  vg <- m[[paste0(p, ".g")]]
  ok <- is.finite(vh) & is.finite(vg)
  if (sum(ok) < 3) next
  d <- vh[ok] - vg[ok]
  pv <- tryCatch(wilcox.test(d)$p.value, error = function(e) NA_real_)
  cat(sprintf("%-12s %9.3f %9.3f %+9.3f %8.3f  (n=%d)\n",
              p, mean(vh[ok]), mean(vg[ok]), mean(d), pv, sum(ok)))
}

cat("\n=== the a < 2 cells (infinite 4th moment of sigma^2; truth a = 3) ===\n")
bad <- tab[tab$a < 2, c("study", "setting", "method", "dataset", "a", "b",
                        "Esig2", "pass")]
print(bad[order(bad$setting, bad$dataset, bad$study), ], row.names = FALSE)

save(tab, file = "recovery_hn_vs_guarded.RData")
cat("\nwrote recovery_hn_vs_guarded.RData\n")
