#########################################################################
# analyze_allparams.R -- recovery of the parameters that the guard's
# diagnostics never recorded: rho, nu, gamma, beta1, beta2 (plus the
# posterior spread of the ones that were), under lambda ~ HN(0, 20).
#
# Truth (setup.R): beta = (10, 0, 0), lambda = 3, rho = 1, nu = 0.5,
# gamma = 0.9, internal tau shape/rate = (3, 8).
# Sampler bounds: rho.upper = 15, nu.upper = 10.
#
#   & $R analyze_allparams.R
#########################################################################

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1]),
    winslash = "/", mustWork = FALSE))
  if (dir.exists(script_dir)) setwd(script_dir)
}

TRUTH <- c(beta0 = 10, beta1 = 0, beta2 = 0, lambda = 3,
           rho = 1, nu = 0.5, gamma = 0.9,
           tau.alpha = 3, tau.beta = 8)
PHI_TRUTH <- list("4" = c(0.80, -0.35), "5" = c(0.15, 0.80))

files <- Sys.glob("results_allparams/ap-*.RData")
if (!length(files)) stop("no all-params results yet")

rows <- list()
for (f in files) {
  e <- new.env(parent = emptyenv())
  load(f, envir = e)
  p <- e$res$post
  base <- data.frame(setting = e$res$setting, method = e$res$method_id,
                     dataset = e$res$dataset)
  for (nm in names(p)) {
    rows[[length(rows) + 1L]] <- cbind(
      base, param = nm,
      mean = p[[nm]][["mean"]], sd = p[[nm]][["sd"]],
      q025 = p[[nm]][["q025"]], q50 = p[[nm]][["q50"]],
      q975 = p[[nm]][["q975"]])
  }
}
tab <- do.call(rbind, rows)
rownames(tab) <- NULL

truth_of <- function(param, setting) {
  if (param %in% names(TRUTH)) return(TRUTH[[param]])
  ph <- PHI_TRUTH[[as.character(setting)]]
  if (grepl("1$", param)) return(ph[1])
  if (grepl("2$", param)) return(ph[2])
  NA_real_
}
tab$truth <- mapply(truth_of, tab$param, tab$setting)
tab$covered <- tab$q025 <= tab$truth & tab$truth <= tab$q975

ord <- c("beta0", "beta1", "beta2", "lambda", "rho", "nu", "gamma",
         "tau.alpha", "tau.beta",
         "phi.z1", "phi.z2", "phi.tau1", "phi.tau2", "phi.w1", "phi.w2")

for (st in sort(unique(tab$setting))) {
  ph <- PHI_TRUTH[[as.character(st)]]
  w <- tab[tab$setting == st, ]
  cells <- nrow(unique(w[, c("method", "dataset")]))
  cat(sprintf("\n==== setting %d (phi truth %.2f, %+.2f), n = %d cells, AR(2) method ====\n",
              st, ph[1], ph[2], cells))
  cat(sprintf("%-10s %7s %8s %8s %8s   %s\n",
              "param", "truth", "post.mean", "bias", "post.sd",
              "95% CI coverage (cells)"))
  for (nm in ord) {
    v <- w[w$param == nm, ]
    if (!nrow(v)) next
    cat(sprintf("%-10s %7.2f %8.3f %+8.3f %8.3f   %d/%d\n",
                nm, v$truth[1], mean(v$mean), mean(v$mean) - v$truth[1],
                mean(v$sd), sum(v$covered), nrow(v)))
  }
}

cat("\n==== per-cell detail for the newly measured parameters ====\n")
new <- tab[tab$param %in% c("rho", "nu", "gamma", "beta1", "beta2"), ]
new <- new[order(new$param, new$setting, new$dataset), ]
cat(sprintf("%-7s %3s %3s %9s %9s %9s %9s %5s\n",
            "param", "s", "d", "truth", "mean", "q025", "q975", "cov"))
for (i in seq_len(nrow(new))) {
  r <- new[i, ]
  cat(sprintf("%-7s %3d %3d %9.2f %9.3f %9.3f %9.3f %5s\n",
              r$param, r$setting, r$dataset, r$truth, r$mean,
              r$q025, r$q975, r$covered))
}

save(tab, file = "recovery_allparams_hn.RData")
cat("\nwrote recovery_allparams_hn.RData\n")
