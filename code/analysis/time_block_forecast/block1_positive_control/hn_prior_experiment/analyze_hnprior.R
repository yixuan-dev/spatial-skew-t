#########################################################################
# analyze_hnprior.R -- readout of the HN-prior experiment.
#
# Per cell: lambda-hat, beta0-hat, sdz_ratio (A'), max predictive SD /
# marginal SD, assertions A/A'/B/C. Cross-referenced against the guarded
# study's attempt-0 (../results/blk1-*.RData): did the HN prior heal the
# cells the guard had to rescue, and did it leave the healthy ones alone?
#
#   & $R analyze_hnprior.R
#########################################################################

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1]),
    winslash = "/", mustWork = FALSE
  ))
  if (dir.exists(script_dir)) setwd(script_dir)
}

files <- Sys.glob("results/hn-*.RData")
if (!length(files)) stop("no results yet")

rows <- list()
for (f in files) {
  e <- new.env(parent = emptyenv()); load(f, envir = e)
  chk <- e$res$chk
  marg <- chk$sd_lead_max / 3 / ifelse(isTRUE(chk$C_spread), 1, 1) # not used
  rows[[f]] <- data.frame(
    setting = chk$setting, method = chk$method, dataset = chk$dataset,
    lambda = chk$lambda, beta0 = chk$beta0, sdz_ratio = chk$sdz_ratio,
    sd_lead_max = chk$sd_lead_max,
    A = chk$A_zconsist, Ap = chk$Aprime_sdz, B = chk$B_truth, C = chk$C_spread,
    pass = chk$pass
  )
}
tab <- do.call(rbind, rows)
tab <- tab[order(tab$setting, tab$method, tab$dataset), ]
rownames(tab) <- NULL

# guarded study attempt-0 status for the same cells
g0 <- list()
for (f in Sys.glob("../results/blk1-*.RData")) {
  e <- new.env(parent = emptyenv()); load(f, envir = e)
  chk <- e$res$chk
  if (is.null(chk) || !("attempt" %in% names(chk))) next
  a0 <- chk[chk$attempt == 0, , drop = FALSE]
  if (nrow(a0) != 1) next
  g0[[f]] <- data.frame(
    setting = a0$setting, method = a0$method, dataset = a0$dataset,
    lam0_N020 = a0$lambda, pass0_N020 = a0$pass,
    attempts_used = nrow(chk), kept_pass = chk$pass[nrow(chk)]
  )
}
g0 <- do.call(rbind, g0)
tab <- merge(tab, g0, by = c("setting", "method", "dataset"),
             all.x = TRUE, sort = FALSE)
tab <- tab[order(tab$setting, tab$method, tab$dataset), ]

num <- c("lambda", "beta0", "sdz_ratio", "sd_lead_max", "lam0_N020")
tab[num] <- lapply(tab[num], function(v) round(v, 2))

cat("==== HN(0,20) prior, single fit per cell, attempt-0 seeds ====\n\n")
print(tab, row.names = FALSE)

hp <- function(x) sprintf("%d/%d", sum(x, na.rm = TRUE), sum(!is.na(x)))
cat("\n---- headline ----\n")
cat("healthy (B & C)      :", hp(tab$pass), "\n")
cat("lambda sign correct  :", hp(tab$lambda > 0), " (must be 40/40 by construction)\n")
cat("A' in healthy band   :", hp(tab$sdz_ratio > 0.9 & tab$sdz_ratio < 1.3), "\n")

cat("\n---- vs guarded study attempt-0 (N(0,20) prior, same seeds) ----\n")
was_bad <- !tab$pass0_N020
cat("cells attempt-0 FAILED under N(0,20):", sum(was_bad, na.rm = TRUE), "\n")
cat("  ...healed by HN prior alone       :", hp(tab$pass[which(was_bad)]), "\n")
cat("cells attempt-0 passed under N(0,20):", sum(!was_bad, na.rm = TRUE), "\n")
cat("  ...still healthy under HN prior   :", hp(tab$pass[which(!was_bad)]), "\n")

flag <- tab[!tab$pass, , drop = FALSE]
if (nrow(flag)) {
  cat("\n---- cells still failing under HN prior ----\n")
  print(flag[, c("setting", "method", "dataset", "lambda", "beta0",
                 "sdz_ratio", "B", "C", "lam0_N020", "kept_pass")],
        row.names = FALSE)
}
