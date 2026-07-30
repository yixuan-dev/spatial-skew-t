#########################################################################
# compare_scores_hnprior.R -- score-level HN vs guarded comparison.
#
# Both studies' forecasts are on disk; no refitting:
#   guarded study : ../results/blk1-<s>-<m>-<d>.RData (kept chains)
#   HN study      : results/hn-<s>-<m>-<d>.RData (single fit, HN prior)
#
# Protocol mirrors guard_effect.R / analyze_blk1.R: CRPS and Brier at the
# validation block's q95 threshold, means over leads 1-5 (primary) and
# 1-15. Readouts:
#   A  per-cell paired scores
#   B  group means by setting x method
#   C  the positive-control headline (AR2 vs iid) recomputed inside each
#      study on its own clean set -- does -31% CRPS at setting 5 survive?
#   D  paired guarded-vs-HN on cells clean in BOTH studies -- what does
#      the prior cost/buy in score?
#
#   & $R compare_scores_hnprior.R
#########################################################################

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1]),
    winslash = "/", mustWork = FALSE))
  if (dir.exists(script_dir)) setwd(script_dir)
}
if (!requireNamespace("scoringRules", quietly = TRUE)) {
  stop("needs scoringRules", call. = FALSE)
}

probs <- c(0.90, 0.92, 0.94, 0.95, 0.96, 0.98, 0.99)
qi95 <- which(probs == 0.95)
lead_short <- 1:5

score_case <- function(yhat, y_val) {   # mirror of analyze_blk1.R
  Hh <- dim(yhat)[3]
  crps_h <- vapply(seq_len(Hh), function(h) {
    mean(scoringRules::crps_sample(y = y_val[, h], dat = t(yhat[, , h])),
         na.rm = TRUE)
  }, numeric(1))
  thr <- quantile(y_val, probs = probs, na.rm = TRUE)
  b95_h <- vapply(seq_len(Hh), function(h) {
    phat <- colMeans(yhat[, , h] > thr[qi95])
    mean((as.numeric(y_val[, h] > thr[qi95]) - phat)^2, na.rm = TRUE)
  }, numeric(1))
  c(crps_s = mean(crps_h[lead_short]), crps_a = mean(crps_h),
    b95_s = mean(b95_h[lead_short]), b95_a = mean(b95_h))
}

grid <- expand.grid(setting = c(4, 5), m = 1:2, d = 1:10,
                    KEEP.OUT.ATTRS = FALSE)
grid <- grid[order(grid$setting, grid$m, grid$d), ]

rows <- list()
for (i in seq_len(nrow(grid))) {
  st <- grid$setting[i]
  m <- grid$m[i]
  d <- grid$d[i]

  eh <- new.env(parent = emptyenv())
  load(sprintf("results/hn-%d-%d-%d.RData", st, m, d), envir = eh)
  sh <- score_case(eh$res$yhat, eh$res$y_val)
  chk_h <- eh$res$chk
  rm(eh)

  eg <- new.env(parent = emptyenv())
  load(sprintf("../results/blk1-%d-%d-%d.RData", st, m, d), envir = eg)
  sg <- score_case(eg$res$yhat, eg$res$y_val)
  chk_g <- eg$res$chk
  chk_g <- chk_g[nrow(chk_g), , drop = FALSE]   # kept attempt
  rm(eg)
  gc()

  rows[[i]] <- data.frame(
    setting = st, method = m, dataset = d,
    lam_hn = chk_h$lambda, pass_hn = chk_h$pass,
    lam_g = chk_g$lambda, pass_g = chk_g$pass,
    attempts_g = chk_g$attempt + 1L,
    crps_hn = sh["crps_s"], crps_g = sg["crps_s"],
    b95_hn = sh["b95_s"], b95_g = sg["b95_s"],
    crps_hn_a = sh["crps_a"], crps_g_a = sg["crps_a"],
    b95_hn_a = sh["b95_a"], b95_g_a = sg["b95_a"])
  cat(sprintf(
    "s%d d%2d m%d: lam %+5.2f/%+5.2f | CRPS(1-5) hn %6.3f g %6.3f | B95(1-5) hn %7.4f g %7.4f\n",
    st, d, m, chk_h$lambda, chk_g$lambda,
    sh["crps_s"], sg["crps_s"], sh["b95_s"], sg["b95_s"]))
}
tab <- do.call(rbind, rows)
rownames(tab) <- NULL
save(tab, file = "scores_hn_vs_guarded.RData")

cat("\n=== B. group means, leads 1-5 (all 40 cells) ===\n")
for (st in c(4, 5)) {
  for (m in 1:2) {
    sub <- tab[tab$setting == st & tab$method == m, ]
    cat(sprintf(
      "setting %d method %d: CRPS hn %6.3f g %6.3f | Brier95 hn %7.4f g %7.4f\n",
      st, m, mean(sub$crps_hn), mean(sub$crps_g),
      mean(sub$b95_hn), mean(sub$b95_g)))
  }
}

cat("\n=== C. positive-control headline inside each study (clean sets) ===\n")
for (st in c(4, 5)) {
  for (study in c("guarded", "hn")) {
    w <- tab[tab$setting == st, ]
    m1 <- w[w$method == 1, ]
    m2 <- w[w$method == 2, ]
    m1 <- m1[order(m1$dataset), ]
    m2 <- m2[order(m2$dataset), ]
    if (study == "guarded") {
      keep <- m1$pass_g & m2$pass_g
      c1 <- m1$crps_g[keep]; c2 <- m2$crps_g[keep]
      b1 <- m1$b95_g[keep]; b2 <- m2$b95_g[keep]
    } else {
      keep <- m1$pass_hn & m2$pass_hn
      c1 <- m1$crps_hn[keep]; c2 <- m2$crps_hn[keep]
      b1 <- m1$b95_hn[keep]; b2 <- m2$b95_hn[keep]
    }
    pw <- tryCatch(
      wilcox.test(c2, c1, paired = TRUE, alternative = "less")$p.value,
      error = function(e) NA_real_)
    pb <- tryCatch(
      wilcox.test(b2, b1, paired = TRUE, alternative = "less")$p.value,
      error = function(e) NA_real_)
    cat(sprintf(
      "s%d %-7s (n=%2d: d%s): CRPS iid %6.3f ar2 %6.3f (%+5.1f%%, p=%.3f) | B95 iid %7.4f ar2 %7.4f (%+5.1f%%, p=%.3f)\n",
      st, study, sum(keep),
      paste(m1$dataset[keep], collapse = ","),
      mean(c1), mean(c2), 100 * (mean(c2) / mean(c1) - 1), pw,
      mean(b1), mean(b2), 100 * (mean(b2) / mean(b1) - 1), pb))
  }
}

cat("\n=== D. paired guarded vs HN, cells clean in BOTH (leads 1-5) ===\n")
both <- tab[tab$pass_hn & tab$pass_g, ]
dc <- both$crps_hn - both$crps_g
db <- both$b95_hn - both$b95_g
cat(sprintf("n = %d cells\n", nrow(both)))
cat(sprintf("CRPS  : hn %6.3f  g %6.3f  (mean diff %+6.3f, Wilcoxon two-sided p=%.3f)\n",
            mean(both$crps_hn), mean(both$crps_g), mean(dc),
            wilcox.test(dc)$p.value))
cat(sprintf("Brier : hn %7.4f  g %7.4f  (mean diff %+7.4f, Wilcoxon two-sided p=%.3f)\n",
            mean(both$b95_hn), mean(both$b95_g), mean(db),
            wilcox.test(db)$p.value))

cat("\nwrote scores_hn_vs_guarded.RData\n")
