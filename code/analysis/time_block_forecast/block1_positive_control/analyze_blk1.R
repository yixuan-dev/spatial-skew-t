#########################################################################
# analyze_blk1.R -- paired AR(2) vs i.i.d. comparison for the block-1
# positive control. Reads results/blk1-<setting>-<method>-<dataset>.RData
# written by run_block1_guarded.R.
#
# PRIMARY endpoint (pre-registered): per setting, the paired difference
# (AR2 - iid) of q95 Brier and of CRPS, averaged over leads 1..5, one-sided
# paired t + Wilcoxon with H1: AR2 better (difference < 0). Short leads
# carry the AR(2) signal (kappa_z(h) decays geometrically).
# SECONDARY: leads 1..15 pooled, per-lead curves, dose response
# (setting 4 vs 5), per-dataset detail, sensitivity including flagged fits.
#
#   & $R analyze_blk1.R
#########################################################################

rm(list = ls())
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1]),
    winslash = "/", mustWork = FALSE))
  if (dir.exists(script_dir)) setwd(script_dir)
}
if (!requireNamespace("scoringRules", quietly = TRUE)) {
  stop("analyze_blk1.R requires the 'scoringRules' package.", call. = FALSE)
}

cli <- commandArgs(trailingOnly = TRUE)
getflag <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), cli, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
settings <- eval(parse(text = paste0("c(", getflag("settings", "4,5"), ")")))
setting_labels <- c(`1` = "iid", `2` = "weak", `3` = "moderate",
                    `4` = "strong", `5` = "nearunit",
                    `6` = "lm35", `7` = "lm45")
suffix <- paste0("-set", paste(settings, collapse = "_"))
datasets <- 1:10
methods <- 1:2                     # 1 = iid, 2 = AR(2)
H <- 15L
probs <- c(0.90, 0.92, 0.94, 0.95, 0.96, 0.98, 0.99)
lead_short <- 1:5                  # primary-endpoint lead window
qi95 <- which(probs == 0.95)

# ---- scoring (mirror of simstudy/scores.R::score_block) ----------------
score_block <- function(yhat, y_val) {
  ns <- dim(yhat)[2]
  Hh <- dim(yhat)[3]
  crps_h <- numeric(Hh)
  for (h in seq_len(Hh)) {
    crps_h[h] <- mean(scoringRules::crps_sample(
      y = y_val[, h], dat = t(yhat[, , h])), na.rm = TRUE)
  }
  thr <- quantile(y_val, probs = probs, na.rm = TRUE)
  brier_h <- matrix(NA_real_, Hh, length(probs))
  for (h in seq_len(Hh)) {
    for (q in seq_along(probs)) {
      phat <- colMeans(yhat[, , h] > thr[q])
      ind <- as.numeric(y_val[, h] > thr[q])
      brier_h[h, q] <- mean((ind - phat)^2, na.rm = TRUE)
    }
  }
  list(crps = crps_h, brier = brier_h)
}

# ---- load + score -------------------------------------------------------
crps <- array(NA_real_, c(H, length(datasets), 2, length(settings)),
              dimnames = list(NULL, datasets, c("iid", "ar2"), settings))
brier <- array(NA_real_, c(H, length(probs), length(datasets), 2, length(settings)),
               dimnames = list(NULL, probs, datasets, c("iid", "ar2"), settings))
diag_rows <- list()

for (si in seq_along(settings)) {
  for (di in seq_along(datasets)) {
    for (mi in seq_along(methods)) {
      f <- sprintf("results/blk1-%d-%d-%d.RData",
                   settings[si], methods[mi], datasets[di])
      if (!file.exists(f)) {
        cat("missing:", f, "\n")
        next
      }
      env <- new.env(parent = emptyenv())
      load(f, envir = env)
      sc <- score_block(env$res$yhat, env$res$y_val)
      crps[, di, mi, si] <- sc$crps
      brier[, , di, mi, si] <- sc$brier
      diag_rows[[length(diag_rows) + 1]] <- env$res$chk
      rm(env)
      gc()
    }
  }
  cat(sprintf("setting %d scored\n", settings[si]))
}
diag <- do.call(rbind, diag_rows)

# datasets excluded from the primary analysis: either method flagged
excluded <- unique(diag$dataset[!diag$pass])
excl_by_setting <- lapply(settings, function(s)
  unique(diag$dataset[diag$setting == s & !diag$pass]))
names(excl_by_setting) <- settings

# ---- paired tests -------------------------------------------------------
one_sided <- function(dvec) {
  # dvec = per-dataset paired differences (AR2 - iid); H1: mean < 0
  dvec <- dvec[!is.na(dvec)]
  if (length(dvec) < 3) return(c(n = length(dvec), mean = mean(dvec),
                                 t_p = NA, w_p = NA))
  c(n = length(dvec), mean = mean(dvec),
    t_p = t.test(dvec, alternative = "less")$p.value,
    w_p = wilcox.test(dvec, alternative = "less", exact = FALSE)$p.value)
}

report_setting <- function(si, use_sets_idx, tag) {
  s <- settings[si]
  # per-dataset averages over lead windows
  b95_short <- apply(brier[lead_short, qi95, use_sets_idx, , si, drop = FALSE],
                     3:4, mean)
  b95_all <- apply(brier[, qi95, use_sets_idx, , si, drop = FALSE], 3:4, mean)
  cr_short <- apply(crps[lead_short, use_sets_idx, , si, drop = FALSE], 2:3, mean)
  cr_all <- apply(crps[, use_sets_idx, , si, drop = FALSE], 2:3, mean)

  cat(sprintf("\n--- setting %d [%s: n=%d datasets] ---\n",
              s, tag, length(use_sets_idx)))
  pr <- function(nm, m) {
    d <- m[, 2] - m[, 1]
    st <- one_sided(d)
    cat(sprintf(
      "%-24s iid %8.4f  ar2 %8.4f  gap %+9.5f  t(p)=%.3f  W(p)=%.3f\n",
      nm, mean(m[, 1]), mean(m[, 2]), st["mean"], st["t_p"], st["w_p"]))
    invisible(c(iid = mean(m[, 1]), ar2 = mean(m[, 2]), st))
  }
  r1 <- pr("Brier q95, leads 1-5*", b95_short)
  r2 <- pr("CRPS,      leads 1-5*", cr_short)
  r3 <- pr("Brier q95, leads 1-15", b95_all)
  r4 <- pr("CRPS,      leads 1-15", cr_all)
  list(b95_short = r1, cr_short = r2, b95_all = r3, cr_all = r4)
}

cat("\n==================== PRIMARY (guarded fits only) ====================")
cat("\n(* = pre-registered primary endpoint; one-sided H1: AR2 better)\n")
primary <- list()
for (si in seq_along(settings)) {
  keep <- setdiff(seq_along(datasets),
                  match(excl_by_setting[[as.character(settings[si])]], datasets))
  primary[[si]] <- report_setting(si, keep, "clean")
}

if (length(excluded)) {
  cat("\n================= SENSITIVITY (all fits incl. flagged) ==============\n")
  for (si in seq_along(settings)) report_setting(si, seq_along(datasets), "all")
} else {
  cat("\n(no flagged fits -> sensitivity run identical to primary)\n")
}

# ---- per-dataset detail -------------------------------------------------
cat("\n==================== per-dataset detail (Brier q95, leads 1-5) ======\n")
cat(sprintf("%-8s %-5s | %9s %9s %9s | %s\n", "setting", "dset",
            "iid", "ar2", "gap", "flags (attempts)"))
for (si in seq_along(settings)) {
  for (di in seq_along(datasets)) {
    b <- apply(brier[lead_short, qi95, di, , si, drop = FALSE], 4, mean)
    dsub <- diag[diag$setting == settings[si] & diag$dataset == datasets[di], ]
    fl <- paste(sprintf("m%d:%s(%d)", dsub$method,
                        ifelse(dsub$pass, "ok", "FLAG"), dsub$attempt + 1L),
                collapse = " ")
    cat(sprintf("%-8d %-5d | %9.4f %9.4f %+9.5f | %s\n",
                settings[si], datasets[di], b[1], b[2], b[2] - b[1], fl))
  }
}

# ---- A' distribution (first field data for the ridge-transverse check) --
cat("\nA' (sdz_ratio) by method x pass:\n")
print(aggregate(sdz_ratio ~ setting + method + pass, diag,
                function(v) round(c(mean = mean(v), min = min(v)), 2)))

# ---- plot: per-lead curves ----------------------------------------------
dir.create("output/plots", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
pdf(sprintf("output/plots/blk1_ar2_vs_iid%s.pdf", suffix), width = 10, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1))
for (si in seq_along(settings)) {
  keep <- setdiff(seq_along(datasets),
                  match(excl_by_setting[[as.character(settings[si])]], datasets))
  for (metric in c("brier", "crps")) {
    if (metric == "brier") {
      m_iid <- apply(brier[, qi95, keep, 1, si, drop = FALSE], 1, mean)
      m_ar2 <- apply(brier[, qi95, keep, 2, si, drop = FALSE], 1, mean)
      ylab <- "Brier q95"
    } else {
      m_iid <- apply(crps[, keep, 1, si, drop = FALSE], 1, mean)
      m_ar2 <- apply(crps[, keep, 2, si, drop = FALSE], 1, mean)
      ylab <- "CRPS"
    }
    ylim <- range(m_iid, m_ar2)
    plot(1:H, m_iid, type = "b", pch = 16, col = "#3465a4", lwd = 2,
         xlab = "lead h", ylab = ylab, ylim = ylim,
         main = sprintf("setting %d (%s), clean n=%d", settings[si],
                        setting_labels[as.character(settings[si])],
                        length(keep)))
    lines(1:H, m_ar2, type = "b", pch = 17, col = "#d1352b", lwd = 2)
    abline(v = max(lead_short) + 0.5, lty = 3, col = "gray60")
    legend("bottomright", bty = "n", pch = c(16, 17), lwd = 2,
           col = c("#3465a4", "#d1352b"), legend = c("iid", "AR(2)"))
  }
}
dev.off()

# ---- tables -------------------------------------------------------------
gap_tab <- do.call(rbind, lapply(seq_along(settings), function(si) {
  keep <- setdiff(seq_along(datasets),
                  match(excl_by_setting[[as.character(settings[si])]], datasets))
  do.call(rbind, lapply(seq_along(probs), function(qi) {
    bs <- apply(brier[lead_short, qi, keep, , si, drop = FALSE], 3:4, mean)
    ba <- apply(brier[, qi, keep, , si, drop = FALSE], 3:4, mean)
    data.frame(setting = settings[si], prob = probs[qi], n = length(keep),
               gap_lead1to5 = mean(bs[, 2] - bs[, 1]),
               t_p_1to5 = one_sided(bs[, 2] - bs[, 1])["t_p"],
               gap_lead1to15 = mean(ba[, 2] - ba[, 1]),
               t_p_1to15 = one_sided(ba[, 2] - ba[, 1])["t_p"])
  }))
}))
write.csv(gap_tab, sprintf("output/tables/blk1_paired%s.csv", suffix),
          row.names = FALSE)
write.csv(diag, sprintf("output/tables/blk1_diag%s.csv", suffix),
          row.names = FALSE)
cat(sprintf("\nwrote output/plots/blk1_ar2_vs_iid%s.pdf and tables%s\n",
            suffix, suffix))
