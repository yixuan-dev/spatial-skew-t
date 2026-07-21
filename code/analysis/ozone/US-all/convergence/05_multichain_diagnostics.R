# convergence/05_multichain_diagnostics.R
# ---------------------------------------------------------------------------
# True multi-chain diagnostics for the settings rerun by 04_multichain_run.R.
# Chain 1 is the production fold-1 chain (pulled from the 01 cache); chains
# 2-4 are the overdispersed reruns. Computes 4-chain rank-normalized split
# R-hat (each chain halved -> 8 pseudo-chains), coda multi-chain ESS, and
# per-chain means/SDs, plus overlaid trace and density figures.
#
# If chains differ in length (e.g. dev-mode reruns against the 5000-draw
# production chain), all chains are truncated to the shortest with a loud
# warning; results are then smoke-test only.
#
# Usage: Rscript convergence/05_multichain_diagnostics.R [setting ...]
#        (default: every setting with at least one multichain-*.rds)
# Output: output/us-all/tables/convergence_multichain.csv
#         output/us-all/plots/convergence/trace-multichain-<N>-fold1.png
#         output/us-all/plots/convergence/density-multichain-<N>-fold1.png
# ---------------------------------------------------------------------------

rm(list = ls())
.this <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), winslash = "/", mustWork = FALSE)) else "."
})
setwd(file.path(.this, ".."))
source("convergence/00_conv_lib.R")

paths <- conv_paths()

cli <- suppressWarnings(as.integer(commandArgs(trailingOnly = TRUE)))
cli <- cli[!is.na(cli)]
targets <- if (length(cli)) cli else {
  found <- list.files(paths$cache, pattern = "^multichain-\\d+-fold1-chain\\d\\.rds$")
  sort(unique(as.integer(sub("^multichain-(\\d+)-.*$", "\\1", found))))
}
if (!length(targets)) stop("no multichain-*.rds found; run 04_multichain_run.R first", call. = FALSE)

chain_cols <- c("black", "#0072B2", "#D55E00", "#009E73")

all_rows <- list()
for (s in targets) {
  slim <- load_setting_slim(s, paths)
  chains_by_id <- list("1" = slim$folds[[1]])
  meta <- list()
  for (cid in 2:4) {
    mf <- multichain_file(s, cid, paths)
    if (file.exists(mf)) {
      mc <- readRDS(mf)
      chains_by_id[[as.character(cid)]] <- mc$chains
      meta[[as.character(cid)]] <- mc
    }
  }
  M <- length(chains_by_id)
  if (M < 2) {
    cat(sprintf("setting %d: only the production chain available, skipping\n", s))
    next
  }
  dev_reruns <- any(vapply(meta, function(m) !identical(m$run_mode, "prod"), TRUE))
  if (dev_reruns) warning(sprintf(
    "setting %d: some reruns are dev-mode; treat results as smoke test only", s))

  # align parameters and lengths across chains
  common <- Reduce(intersect, lapply(chains_by_id, colnames))
  common <- colnames(chains_by_id[[1]])[colnames(chains_by_id[[1]]) %in% common]
  nmin <- min(vapply(chains_by_id, nrow, 0L))
  nmax <- max(vapply(chains_by_id, nrow, 0L))
  if (nmin != nmax) warning(sprintf(
    "setting %d: unequal chain lengths (%d..%d), truncating to the last %d draws",
    s, nmin, nmax, nmin))
  chains_by_id <- lapply(chains_by_id, function(m) {
    m[(nrow(m) - nmin + 1L):nrow(m), common, drop = FALSE]
  })
  ids <- names(chains_by_id)
  cat(sprintf("setting %d: %d chains x %d draws x %d params\n",
              s, M, nmin, length(common)))

  # ---- diagnostics table ----
  # rhat: all chains incl. the production chain (chain 1). rhat_reruns: the
  # fresh overdispersed chains only -- separates "do the reruns agree with
  # each other" from "do they agree with the production chain", which
  # matters when the production chain is stuck (e.g. the frozen-z runs).
  rows <- lapply(common, function(pn) {
    mat <- do.call(cbind, lapply(chains_by_id, function(m) m[, pn]))
    degenerate <- any(apply(mat, 2, sd) == 0)
    r <- if (degenerate) c(rhat_bulk = NA_real_, rhat_tail = NA_real_, rhat = NA_real_)
         else rhat_rank_norm(split_chains_half(mat))
    rerun_idx <- which(ids != "1")
    r_rr <- if (length(rerun_idx) >= 2 &&
                all(apply(mat[, rerun_idx, drop = FALSE], 2, sd) > 0)) {
      rhat_rank_norm(split_chains_half(mat[, rerun_idx, drop = FALSE]))[["rhat"]]
    } else NA_real_
    ess_total <- if (degenerate) NA_real_
                 else ess_coda_multi(lapply(seq_len(M), function(j) mat[, j]))
    stats <- unlist(lapply(seq_len(M), function(j) {
      v <- c(mean(mat[, j]), sd(mat[, j]))
      names(v) <- paste0(c("mean_chain", "sd_chain"), ids[j])
      v
    }))
    df <- data.frame(setting = s, param = pn, nchains = M, ndraws = nmin,
                     rhat_bulk = r[["rhat_bulk"]], rhat_tail = r[["rhat_tail"]],
                     rhat = r[["rhat"]], rhat_reruns = r_rr, ess_total = ess_total,
                     degenerate = degenerate, stringsAsFactors = FALSE)
    cbind(df, as.data.frame(as.list(stats)))
  })
  rows <- do.call(rbind, rows)
  rows$flag <- classify_flag(rows$rhat, rows$ess_total, rows$degenerate)
  all_rows[[length(all_rows) + 1L]] <- rows

  # ---- overlay trace + density figures ----
  P <- length(common)
  nc <- ceiling(sqrt(P)); nr <- ceiling(P / nc)
  iter_idx <- seq_len(nmin) + if (nmin == 5000L) 25000L else 0L
  line_cols <- adjustcolor(chain_cols[seq_len(M)], alpha.f = 0.65)

  out_tr <- file.path(paths$plots, sprintf("trace-multichain-%d-fold1.png", s))
  png(out_tr, width = 480 * nc, height = 340 * nr, res = 110)
  par(mfrow = c(nr, nc), mar = c(3.2, 3.2, 2.6, 0.8), mgp = c(2.1, 0.7, 0),
      oma = c(0, 0, 2.2, 0))
  for (pn in common) {
    mat <- do.call(cbind, lapply(chains_by_id, function(m) m[, pn]))
    row <- rows[rows$param == pn, ]
    matplot(iter_idx, mat, type = "l", lty = 1, lwd = 0.4, col = line_cols,
            xlab = "iteration", ylab = pn, main = "")
    title(main = sprintf("%s  Rhat=%.3f  ESS=%.0f", pn, row$rhat, row$ess_total),
          col.main = if (row$flag %in% c("fail", "degenerate")) "red"
                     else if (row$flag == "marginal") "darkorange" else "black",
          cex.main = 1.0)
    if (pn == common[1]) {
      legend("topright", bty = "n", lty = 1, lwd = 1.4, col = chain_cols[seq_len(M)],
             legend = c("production", paste0("chain ", ids[-1])), cex = 0.8)
    }
  }
  mtext(sprintf("us-all setting %d fold 1: %d chains (rank-normalized split R-hat)", s, M),
        side = 3, line = 0.6, outer = TRUE, cex = 1.0)
  dev.off()
  cat(sprintf("wrote %s\n", out_tr))

  out_de <- file.path(paths$plots, sprintf("density-multichain-%d-fold1.png", s))
  png(out_de, width = 480 * nc, height = 340 * nr, res = 110)
  par(mfrow = c(nr, nc), mar = c(3.2, 3.2, 2.6, 0.8), mgp = c(2.1, 0.7, 0),
      oma = c(0, 0, 2.2, 0))
  for (pn in common) {
    dens <- lapply(chains_by_id, function(m) density(m[, pn]))
    xr <- range(unlist(lapply(dens, function(d) d$x)))
    yr <- c(0, max(unlist(lapply(dens, function(d) d$y))))
    plot(NA, xlim = xr, ylim = yr, xlab = pn, ylab = "density", main = pn)
    for (j in seq_len(M)) lines(dens[[j]], col = chain_cols[j], lwd = 1.4)
    if (pn == common[1]) {
      legend("topright", bty = "n", lty = 1, lwd = 1.4, col = chain_cols[seq_len(M)],
             legend = c("production", paste0("chain ", ids[-1])), cex = 0.8)
    }
  }
  mtext(sprintf("us-all setting %d fold 1: per-chain posterior densities", s),
        side = 3, line = 0.6, outer = TRUE, cex = 1.0)
  dev.off()
  cat(sprintf("wrote %s\n", out_de))
}

if (length(all_rows)) {
  tab <- do.call(rbind, all_rows)
  num_cols <- setdiff(names(tab)[vapply(tab, is.numeric, TRUE)], c("setting", "nchains", "ndraws"))
  tab[num_cols] <- lapply(tab[num_cols], function(x) round(x, 4))
  out_csv <- file.path(paths$tables, "convergence_multichain.csv")
  write.csv(tab, out_csv, row.names = FALSE)
  cat(sprintf("\nwrote %s\n", out_csv))
  print(tab[, c("setting", "param", "rhat", "ess_total", "flag")], row.names = FALSE)
}
