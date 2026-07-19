#########################################################################
# Full relative-Brier K profile, per setting x method.
#
# For every fitted setting (2 is K=0-only and excluded) and every method,
# compute the relative Brier score R(m, K, d) = BS(m, K, d) / BS(m, 0, d)
# averaged over datasets, at EVERY K on the grid -- the per-K companion to
# the best-over-K cell in tab:brier of the report.
#
# Output:
#   tex/nonsta_simstudy/brier_byK_table.tex  (LaTeX fragment, tab:brierK)
#   output/brier_byK_table.csv               (long format, mean + SE + n)
#########################################################################

if (basename(getwd()) != "simstudy") setwd("d:/Github/spatial-skew-t/code/analysis/simstudy")

# Auto-detect scored settings (2 is K=0-only, skipped); non-t block (>=13)
# appears as its Phase C caches arrive, so this is re-runnable incrementally.
candidate_settings <- c(1L, 3:21)
settings <- candidate_settings[vapply(candidate_settings, function(s)
  file.exists(sprintf("output/results/scores%d_nonsta.RData", s)), TRUE)]
first_nont <- 13L
K.all    <- c(3, 5, 8, 10, 12, 15, 20, 25, 30, 40, 50)
mlab.tex <- c("Gaussian", "Skew-$t$ K1", "Sym-$t$ K1 T", "Skew-$t$ K5", "Sym-$t$ K5 T")
mlab.csv <- c("Gaussian", "SkewT-K1", "SymT-K1-T", "SkewT-K5", "SymT-K5-T")

to3 <- function(a) if (length(dim(a)) == 4) apply(a, c(2, 3, 4), mean, na.rm = TRUE) else a

# per-dataset relative-to-own-K0 ratios: list(M = mean, S = SE, Ks, n)
relK <- function(a) {
  Ks <- as.numeric(dimnames(a)[[3]]); k0 <- which(Ks == 0)
  keep <- which(apply(a, 1, function(v) any(!is.na(v))))
  a <- a[keep, , , drop = FALSE]
  R <- sweep(a, c(1, 2), a[, , k0], "/")
  list(M  = apply(R, c(2, 3), mean, na.rm = TRUE),
       S  = apply(R, c(2, 3), function(v) sd(v, na.rm = TRUE) / sqrt(sum(!is.na(v)))),
       Ks = Ks, n = length(keep))
}

# Load a setting's brier array, appending the K=30/40/50 extension cache when
# present. Dataset dims differ between caches (setting 1: 50 vs 10 rows), so
# rows are aligned by the `datasets` id vector saved in each cache.
load_brier <- function(s) {
  e <- new.env()
  load(sprintf("output/results/scores%d_nonsta.RData", s), envir = e)
  a <- to3(e$brier.score)
  ids_a <- if (exists("datasets", envir = e)) e$datasets else seq_len(dim(a)[1])
  stopifnot(length(ids_a) == dim(a)[1])
  ext_file <- sprintf("output/results/scores%d_nonsta_K3050.RData", s)
  if (file.exists(ext_file)) {
    ee <- new.env(); load(ext_file, envir = ee)
    b <- to3(ee$brier.score)
    ids_b <- if (exists("datasets", envir = ee)) ee$datasets else seq_len(dim(b)[1])
    stopifnot(identical(dimnames(a)[[2]], dimnames(b)[[2]]),
              length(ids_b) == dim(b)[1], all(ids_b %in% ids_a))
    ab <- array(NA_real_, dim = c(dim(a)[1], dim(a)[2], dim(a)[3] + dim(b)[3]),
                dimnames = list(NULL, dimnames(a)[[2]],
                                c(dimnames(a)[[3]], dimnames(b)[[3]])))
    ab[, , seq_len(dim(a)[3])] <- a
    ab[match(ids_b, ids_a), , dim(a)[3] + seq_len(dim(b)[3])] <- b
    a <- ab
  }
  a
}

res <- list(); csv <- NULL
for (s in settings) {
  r <- relK(load_brier(s))
  res[[as.character(s)]] <- r
  for (m in 1:5) for (j in seq_along(r$Ks)) {
    if (r$Ks[j] == 0) next
    csv <- rbind(csv, data.frame(setting = s, method = mlab.csv[m], K = r$Ks[j],
                                 rel.brier = r$M[m, j], se = r$S[m, j], n = r$n))
  }
}
write.csv(csv, "output/brier_byK_table.csv", row.names = FALSE)

# ---- LaTeX fragment ---------------------------------------------------
fmt <- function(r, m, K) {
  j <- match(K, r$Ks)
  if (is.na(j) || is.na(r$M[m, j])) return("--")
  star <- !is.na(r$S[m, j]) && abs(r$M[m, j] - 1) > 2 * r$S[m, j]
  sprintf("$%.3f%s$", r$M[m, j], if (star) "^{*}" else "")
}

ncol.tot <- 2 + length(K.all)

# Build one K-profile table over a chosen subset of method panels AND a chosen
# subset of settings.  The full grid is split three ways so each table fits a
# page: the K1-knot methods are split by setting block (skew-t settings vs the
# non-t extension), and the K5-knot methods form a third table.
build_brierK_table <- function(methods, sett, label, panel_group, block_note) {
  cap <- paste0(
    "\\caption{\\textbf{Relative Brier score at every $\\Kmrts$}, ",
    panel_group, ", ", block_note, " (vs the same method's own $K=0$, mean ",
    "over datasets; $^{*}$ marks $|\\bar R - 1| > 2\\,\\widehat{\\mathrm{SE}}$). ",
    "A dash means that $K$ was not fitted for that setting (settings 13--21 use ",
    "only $K \\in \\{30,40,50\\}$). Isolated large values of Sym-$t$ K5 T are the ",
    "exponential-fallback outliers discussed in the Caveats and are reported ",
    "as-is.}")
  out <- c(
    "\\begin{table}[p]",
    "\\centering\\footnotesize",
    "\\setlength{\\tabcolsep}{3.5pt}",
    cap,
    sprintf("\\label{%s}", label),
    paste0("\\begin{tabular}{cc", paste(rep("c", length(K.all)), collapse = ""), "}"),
    "\\toprule",
    paste0("ID & $n$ & ", paste0("$K{=}", K.all, "$", collapse = " & "), " \\\\")
  )
  for (m in methods) {
    out <- c(out, "\\midrule",
             sprintf("\\multicolumn{%d}{l}{\\emph{%s}} \\\\[1pt]", ncol.tot, mlab.tex[m]))
    for (s in sett) {
      r <- res[[as.character(s)]]
      out <- c(out, paste0(s, " & ", r$n, " & ",
                           paste(vapply(K.all, function(K) fmt(r, m, K), ""), collapse = " & "),
                           " \\\\"))
    }
  }
  c(out, "\\bottomrule", "\\end{tabular}", "\\end{table}")
}

sk_sett   <- settings[settings <  first_nont]   # skew-t settings (1, 3-12)
nont_sett <- settings[settings >= first_nont]   # non-t extension (13-21)

# Table 6 split in two by setting block; table for K5-knot methods unchanged.
writeLines(
  build_brierK_table(1:3, sk_sett, "tab:brierK",
                     "K1-knot methods (Gaussian, Skew-$t$ K1, Sym-$t$ K1 T)",
                     "skew-$t$ settings"),
  "../../../tex/nonsta_simstudy/brier_byK_table.tex")
writeLines(
  build_brierK_table(1:3, nont_sett, "tab:brierKc",
                     "K1-knot methods (Gaussian, Skew-$t$ K1, Sym-$t$ K1 T)",
                     "non-$t$ extension (Gaussian-GP errors)"),
  "../../../tex/nonsta_simstudy/brier_byK_table_c.tex")
writeLines(
  build_brierK_table(4:5, settings, "tab:brierKb",
                     "K5-knot methods (Skew-$t$ K5, Sym-$t$ K5 T)",
                     "all settings"),
  "../../../tex/nonsta_simstudy/brier_byK_table_b.tex")
cat("written: brier_byK_table.tex (K1/skew-t) + brier_byK_table_c.tex",
    "(K1/non-t) + brier_byK_table_b.tex (K5/all) +",
    "output/brier_byK_table.csv (", nrow(csv), "rows )\n")

# sanity anchors (values already reported elsewhere)
chk <- function(s, m, K) round(res[[as.character(s)]]$M[m, match(K, res[[as.character(s)]]$Ks)], 3)
cat(sprintf("check: s11 m2 K12 = %.3f (expect ~1.163)\n", chk(11, 2, 12)))
cat(sprintf("check: s12 m2 K50 = %.3f (expect ~0.925)\n", chk(12, 2, 50)))
cat(sprintf("check: s12 m3 K50 = %.3f (expect ~0.880)\n", chk(12, 3, 50)))
