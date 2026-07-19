#########################################################################
# Best-over-K score tables for the report (tab:energy / tab:brier /
# tab:vario in tex/nonsta_simstudy/nonsta_simstudy.tex).
#
# For every setting and method: min over K > 0 of the relative score
# R(m,K) = mean_d S(m,K,d)/S(m,0,d), with the minimising K, an asterisk
# when |R - 1| > 2 SE, and bold when R < 0.75.  The K grid is whatever is
# on disk: base cache scores{s}_nonsta.RData, extended by
# scores{s}_nonsta_K3050.RData when that file exists -- so re-running this
# script after new K fits are scored updates the tables automatically.
#
# Output: tex/nonsta_simstudy/bestK_{energy,brier,vario}.tex  (tabular only)
#########################################################################

if (basename(getwd()) != "simstudy") setwd("d:/Github/spatial-skew-t/code/analysis/simstudy")

# Auto-detect which settings have been scored (setting 2 is K=0-only, skipped).
# Settings 13-21 are the non-t extension and appear only once their Phase C
# score cache exists, so this script can be re-run incrementally.
candidate_settings <- c(1L, 3:21)
settings <- candidate_settings[vapply(candidate_settings, function(s)
  file.exists(sprintf("output/results/scores%d_nonsta.RData", s)), TRUE)]
first_nont <- 13L                      # non-t block starts here (Gaussian-GP)
slab <- c(
  "1"  = "mean (2 bumps, sd 3.7)",
  "3"  = "dependence (low-rank)",
  "4"  = "dependence (PS, repl.\\ $C$)",
  "5"  = "dependence (PS, additive)",
  "6"  = "dependence ($\\rho(\\bs)$, repl.\\ $C$)",
  "7"  = "dependence ($\\rho(\\bs)$, additive)",
  "8"  = "dependence (Fuentes, repl.\\ $C$)",
  "9"  = "tails ($g$-and-$h$)",
  "10" = "skewness ($\\lambda(\\bs)$)",
  "11" = "mean (2 bumps, sd 11)",
  "12" = "mean (6 bumps, sd 11)",
  "13" = "mean, Franke (multiscale)",
  "14" = "mean, curved ridge",
  "15" = "mean, annulus (Ricker)",
  "16" = "mean, sigmoid front",
  "17" = "mean, GP draw",
  "18" = "mean, Franke + sinh-arcsinh err.",
  "19" = "mean, poly2 (bowl+saddle)",
  "20" = "mean, transform mix",
  "21" = "mean, non-smooth (hinges)"
)
score_files <- c(energy = "energy.score", brier = "brier.score", vario = "vario.score")

to3 <- function(a) if (length(dim(a)) == 4) apply(a, c(2, 3, 4), mean, na.rm = TRUE) else a

# The dataset dimension differs between caches (e.g. setting 1: 50 rows in
# the base cache, 10 in the K3050 one), so rows MUST be aligned by dataset id
# (the `datasets` vector saved in each cache) -- plain abind would recycle.
load_score <- function(s, nm) {
  e <- new.env()
  load(sprintf("output/results/scores%d_nonsta.RData", s), envir = e)
  a <- to3(get(nm, envir = e))
  ids_a <- if (exists("datasets", envir = e)) e$datasets else seq_len(dim(a)[1])
  stopifnot(length(ids_a) == dim(a)[1])
  ext_file <- sprintf("output/results/scores%d_nonsta_K3050.RData", s)
  if (file.exists(ext_file)) {
    ee <- new.env(); load(ext_file, envir = ee)
    b <- to3(get(nm, envir = ee))
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

relK <- function(a) {
  Ks <- as.numeric(dimnames(a)[[3]]); k0 <- which(Ks == 0)
  keep <- which(apply(a, 1, function(v) any(!is.na(v))))
  a <- a[keep, , , drop = FALSE]
  R <- sweep(a, c(1, 2), a[, , k0], "/")
  list(M  = apply(R, c(2, 3), mean, na.rm = TRUE),
       S  = apply(R, c(2, 3), function(v) sd(v, na.rm = TRUE) / sqrt(sum(!is.na(v)))),
       Ks = Ks, n = length(keep))
}

best_cell <- function(r, m) {
  jpos <- which(r$Ks > 0 & !is.na(r$M[m, ]))
  if (length(jpos) == 0) return("--")
  j <- jpos[which.min(r$M[m, jpos])]
  val <- r$M[m, j]
  star <- !is.na(r$S[m, j]) && abs(val - 1) > 2 * r$S[m, j]
  core <- sprintf("%.3f%s", val, if (star) "^{*}" else "")
  if (val < 0.75) core <- sprintf("\\mathbf{%s}", core)
  sprintf("$%s$ {\\tiny(%d)}", core, as.integer(r$Ks[j]))
}

header <- paste0(
  "\\setlength{\\tabcolsep}{3.7pt}\n",
  "\\begin{tabular}{clcccccc}\n\\toprule\n",
  "ID & non-stationarity & $n$ & Gaussian & Skew-$t$ K1 & Sym-$t$ K1 T & ",
  "Skew-$t$ K5 & Sym-$t$ K5 T \\\\\n\\midrule"
)

for (sc in names(score_files)) {
  L <- header
  nont_marked <- FALSE
  for (s in settings) {
    if (!nont_marked && s >= first_nont) {
      L <- c(L, "\\midrule",
             paste0("\\multicolumn{8}{l}{\\emph{Non-$t$ extension ",
                    "(Gaussian-GP errors; skew-$t$ is misspecified)}} \\\\"))
      nont_marked <- TRUE
    }
    r <- relK(load_score(s, score_files[[sc]]))
    cells <- vapply(1:5, function(m) best_cell(r, m), "")
    L <- c(L, paste0(s, " & ", slab[[as.character(s)]], " & ", r$n, " & ",
                     paste(cells, collapse = " & "), " \\\\"))
  }
  L <- c(L, "\\bottomrule", "\\end{tabular}")
  out <- sprintf("../../../tex/nonsta_simstudy/bestK_%s.tex", sc)
  writeLines(L, out)
  cat("written:", out, "\n")
}
