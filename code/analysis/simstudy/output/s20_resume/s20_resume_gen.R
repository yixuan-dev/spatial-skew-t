# Generate the minimal set of run-settings.R calls covering exactly the
# missing setting-20 cells (grouped by K; methods with identical missing
# dataset-sets merged).  Emits PowerShell lines to s20_resume_calls.txt.
setwd("d:/Github/spatial-skew-t/code/analysis/simstudy")
miss <- read.csv("output/s20_missing.csv")

spec <- function(v) {                       # compact R index spec for CLI
  v <- sort(unique(v))
  if (length(v) == 1) return(as.character(v))
  if (all(diff(v) == 1)) return(sprintf("%d:%d", v[1], v[length(v)]))
  sprintf("c(%s)", paste(v, collapse = ","))
}

lines <- character(0)
for (K in sort(unique(miss$K))) {
  sub <- miss[miss$K == K, ]
  # group methods by their missing dataset-set
  key <- vapply(split(sub$d, sub$m), function(x) paste(sort(x), collapse = ","), "")
  for (dsig in unique(key)) {
    ms <- as.integer(names(key)[key == dsig])
    ds <- as.integer(strsplit(dsig, ",")[[1]])
    kpart <- if (K == 0) "" else sprintf(' "%d"', K)
    lines <- c(lines, sprintf(
      'Rscript run-settings.R --data=simdata_nonsta.RData --setting=20 "%s" 8 2 "%s"%s',
      spec(ds), spec(ms), kpart))
  }
}
writeLines(lines, "output/s20_resume_calls.txt")
cat(sprintf("%d calls covering %d missing cells\n", length(lines), nrow(miss)))
cat(head(lines, 8), sep = "\n")
