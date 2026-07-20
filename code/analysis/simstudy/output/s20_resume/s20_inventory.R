# Inventory setting-20 fits vs the 600-cell target grid; validate the newest
# files (possible mid-write corruption at the pause); emit missing combos.
setwd("d:/Github/spatial-skew-t/code/analysis/simstudy")
datasets <- 1:10; methods <- 1:5
Ks <- c(0, 3, 5, 8, 10, 12, 15, 20, 25, 30, 40, 50)

fname <- function(m, d, K) {
  if (K == 0) sprintf("results_nonsta/20-%d-%d.RData", m, d)
  else        sprintf("results_nonsta/20-%d-%d-K%d.RData", m, d, K)
}

grid <- expand.grid(m = methods, d = datasets, K = Ks)
grid$file <- mapply(fname, grid$m, grid$d, grid$K)
grid$exists <- file.exists(grid$file)

# validate the 12 newest existing files (workers interrupted mid-save)
ex <- grid[grid$exists, ]
mt <- file.mtime(ex$file)
newest <- ex$file[order(mt, decreasing = TRUE)][1:min(12, nrow(ex))]
bad <- character(0)
for (f in newest) {
  ok <- tryCatch({ e <- new.env(); load(f, envir = e); !is.null(e$fit.1) },
                 error = function(err) FALSE)
  if (!ok) bad <- c(bad, f)
}
if (length(bad)) {
  cat("CORRUPT (treat as missing):\n"); print(bad)
  grid$exists[grid$file %in% bad] <- FALSE
} else cat("newest 12 files all load cleanly\n")

miss <- grid[!grid$exists, c("m", "d", "K")]
cat(sprintf("existing valid: %d / 600; missing: %d\n", sum(grid$exists), nrow(miss)))

# group missing by (K-set signature per dataset) is complex; instead group by
# dataset: for each dataset list missing (m, K) pairs
for (d in datasets) {
  sub <- miss[miss$d == d, ]
  if (nrow(sub) == 0) next
  cat(sprintf("dataset %d: %d missing -> methods %s | Ks %s\n", d, nrow(sub),
              paste(sort(unique(sub$m)), collapse = ","),
              paste(sort(unique(sub$K)), collapse = ",")))
}
write.csv(miss, "output/s20_missing.csv", row.names = FALSE)
cat("written: output/s20_missing.csv\n")
