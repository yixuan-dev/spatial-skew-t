#########################################################################
# merge_score_caches.R -- merge scores<setting>_nonsta score caches along
# the dataset dimension (n=10 cache + n=30 chunk caches -> one n=30 cache).
#
# Usage (from code/analysis/simstudy):
#   Rscript DESKTOP-61SBCCI/merge_score_caches.R <out.RData> <in1.RData> <in2.RData> ...
#
# Object handling (classified dynamically, not by a hardcoded name list):
#   - arrays with a named "dataset" dimension -> bound along that dimension,
#     then sorted by dataset id
#   - `datasets` integer vector              -> concatenated + sorted
#   - `data_path`                            -> kept from the first cache
#     (machine-specific absolute path)
#   - everything else (setting, methods, mrts_ks, probs, intervals, vs_p,
#     es_max_draws, data_suffix, fits_dir, ...) -> must be identical across
#     all inputs, else the merge aborts
#
# Schema evolution: a dataset-array present in only SOME inputs (e.g.
# `crps.score`, added to scores.R after the setting-16 n=10 fits were
# already deleted) is NA-filled for the datasets of the inputs that lack
# it. A NON-dataset object missing from some inputs is still an error.
#
# Built-in regression check: after merging, every input cache's cells must
# be reproduced exactly (tolerance 0) in the merged output; on failure the
# output file is removed and the script exits non-zero. With the n=10
# cache as an input this IS the "merged 1:10 == n10.bak" gate.
#########################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("usage: merge_score_caches.R <out.RData> <in1.RData> <in2.RData> ...",
       call. = FALSE)
}
out_file <- args[1]
in_files <- args[-1]
for (f in in_files) if (!file.exists(f)) stop(sprintf("input not found: %s", f), call. = FALSE)

envs <- lapply(in_files, function(f) { e <- new.env(); load(f, envir = e); e })

nms <- sort(unique(unlist(lapply(envs, ls))))

is_dataset_array <- function(o) {
  is.array(o) && !is.null(names(dimnames(o))) && "dataset" %in% names(dimnames(o))
}

# which inputs carry each object; partial presence is allowed only for
# dataset-arrays (schema evolution, e.g. crps.score) -- they get NA-filled
have_idx <- lapply(nms, function(nm) {
  which(vapply(envs, function(e) exists(nm, envir = e, inherits = FALSE), TRUE))
})
names(have_idx) <- nms
for (nm in nms) {
  hv <- have_idx[[nm]]
  if (length(hv) == length(envs)) next
  o1 <- get(nm, envir = envs[[hv[1]]])
  if (!is_dataset_array(o1)) {
    stop(sprintf("non-dataset object `%s` is missing from: %s",
                 nm, paste(in_files[-hv], collapse = ", ")), call. = FALSE)
  }
  cat(sprintf("note: `%s` absent from %d input(s) -- NA-filled for their datasets\n",
              nm, length(envs) - length(hv)))
}

# expand a dataset-array to an env that lacks it: same dims, NA values,
# dataset dimnames = that env's datasets
na_like <- function(template, ds) {
  along <- which(names(dimnames(template)) == "dataset")
  d <- dim(template); d[along] <- length(ds)
  dn <- dimnames(template); dn[[along]] <- as.character(ds)
  array(NA_real_, dim = d, dimnames = dn)
}

bind_dataset_arrays <- function(arrs, nm) {
  dn1 <- dimnames(arrs[[1]])
  along <- which(names(dn1) == "dataset")
  nd <- length(dim(arrs[[1]]))
  perm <- c(along, setdiff(seq_len(nd), along))
  ap_list <- lapply(arrs, aperm, perm = perm)
  ref_rest <- dimnames(ap_list[[1]])[-1]
  for (i in seq_along(ap_list)[-1]) {
    if (!identical(dimnames(ap_list[[i]])[-1], ref_rest)) {
      stop(sprintf("%s: non-dataset dimnames differ between inputs 1 and %d", nm, i),
           call. = FALSE)
    }
  }
  ds_names <- unlist(lapply(ap_list, function(a) dimnames(a)[[1]]))
  if (anyDuplicated(ds_names)) {
    stop(sprintf("%s: duplicate dataset ids across inputs: %s", nm,
                 paste(unique(ds_names[duplicated(ds_names)]), collapse = ",")),
         call. = FALSE)
  }
  big <- do.call(rbind, lapply(ap_list, function(a) matrix(a, nrow = dim(a)[1])))
  ord <- order(as.integer(ds_names))
  big <- big[ord, , drop = FALSE]
  out <- array(big, dim = c(length(ds_names), dim(ap_list[[1]])[-1]))
  dimnames(out) <- c(list(dataset = ds_names[ord]), ref_rest)
  aperm(out, order(perm))   # order(perm) is the inverse permutation
}

merged <- new.env()
for (nm in nms) {
  hv <- have_idx[[nm]]
  o1 <- get(nm, envir = envs[[hv[1]]])
  if (is_dataset_array(o1)) {
    objs <- lapply(seq_along(envs), function(i) {
      if (i %in% hv) get(nm, envir = envs[[i]])
      else na_like(o1, get("datasets", envir = envs[[i]]))
    })
    assign(nm, bind_dataset_arrays(objs, nm), envir = merged)
  } else if (nm == "datasets") {
    objs <- lapply(envs, get, x = nm)
    ds <- unlist(objs)
    if (anyDuplicated(ds)) stop("`datasets` overlap across inputs", call. = FALSE)
    assign(nm, sort(as.integer(ds)), envir = merged)
  } else if (nm == "data_path") {
    assign(nm, o1, envir = merged)
  } else {
    for (i in hv[-1]) {
      oi <- get(nm, envir = envs[[i]])
      if (!isTRUE(all.equal(o1, oi, tolerance = 0))) {
        stop(sprintf("metadata `%s` differs between %s and %s -- refusing to merge",
                     nm, in_files[hv[1]], in_files[i]), call. = FALSE)
      }
    }
    assign(nm, o1, envir = merged)
  }
}

# consistency: datasets vector == dataset dimnames of every merged array
ds_vec <- as.character(get("datasets", envir = merged))
for (nm in nms) {
  o <- get(nm, envir = merged)
  if (is_dataset_array(o)) {
    if (!identical(dimnames(o)[["dataset"]], ds_vec)) {
      stop(sprintf("`%s` dataset dimnames disagree with `datasets` vector", nm),
           call. = FALSE)
    }
  }
}

# regression: every input cache must be reproduced exactly in the merge
subset_by_datasets <- function(o, ds) {
  along <- which(names(dimnames(o)) == "dataset")
  idx <- rep(list(quote(expr = )), length(dim(o)))
  idx[[along]] <- match(ds, dimnames(o)[[along]])
  do.call(`[`, c(list(o), idx, list(drop = FALSE)))
}
for (i in seq_along(envs)) {
  ds_i <- as.character(get("datasets", envir = envs[[i]]))
  for (nm in nms) {
    if (!(i %in% have_idx[[nm]])) next
    o_in <- get(nm, envir = envs[[i]])
    if (!is_dataset_array(o_in)) next
    o_sub <- subset_by_datasets(get(nm, envir = merged), ds_i)
    ok <- isTRUE(all.equal(as.vector(o_sub), as.vector(o_in), tolerance = 0)) &&
      identical(dimnames(o_sub)[["dataset"]], ds_i)
    if (!ok) {
      if (file.exists(out_file)) file.remove(out_file)
      stop(sprintf("REGRESSION FAILED: merged `%s` does not reproduce %s", nm, in_files[i]),
           call. = FALSE)
    }
  }
}
cat(sprintf("regression: all %d inputs reproduced exactly in the merge\n", length(envs)))

save(list = nms, envir = merged, file = out_file)
cat(sprintf("merged %d caches -> %s (datasets %s, n=%d)\n",
            length(in_files), out_file,
            paste(range(get("datasets", envir = merged)), collapse = ".."),
            length(get("datasets", envir = merged))))
