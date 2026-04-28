# autoselect.R
#
# Scores each candidate setting on the validation split by loading the
# corresponding fit from fits/<setting>.RData (produced by run_settings_val.R),
# computing Brier score against Y[val_sites, ], and selecting the best setting.
#
# How to run (PowerShell):
#   cd D:\Github\spatial-skew-t\code\analysis\ozone\US-all-auto
#   $env:US_ALL_AUTOSELECT_SETTINGS    = '111,112,204:206'
#   $env:US_ALL_AUTOSELECT_TARGET_PROBS = '0.90,0.95,0.98,0.99,0.995'
#   $env:US_ALL_AUTOSELECT_OBJECTIVE   = 'each'
#   Rscript autoselect.R
#
# Prerequisites:
#   Rscript us-all-setup-auto.R    (produces us-all-setup-auto.RData)
#   Rscript run_settings_val.R     (produces fits/val-<setting>.RData)
#
# Environment variables:
#   US_ALL_AUTOSELECT_SETTINGS     Range tokens, e.g. "111,112,204:206"  [required]
#   US_ALL_AUTOSELECT_TARGET_PROBS Comma-separated probabilities          [default "0.90,0.95,0.98,0.99,0.995"]
#   US_ALL_AUTOSELECT_OBJECTIVE    "each" or "mean"        [default "each"]
#   US_ALL_SUMMARY_DRAWS           Posterior draw cap      [default 5000]
#   US_ALL_VAL_RESULTS_DIR         Override fits/ dir      [default: auto]
#   US_ALL_AUTOSELECT_SETUP_FILE   Path to setup .RData   [default: auto]

rm(list = ls())

# ---- Locate this script ----

.this_script_dir <- local({
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", script_arg[1]),
      winslash = "/", mustWork = FALSE
    ))
  } else {
    tryCatch(
      dirname(normalizePath(sys.frame(1)$ofile,
        winslash = "/", mustWork = FALSE
      )),
      error = function(e) "."
    )
  }
})
setwd(.this_script_dir)

.us_all_auto_root <- normalizePath(.this_script_dir, winslash = "/", mustWork = FALSE)
.r_root <- normalizePath(file.path(.us_all_auto_root, "../../../R"),
  winslash = "/", mustWork = FALSE
)

# Source BrierScore — identical in ar2/ and R/; prefer ar2/ when present.
source(local({
  ar2 <- file.path(.r_root, "ar2", "auxfunctions.R")
  if (file.exists(ar2)) ar2 else file.path(.r_root, "auxfunctions.R")
}))

# ---- Helpers ----

parse_env_int <- function(name, default) {
  val <- suppressWarnings(as.integer(trimws(Sys.getenv(name, unset = ""))))
  if (length(val) == 1L && is.finite(val) && val > 0L) val else as.integer(default)
}

expand_range_tokens <- function(raw_str) {
  tokens <- unlist(strsplit(raw_str, ",", fixed = TRUE))
  out <- integer(0L)
  for (tok in tokens) {
    tok <- trimws(tok)
    if (!nzchar(tok)) next
    if (grepl("^[0-9]+:[0-9]+$", tok)) {
      parts <- as.integer(strsplit(tok, ":", fixed = TRUE)[[1]])
      if (length(parts) == 2L && !anyNA(parts)) {
        step <- if (parts[1L] <= parts[2L]) 1L else -1L
        out <- c(out, seq.int(parts[1L], parts[2L], by = step))
      }
    } else {
      val <- suppressWarnings(as.integer(tok))
      if (!is.na(val)) out <- c(out, val)
    }
  }
  sort(unique(out))
}

# ---- Parse environment variables ----

settings_raw <- trimws(Sys.getenv("US_ALL_AUTOSELECT_SETTINGS", unset = ""))
if (!nzchar(settings_raw)) {
  stop(
    "US_ALL_AUTOSELECT_SETTINGS must be set.\n",
    "Example:  $env:US_ALL_AUTOSELECT_SETTINGS = '111,112,204:206'"
  )
}
candidate_settings <- expand_range_tokens(settings_raw)
if (length(candidate_settings) == 0L) {
  stop("No valid settings parsed from US_ALL_AUTOSELECT_SETTINGS='", settings_raw, "'")
}

target_probs <- suppressWarnings(
  as.numeric(unlist(strsplit(
    trimws(Sys.getenv("US_ALL_AUTOSELECT_TARGET_PROBS", unset = "0.90,0.95,0.98,0.99,0.995")),
    ",",
    fixed = TRUE
  )))
)
target_probs <- sort(unique(
  target_probs[is.finite(target_probs) & target_probs > 0 & target_probs < 1]
))
if (length(target_probs) == 0L) {
  stop("No valid probabilities in US_ALL_AUTOSELECT_TARGET_PROBS.")
}

objective <- trimws(Sys.getenv("US_ALL_AUTOSELECT_OBJECTIVE", unset = "each"))
if (!objective %in% c("each", "mean")) {
  message("Unknown objective '", objective, "'; falling back to 'each'.")
  objective <- "each"
}

summary_draws <- parse_env_int("US_ALL_SUMMARY_DRAWS", 5000L)
fits_dir_override <- trimws(Sys.getenv("US_ALL_VAL_RESULTS_DIR", unset = ""))
fits_dir <- if (nzchar(fits_dir_override)) {
  normalizePath(fits_dir_override, winslash = "/", mustWork = FALSE)
} else {
  file.path(.us_all_auto_root, "fits")
}

# ---- Load setup ----

setup_file_override <- trimws(Sys.getenv("US_ALL_AUTOSELECT_SETUP_FILE", unset = ""))
setup_path <- if (nzchar(setup_file_override)) {
  normalizePath(setup_file_override, winslash = "/", mustWork = FALSE)
} else {
  file.path(.us_all_auto_root, "us-all-setup-auto.RData")
}
if (!file.exists(setup_path)) {
  stop(
    "Setup file not found: ", setup_path,
    "\nRun:  Rscript us-all-setup-auto.R"
  )
}
setup_env <- new.env(parent = emptyenv())
load(setup_path, envir = setup_env)

required <- c("Y", "split.lst")
missing <- required[!vapply(required, exists, logical(1),
  envir = setup_env, inherits = FALSE
)]
if (length(missing) > 0L) {
  stop(
    "us-all-setup-auto.RData is missing: ", paste(missing, collapse = ", "),
    "\nRe-run us-all-setup-auto.R."
  )
}

Y <- get("Y", envir = setup_env)
split.lst <- get("split.lst", envir = setup_env)

val_sites <- split.lst$val
n_train <- length(split.lst$train)
n_val <- length(split.lst$val)
n_test <- length(split.lst$test)

# ---- Output directories ----

.out_root <- file.path("output", "us-all-auto")
for (.d in file.path(.out_root, c(".", "results", "tables"))) {
  dir.create(.d, recursive = TRUE, showWarnings = FALSE)
}

out_path <- function(..., subdir) file.path(.out_root, subdir, ...)

# ---- Load settings metadata (for output decoration) ----

settings <- local({
  candidates <- c(
    file.path(.us_all_auto_root, "settings-auto.csv"),
    file.path(.us_all_auto_root, "../US-all/settings.csv")
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0L) {
    return(NULL)
  }
  df <- read.csv(hit[1L], stringsAsFactors = FALSE)
  if ("setting" %in% names(df) && !("setting_num" %in% names(df))) {
    df$setting_num <- suppressWarnings(as.integer(df$setting))
  }
  df
})

cat("=== US-all-auto: autoselect ===\n")
cat(sprintf("Setup              : %s\n", setup_path))
cat(sprintf("Fits directory     : %s\n", fits_dir))
cat(sprintf(
  "Site partition     : %d train | %d val | %d test\n",
  n_train, n_val, n_test
))
cat(sprintf(
  "Candidate settings : %d  [%s]\n",
  length(candidate_settings), paste(candidate_settings, collapse = ", ")
))
cat(sprintf("Target probs       : %s\n", paste(target_probs, collapse = ", ")))
cat(sprintf("Objective          : %s\n", objective))
cat(sprintf("Summary draws      : %d\n\n", summary_draws))

# ---- Compute exceedance thresholds from the full Y ----
# Use all rows of Y (all 800 sites x all time) to avoid split bias.

thresholds <- quantile(Y, probs = target_probs, na.rm = TRUE)
cat("Exceedance thresholds (full dataset):\n")
for (k in seq_along(target_probs)) {
  cat(sprintf("  p = %.4f  ->  %.4f\n", target_probs[k], thresholds[k]))
}
cat("\n")

Y_val <- Y[val_sites, ] # n_val x ntime — observed values at val sites

# ---- Score all candidate settings ----

cat("Scoring", length(candidate_settings), "settings on", n_val, "val sites...\n")

val_brier <- matrix(NA_real_,
  nrow = length(candidate_settings), ncol = length(target_probs),
  dimnames = list(
    as.character(candidate_settings),
    as.character(target_probs)
  )
)
val_status <- character(length(candidate_settings))
names(val_status) <- as.character(candidate_settings)
val_elapsed <- rep(NA_real_, length(candidate_settings))
names(val_elapsed) <- as.character(candidate_settings)

for (ii in seq_along(candidate_settings)) {
  sid <- candidate_settings[ii]
  fpath <- file.path(fits_dir, sprintf("val-%d.RData", sid))

  cat(sprintf(
    "  [%d/%d] setting %d: %s\n",
    ii, length(candidate_settings), sid, basename(fpath)
  ))

  if (!file.exists(fpath)) {
    val_status[ii] <- "missing_file"
    next
  }

  load_env <- new.env(parent = emptyenv())
  ok <- isTRUE(tryCatch(
    {
      load(fpath, envir = load_env)
      TRUE
    },
    error = function(e) FALSE
  ))
  if (!ok) {
    val_status[ii] <- "load_error"
    next
  }

  if (!exists("fit", envir = load_env, inherits = FALSE)) {
    val_status[ii] <- "no_fit_object"
    next
  }

  fit <- get("fit", envir = load_env, inherits = FALSE)
  if (exists("runtime_info", envir = load_env, inherits = FALSE)) {
    ri <- get("runtime_info", envir = load_env, inherits = FALSE)
    if (is.numeric(ri$elapsed_sec) && length(ri$elapsed_sec) == 1L) {
      val_elapsed[ii] <- ri$elapsed_sec
    }
  }
  if (is.null(fit$yp) || length(dim(fit$yp)) != 3L) {
    val_status[ii] <- "bad_yp"
    next
  }

  result <- tryCatch(
    {
      n_draws <- dim(fit$yp)[1L]
      draw_idx <- if (n_draws <= summary_draws) {
        seq_len(n_draws)
      } else {
        sort(sample.int(n_draws, summary_draws))
      }
      BrierScore(
        preds = fit$yp[draw_idx, , , drop = FALSE],
        thresholds = thresholds, validate = Y_val
      )
    },
    error = function(e) {
      message("  Scoring error for setting ", sid, ": ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(result)) {
    val_status[ii] <- "scoring_error"
  } else {
    val_brier[ii, ] <- result
    val_status[ii] <- "ok"
  }
}

# ---- Select best setting(s) ----

top_n_report <- 4L

if (objective == "each") {
  top_by_prob <- matrix(NA_integer_,
    nrow = length(target_probs), ncol = top_n_report,
    dimnames = list(
      as.character(target_probs),
      paste0("rank", seq_len(top_n_report))
    )
  )
  n_valid_by_prob <- integer(length(target_probs))
  for (k in seq_along(target_probs)) {
    scores_k <- val_brier[, k]
    valid_k <- is.finite(scores_k)
    n_valid_by_prob[k] <- sum(valid_k)
    if (any(valid_k)) {
      ord_k <- order(scores_k[valid_k])
      take <- seq_len(min(top_n_report, sum(valid_k)))
      top_by_prob[k, take] <- candidate_settings[valid_k][ord_k[take]]
    }
  }
  best_by_prob <- top_by_prob[, 1L]
  if (all(is.na(best_by_prob))) {
    stop(
      "No candidate settings produced valid Brier scores.\n",
      "Check that fits/val-<setting>.RData files exist and contain fit$yp."
    )
  }
  cat("\n=== Validation selection result (each) ===\n")
  for (k in seq_along(target_probs)) {
    cat(sprintf(
      "  p = %.4f  (valid: %d settings)\n",
      target_probs[k], n_valid_by_prob[k]
    ))
    for (r in seq_len(top_n_report)) {
      sid <- top_by_prob[k, r]
      if (is.na(sid)) next
      cat(sprintf(
        "    rank %d: setting %d  (Brier = %.6f)\n",
        r, sid, val_brier[as.character(sid), k]
      ))
    }
  }
  cat("\n")
} else {
  selection_score <- rowMeans(val_brier, na.rm = TRUE)
  valid_mask <- is.finite(selection_score)
  if (!any(valid_mask)) {
    stop(
      "No candidate settings produced valid Brier scores.\n",
      "Check that fits/val-<setting>.RData files exist and contain fit$yp."
    )
  }
  valid_ord <- order(selection_score[valid_mask])
  valid_settings <- candidate_settings[valid_mask]
  take <- seq_len(min(top_n_report, sum(valid_mask)))
  top_settings <- valid_settings[valid_ord[take]]
  best_setting <- top_settings[1L]
  best_val_brier <- val_brier[as.character(best_setting), , drop = TRUE]
  cat("\n=== Validation selection result (mean) ===\n")
  cat(sprintf("Candidates scored: %d / %d\n", sum(valid_mask), length(candidate_settings)))
  for (r in seq_along(top_settings)) {
    sid <- top_settings[r]
    score <- selection_score[which(candidate_settings == sid)]
    cat(sprintf("  rank %d: setting %d  (mean Brier = %.6f)\n", r, sid, score))
  }
  cat("\n")
}

# ---- Build output table ----

if (objective == "each") {
  rank_ref <- val_brier[, 1L]
  val_df <- data.frame(
    setting = candidate_settings,
    status = val_status,
    rank = rank(ifelse(is.finite(rank_ref), rank_ref, Inf), ties.method = "min"),
    elapsed_sec = val_elapsed,
    stringsAsFactors = FALSE
  )
  for (k in seq_along(target_probs)) {
    p_tag <- sprintf("p%.4f", target_probs[k])
    scores_k <- val_brier[, k]
    val_df[[paste0("val_brier_", p_tag)]] <- scores_k
    val_df[[paste0("rank_", p_tag)]] <- ifelse(
      is.finite(scores_k),
      rank(ifelse(is.finite(scores_k), scores_k, Inf), ties.method = "min"),
      NA_integer_
    )
  }
} else {
  val_df <- data.frame(
    setting = candidate_settings,
    status = val_status,
    rank = rank(ifelse(valid_mask, selection_score, Inf), ties.method = "min"),
    elapsed_sec = val_elapsed,
    selection_score = selection_score,
    is_best = candidate_settings == best_setting,
    stringsAsFactors = FALSE
  )
  for (k in seq_along(target_probs)) {
    val_df[[sprintf("val_brier_p%.4f", target_probs[k])]] <- val_brier[, k]
  }
}

if (!is.null(settings) && "setting_num" %in% names(settings)) {
  meta_idx <- match(val_df$setting, settings$setting_num)
  for (col in setdiff(names(settings), c("setting", "setting_num"))) {
    val_df[[col]] <- settings[[col]][meta_idx]
  }
}
val_df <- val_df[order(val_df$rank), ]

# ---- Write outputs ----

write.csv(val_df,
  file = out_path("autoselect_validation_scores_200-200.csv", subdir = "tables"),
  row.names = FALSE
)

save_vars <- c(
  "candidate_settings", "target_probs", "thresholds", "objective",
  "n_train", "n_val", "val_sites", "split.lst",
  "val_brier", "val_status", "val_df", "fits_dir", "summary_draws"
)
save_vars <- c(
  save_vars,
  if (objective == "each") {
    c("top_by_prob", "best_by_prob", "n_valid_by_prob")
  } else {
    c("selection_score", "valid_mask", "top_settings", "best_setting", "best_val_brier")
  }
)
save(
  list = save_vars,
  file = out_path("us-all-auto-results-200-200.RData", subdir = "results")
)

cat("Outputs written:\n")
cat(" -", out_path("autoselect_validation_scores_200-200.csv", subdir = "tables"), "\n")
cat(" -", out_path("us-all-auto-results-200-200.RData", subdir = "results"), "\n")
