rm(list = ls())
library(compiler)
enableJIT(3)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
    script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = FALSE)
    script_dir <- dirname(script_path)
    if (dir.exists(script_dir)) {
        setwd(script_dir)
    }
}

load("us-all-setup.RData")
source("../../../R/auxfunctions.R")
settings <- read.csv("settings.csv", stringsAsFactors = FALSE)
settings$setting_num <- suppressWarnings(as.integer(settings$setting))

# Keep the Morris baseline frozen for reproducibility.
done_morris <- c(1:5, 7:9, 11:13, 15:17, 33:36, 38:41, 43:46, 51:74)

# Proposed model namespace (AR2 lane). The settings file already marks these rows.
done_proposed <- settings$setting_num[tolower(settings$ar2) == "yes"]
done_proposed <- sort(unique(done_proposed[!is.na(done_proposed)]))

all_requested <- sort(unique(c(done_morris, done_proposed)))

# Use exactly the same scoring basis as us-all-results.R.
probs <- c(0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995)
thresholds <- quantile(Y, probs = probs, na.rm = TRUE)
nsets <- length(cv.lst)

max_setting <- max(all_requested)
quant.score <- array(NA_real_, dim = c(length(probs), nsets, max_setting))
brier.score <- array(NA_real_, dim = c(length(thresholds), nsets, max_setting))

available_settings <- integer(0)
skipped_missing_file <- integer(0)
skipped_bad_contract <- integer(0)
skipped_scoring_error <- integer(0)

validate_pred_contract <- function(pred, validate) {
    if (is.null(pred)) {
        return(FALSE)
    }
    if (length(dim(pred)) != 3) {
        return(FALSE)
    }
    pred_dims <- dim(pred)
    val_dims <- dim(validate)

    # Expected contract from existing scorer usage:
    # pred: [iter, n_validation_sites, n_time]
    # validate: [n_validation_sites, n_time]
    if (length(val_dims) != 2) {
        return(FALSE)
    }
    if (pred_dims[2] != val_dims[1] || pred_dims[3] != val_dims[2]) {
        return(FALSE)
    }
    TRUE
}

for (i in all_requested) {
    file <- sprintf("results/us-all-%d.RData", i)
    cat("start file", file, "\n")

    if (!file.exists(file)) {
        skipped_missing_file <- c(skipped_missing_file, i)
        cat("skip file (missing):", file, "\n")
        next
    }

    load(file)
    if (!exists("fit")) {
        skipped_bad_contract <- c(skipped_bad_contract, i)
        cat("skip setting", i, ": object 'fit' not found\n")
        next
    }

    has_bad_fold <- FALSE
    has_scoring_error <- FALSE
    score_error_message <- ""
    for (d in seq_len(nsets)) {
        if (is.null(fit[[d]])) {
            has_bad_fold <- TRUE
            break
        }
        val.idx <- cv.lst[[d]]
        validate <- Y[val.idx, ]
        pred.d <- fit[[d]]$yp[, , ]

        if (i != 2 && !validate_pred_contract(pred.d, validate)) {
            has_bad_fold <- TRUE
            break
        }

        trans <- (i == 2)
        score_try <- tryCatch(
            {
                list(
                    quant = QuantScore(pred.d, probs, validate, trans = trans),
                    brier = BrierScore(pred.d, thresholds, validate, trans = trans)
                )
            },
            error = function(e) {
                score_error_message <<- conditionMessage(e)
                NULL
            }
        )

        if (is.null(score_try)) {
            has_scoring_error <- TRUE
            break
        }

        quant.score[, d, i] <- score_try$quant
        brier.score[, d, i] <- score_try$brier
    }

    if (has_scoring_error) {
        skipped_scoring_error <- c(skipped_scoring_error, i)
        cat("skip setting", i, ": scoring error ->", score_error_message, "\n")
        next
    }

    if (has_bad_fold) {
        skipped_bad_contract <- c(skipped_bad_contract, i)
        cat("skip setting", i, ": fit[[d]]$yp contract mismatch\n")
        next
    }

    available_settings <- c(available_settings, i)
    cat("finish file", file, "\n")
}

available_settings <- sort(unique(available_settings))
skipped_missing_file <- sort(unique(skipped_missing_file))
skipped_bad_contract <- sort(unique(skipped_bad_contract))
skipped_scoring_error <- sort(unique(skipped_scoring_error))

if (!1 %in% available_settings) {
    stop("Gaussian baseline setting 1 is missing; cannot compute relative scores.")
}

quant.score.mean <- matrix(NA_real_, nrow = max_setting, ncol = length(probs))
brier.score.mean <- matrix(NA_real_, nrow = max_setting, ncol = length(thresholds))
quant.score.se <- matrix(NA_real_, nrow = max_setting, ncol = length(probs))
brier.score.se <- matrix(NA_real_, nrow = max_setting, ncol = length(thresholds))

for (i in available_settings) {
    quant.score.mean[i, ] <- apply(quant.score[, , i, drop = FALSE], 1, mean, na.rm = TRUE)
    quant.score.se[i, ] <- apply(quant.score[, , i, drop = FALSE], 1, sd, na.rm = TRUE) / sqrt(nsets)
    brier.score.mean[i, ] <- apply(brier.score[, , i, drop = FALSE], 1, mean, na.rm = TRUE)
    brier.score.se[i, ] <- apply(brier.score[, , i, drop = FALSE], 1, sd, na.rm = TRUE) / sqrt(nsets)
}

bs.mean.ref.gau <- matrix(NA_real_, nrow = max_setting, ncol = length(thresholds))
qs.mean.ref.gau <- matrix(NA_real_, nrow = max_setting, ncol = length(probs))
for (i in available_settings) {
    bs.mean.ref.gau[i, ] <- brier.score.mean[i, ] / brier.score.mean[1, ]
    qs.mean.ref.gau[i, ] <- quant.score.mean[i, ] / quant.score.mean[1, ]
}

settings_ext <- settings
settings_ext$is_baseline <- settings_ext$setting %in% done_morris
settings_ext$is_proposed <- settings_ext$setting %in% done_proposed
settings_ext$is_available <- settings_ext$setting %in% available_settings

# Full table (one row per setting x quantile x metric)
full_rows <- list()
row_id <- 1
for (i in available_settings) {
    setting_meta <- settings_ext[settings_ext$setting == i, , drop = FALSE]
    for (q in seq_along(probs)) {
        full_rows[[row_id]] <- data.frame(
            setting = i,
            metric = "brier",
            quantile = probs[q],
            score_mean = brier.score.mean[i, q],
            score_se = brier.score.se[i, q],
            rel_to_gaussian = bs.mean.ref.gau[i, q],
            method = if (nrow(setting_meta) > 0) setting_meta$method[1] else NA,
            knots = if (nrow(setting_meta) > 0) setting_meta$knots[1] else NA,
            thresh = if (nrow(setting_meta) > 0) setting_meta$thresh[1] else NA,
            CMAQ = if (nrow(setting_meta) > 0) setting_meta$CMAQ[1] else NA,
            TS = if (nrow(setting_meta) > 0) setting_meta$TS[1] else NA,
            ar2 = if (nrow(setting_meta) > 0) setting_meta$ar2[1] else NA,
            is_baseline = i %in% done_morris,
            is_proposed = i %in% done_proposed,
            stringsAsFactors = FALSE
        )
        row_id <- row_id + 1

        full_rows[[row_id]] <- data.frame(
            setting = i,
            metric = "quantile",
            quantile = probs[q],
            score_mean = quant.score.mean[i, q],
            score_se = quant.score.se[i, q],
            rel_to_gaussian = qs.mean.ref.gau[i, q],
            method = if (nrow(setting_meta) > 0) setting_meta$method[1] else NA,
            knots = if (nrow(setting_meta) > 0) setting_meta$knots[1] else NA,
            thresh = if (nrow(setting_meta) > 0) setting_meta$thresh[1] else NA,
            CMAQ = if (nrow(setting_meta) > 0) setting_meta$CMAQ[1] else NA,
            TS = if (nrow(setting_meta) > 0) setting_meta$TS[1] else NA,
            ar2 = if (nrow(setting_meta) > 0) setting_meta$ar2[1] else NA,
            is_baseline = i %in% done_morris,
            is_proposed = i %in% done_proposed,
            stringsAsFactors = FALSE
        )
        row_id <- row_id + 1
    }
}
comparison_full_table <- do.call(rbind, full_rows)
write.csv(comparison_full_table, "comparison_full_table.csv", row.names = FALSE)

# Top-2 table at target quantiles for brier score.
target_quantiles <- c(0.95, 0.98, 0.99, 0.995)
top2_rows <- list()
row_id <- 1
for (q in target_quantiles) {
    q_idx <- which(abs(probs - q) < 1e-12)
    if (length(q_idx) == 0) {
        next
    }

    rel_vec <- bs.mean.ref.gau[available_settings, q_idx]
    ord <- order(rel_vec, decreasing = FALSE, na.last = NA)
    n_take <- min(2, length(ord))

    if (n_take > 0) {
        winners <- available_settings[ord[seq_len(n_take)]]
        for (rank_idx in seq_along(winners)) {
            w <- winners[rank_idx]
            meta <- settings_ext[settings_ext$setting == w, , drop = FALSE]
            top2_rows[[row_id]] <- data.frame(
                quantile = q,
                rank = rank_idx,
                setting = w,
                rel_brier_to_gaussian = bs.mean.ref.gau[w, q_idx],
                method = if (nrow(meta) > 0) meta$method[1] else NA,
                knots = if (nrow(meta) > 0) meta$knots[1] else NA,
                thresh = if (nrow(meta) > 0) meta$thresh[1] else NA,
                CMAQ = if (nrow(meta) > 0) meta$CMAQ[1] else NA,
                TS = if (nrow(meta) > 0) meta$TS[1] else NA,
                ar2 = if (nrow(meta) > 0) meta$ar2[1] else NA,
                is_baseline = w %in% done_morris,
                is_proposed = w %in% done_proposed,
                stringsAsFactors = FALSE
            )
            row_id <- row_id + 1
        }
    }
}
comparison_top2 <- if (length(top2_rows) > 0) do.call(rbind, top2_rows) else data.frame()
write.csv(comparison_top2, "comparison_top2.csv", row.names = FALSE)

# Same-basis paired table: match proposed rows to baseline rows by method/knots/thresh/CMAQ/TS.
base_meta <- settings_ext[settings_ext$setting %in% done_morris, ]
prop_meta <- settings_ext[settings_ext$setting %in% done_proposed, ]
paired_rows <- list()
row_id <- 1

for (k in seq_len(nrow(prop_meta))) {
    p <- prop_meta[k, ]
    matches <- base_meta[
        base_meta$method == p$method &
            base_meta$knots == p$knots &
            base_meta$thresh == p$thresh &
            base_meta$CMAQ == p$CMAQ &
            base_meta$TS == p$TS,
    ]

    if (nrow(matches) == 0) {
        next
    }

    b_setting <- matches$setting[1]
    p_setting <- p$setting

    if (!(b_setting %in% available_settings && p_setting %in% available_settings)) {
        next
    }

    for (q in seq_along(probs)) {
        paired_rows[[row_id]] <- data.frame(
            quantile = probs[q],
            baseline_setting = b_setting,
            proposed_setting = p_setting,
            baseline_rel_brier = bs.mean.ref.gau[b_setting, q],
            proposed_rel_brier = bs.mean.ref.gau[p_setting, q],
            delta_rel_brier = bs.mean.ref.gau[p_setting, q] - bs.mean.ref.gau[b_setting, q],
            baseline_rel_quant = qs.mean.ref.gau[b_setting, q],
            proposed_rel_quant = qs.mean.ref.gau[p_setting, q],
            delta_rel_quant = qs.mean.ref.gau[p_setting, q] - qs.mean.ref.gau[b_setting, q],
            method = p$method,
            knots = p$knots,
            thresh = p$thresh,
            CMAQ = p$CMAQ,
            TS = p$TS,
            stringsAsFactors = FALSE
        )
        row_id <- row_id + 1
    }
}

comparison_paired_same_basis <- if (length(paired_rows) > 0) do.call(rbind, paired_rows) else data.frame()
write.csv(comparison_paired_same_basis, "comparison_paired_same_basis.csv", row.names = FALSE)

save(
    list = c(
        "done_morris", "done_proposed", "available_settings",
        "skipped_missing_file", "skipped_bad_contract",
        "skipped_scoring_error",
        "probs", "thresholds", "quant.score", "brier.score",
        "quant.score.mean", "brier.score.mean", "quant.score.se", "brier.score.se",
        "bs.mean.ref.gau", "qs.mean.ref.gau",
        "comparison_full_table", "comparison_top2", "comparison_paired_same_basis"
    ),
    file = "us-all-results-proposed.RData"
)

cat("\n=== us-all-results-proposed summary ===\n")
cat("Requested settings:", length(all_requested), "\n")
cat("Available settings:", length(available_settings), "\n")
cat("Missing result files:", length(skipped_missing_file), "\n")
if (length(skipped_missing_file) > 0) {
    cat("  ->", paste(skipped_missing_file, collapse = ", "), "\n")
}
cat("Contract mismatch settings:", length(skipped_bad_contract), "\n")
if (length(skipped_bad_contract) > 0) {
    cat("  ->", paste(skipped_bad_contract, collapse = ", "), "\n")
}
cat("Scoring error settings:", length(skipped_scoring_error), "\n")
if (length(skipped_scoring_error) > 0) {
    cat("  ->", paste(skipped_scoring_error, collapse = ", "), "\n")
}
cat("Outputs written:\n")
cat("- us-all-results-proposed.RData\n")
cat("- comparison_full_table.csv\n")
cat("- comparison_top2.csv\n")
cat("- comparison_paired_same_basis.csv\n")
