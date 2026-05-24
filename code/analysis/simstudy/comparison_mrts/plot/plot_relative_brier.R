rm(list = ls())

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
    script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = FALSE)
    script_dir <- dirname(script_path)
    if (dir.exists(script_dir)) {
        setwd(script_dir)
    }
}

simstudy_dir <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
setwd(simstudy_dir)

source("../../R/auxfunctions.R")
load("simdata.RData")

parse_index_expr <- function(expr_str, arg_name) {
    expr <- trimws(expr_str)
    if (!nzchar(expr)) {
        stop(sprintf("%s cannot be empty", arg_name), call. = FALSE)
    }
    if (grepl("^\\(.*\\)$", expr)) {
        expr <- paste0("c", expr)
    }

    values <- tryCatch(
        eval(parse(text = expr), envir = baseenv()),
        error = function(e) {
            stop(sprintf("Invalid %s expression '%s'", arg_name, expr_str), call. = FALSE)
        }
    )

    if (!is.numeric(values) || length(values) == 0) {
        stop(sprintf("%s must evaluate to a non-empty numeric vector", arg_name), call. = FALSE)
    }

    values_int <- as.integer(values)
    if (any(is.na(values_int)) || any(values_int != values)) {
        stop(sprintf("%s must contain integers only", arg_name), call. = FALSE)
    }

    sort(unique(values_int))
}

build_result_file <- function(results_dir, setting_id, method_id, dataset_id, mrts_k = NA_integer_) {
    tag <- if (is.na(mrts_k)) "" else paste0("-K", as.integer(mrts_k))
    file.path(results_dir, sprintf("%d-%d-%d%s.RData", as.integer(setting_id), as.integer(method_id), as.integer(dataset_id), tag))
}

extract_fit <- function(result_file) {
    if (!file.exists(result_file)) {
        return(NULL)
    }

    e <- new.env(parent = emptyenv())
    load(result_file, envir = e)

    if (exists("fit.1", envir = e, inherits = FALSE)) {
        return(get("fit.1", envir = e, inherits = FALSE))
    }
    if (exists("fit", envir = e, inherits = FALSE)) {
        return(get("fit", envir = e, inherits = FALSE))
    }

    NULL
}

compute_brier <- function(pred, thresholds, validate) {
    score <- tryCatch(
        BrierScore(pred, thresholds, validate),
        error = function(e) NULL
    )

    if (is.null(score)) {
        score <- tryCatch(
            BrierScore(pred, thresholds, validate, trans = TRUE),
            error = function(e) NULL
        )
    }

    if (is.null(score)) {
        return(rep(NA_real_, length(thresholds)))
    }

    as.numeric(score)
}

resolve_plot_k <- function(method_ids, default_k = 20L) {
    k_override_raw <- trimws(Sys.getenv("SIMSTUDY_MRTS_PLOT_K", unset = ""))
    if (nzchar(k_override_raw)) {
        k_override <- parse_index_expr(k_override_raw, "SIMSTUDY_MRTS_PLOT_K")
        if (length(k_override) != 1 || k_override[1] < 1) {
            stop("SIMSTUDY_MRTS_PLOT_K must be a single positive integer", call. = FALSE)
        }
        return(setNames(rep(as.integer(k_override[1]), length(method_ids)), as.character(method_ids)))
    }

    setNames(rep(as.integer(default_k), length(method_ids)), as.character(method_ids))
}

setting_id <- suppressWarnings(as.integer(Sys.getenv("SIMSTUDY_MRTS_PLOT_SETTING", unset = "4")))
if (!is.finite(setting_id) || is.na(setting_id) || setting_id < 1 || setting_id > dim(y)[4]) {
    stop(sprintf("SIMSTUDY_MRTS_PLOT_SETTING must be in 1..%d", dim(y)[4]), call. = FALSE)
}

results_dir <- trimws(Sys.getenv("SIMSTUDY_MRTS_RESULTS_DIR", unset = "results"))
if (!dir.exists(results_dir)) {
    stop(sprintf("Results directory not found: %s", results_dir), call. = FALSE)
}

datasets_raw <- trimws(Sys.getenv("SIMSTUDY_MRTS_DATASETS", unset = ""))
dataset_ids <- if (nzchar(datasets_raw)) {
    ids <- parse_index_expr(datasets_raw, "SIMSTUDY_MRTS_DATASETS")
    ids[ids >= 1 & ids <= dim(y)[3]]
} else {
    seq_len(dim(y)[3])
}
if (length(dataset_ids) == 0) {
    stop("No valid dataset ids selected", call. = FALSE)
}

target_probs <- c(seq(0, 0.9, by = 0.1), 0.95, 0.98, 0.995)
method_ids <- 2:5
method_label_map <- c(
    "2" = "Skew-t, K=1",
    "3" = "t, K=1, q(0.80)",
    "4" = "Skew-t, K=5",
    "5" = "t, K=5, q(0.80)"
)
mrts_k_by_method <- resolve_plot_k(method_ids = method_ids, default_k = 20L)

obs <- c(rep(TRUE, nrow(y) - ntest), rep(FALSE, ntest))

model_specs <- rbind(
    data.frame(panel = "method 2-5", method_id = method_ids, mrts_k = NA_integer_, stringsAsFactors = FALSE),
    data.frame(panel = "method 2-5 + mrts", method_id = method_ids, mrts_k = as.integer(mrts_k_by_method[as.character(method_ids)]), stringsAsFactors = FALSE)
)

n_probs <- length(target_probs)
n_models <- nrow(model_specs)
n_datasets <- length(dataset_ids)

brier_arr <- array(NA_real_, dim = c(n_probs, n_datasets, n_models))
ref_arr <- array(NA_real_, dim = c(n_probs, n_datasets))

missing_rows <- list()
missing_id <- 1L

for (jj in seq_along(dataset_ids)) {
    dataset_id <- dataset_ids[jj]
    thresholds <- quantile(y[, , dataset_id, setting_id], probs = target_probs, na.rm = TRUE, names = FALSE)
    validate <- y[!obs, , dataset_id, setting_id]

    ref_file <- build_result_file(results_dir, setting_id, method_id = 1L, dataset_id = dataset_id, mrts_k = NA_integer_)
    ref_fit <- extract_fit(ref_file)
    if (is.null(ref_fit) || is.null(ref_fit$yp)) {
        missing_rows[[missing_id]] <- data.frame(
            dataset = dataset_id,
            method_id = 1L,
            mrts_k = NA_integer_,
            result_file = ref_file,
            stringsAsFactors = FALSE
        )
        missing_id <- missing_id + 1L
    } else {
        ref_arr[, jj] <- compute_brier(ref_fit$yp, thresholds, validate)
    }

    for (ii in seq_len(n_models)) {
        spec <- model_specs[ii, , drop = FALSE]
        result_file <- build_result_file(
            results_dir = results_dir,
            setting_id = setting_id,
            method_id = spec$method_id[1],
            dataset_id = dataset_id,
            mrts_k = spec$mrts_k[1]
        )
        fit_obj <- extract_fit(result_file)

        if (is.null(fit_obj) || is.null(fit_obj$yp)) {
            missing_rows[[missing_id]] <- data.frame(
                dataset = dataset_id,
                method_id = spec$method_id[1],
                mrts_k = spec$mrts_k[1],
                result_file = result_file,
                stringsAsFactors = FALSE
            )
            missing_id <- missing_id + 1L
            next
        }

        brier_arr[, jj, ii] <- compute_brier(fit_obj$yp, thresholds, validate)
    }
}

mean_ref <- apply(ref_arr, 1, mean, na.rm = TRUE)
mean_model <- apply(brier_arr, c(1, 3), mean, na.rm = TRUE)
relative_model <- sweep(mean_model, 1, mean_ref, FUN = "/")
relative_model[!is.finite(relative_model)] <- NA_real_

log_relative_model <- suppressWarnings(log(relative_model))
log_relative_model[!is.finite(log_relative_model)] <- NA_real_
nonpositive_rbs_count <- sum(is.finite(relative_model) & relative_model <= 0, na.rm = TRUE)

plot_rows <- list()
row_id <- 1L
for (ii in seq_len(n_models)) {
    spec <- model_specs[ii, , drop = FALSE]
    method_label <- unname(method_label_map[as.character(spec$method_id[1])])
    if (is.na(method_label) || !nzchar(method_label)) {
        method_label <- sprintf("method %d", spec$method_id[1])
    }

    label <- if (is.na(spec$mrts_k[1])) {
        method_label
    } else {
        sprintf("%s + mrts(k=%d)", method_label, spec$mrts_k[1])
    }

    plot_rows[[row_id]] <- data.frame(
        quantile = target_probs,
        relative_brier = relative_model[, ii],
        log_relative_brier = log_relative_model[, ii],
        panel = spec$panel[1],
        method_id = spec$method_id[1],
        mrts_k = spec$mrts_k[1],
        line_label = label,
        stringsAsFactors = FALSE
    )
    row_id <- row_id + 1L
}

plot_df <- do.call(rbind, plot_rows)
plot_df$panel <- factor(plot_df$panel, levels = c("method 2-5", "method 2-5 + mrts"))

output_dir <- file.path("comparison_mrts", "plot")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

data_out <- file.path(output_dir, sprintf("setting%d_relative_brier_plot_data.csv", setting_id))
write.csv(plot_df, data_out, row.names = FALSE)

missing_df <- if (length(missing_rows) > 0) do.call(rbind, missing_rows) else data.frame()
missing_out <- file.path(output_dir, sprintf("setting%d_relative_brier_plot_missing_files.csv", setting_id))
write.csv(missing_df, missing_out, row.names = FALSE)

plot_file <- file.path(output_dir, sprintf("setting%d_relative_brier_method2to5_vs_mrts.png", setting_id))

method_colors <- c(
    "2" = "#D55E00",
    "3" = "#0072B2",
    "4" = "#009E73",
    "5" = "#CC79A7"
)
method_pch <- c("2" = 16, "3" = 17, "4" = 15, "5" = 18)

png(filename = plot_file, width = 1800, height = 2200, res = 220)
par(mfrow = c(2, 1), mar = c(5.2, 6.2, 3.8, 1.5), oma = c(0, 0, 2, 0))

for (panel_name in levels(plot_df$panel)) {
    panel_df <- plot_df[plot_df$panel == panel_name, , drop = FALSE]
    panel_df <- panel_df[order(panel_df$method_id, panel_df$quantile), , drop = FALSE]

    y_vals <- panel_df$log_relative_brier[is.finite(panel_df$log_relative_brier)]
    if (length(y_vals) == 0) {
        y_vals <- c(0)
    }
    y_lim <- range(c(y_vals, 0), na.rm = TRUE)
    y_pad <- 0.05 * diff(y_lim)
    if (!is.finite(y_pad) || y_pad == 0) {
        y_pad <- 0.05
    }
    y_lim <- c(y_lim[1] - y_pad, y_lim[2] + y_pad)

    first_method <- sort(unique(panel_df$method_id))[1]
    first_df <- panel_df[panel_df$method_id == first_method, , drop = FALSE]

    plot(
        x = first_df$quantile,
        y = first_df$log_relative_brier,
        type = "o",
        col = method_colors[as.character(first_method)],
        pch = method_pch[as.character(first_method)],
        lwd = 2,
        cex = 1.1,
        ylim = y_lim,
        xaxt = "n",
        xlab = "Threshold quantile",
        ylab = "log(Relative Brier score) (vs method 1)",
        main = panel_name
    )

    axis(1, at = target_probs, labels = format(target_probs, trim = TRUE, scientific = FALSE), cex.axis = 0.9)
    abline(h = 0, lty = 2, col = "gray40")

    all_methods <- sort(unique(panel_df$method_id))
    if (length(all_methods) > 1) {
        for (method_id in all_methods[-1]) {
            method_df <- panel_df[panel_df$method_id == method_id, , drop = FALSE]
            lines(
                x = method_df$quantile,
                y = method_df$log_relative_brier,
                type = "o",
                col = method_colors[as.character(method_id)],
                pch = method_pch[as.character(method_id)],
                lwd = 2,
                cex = 1.1
            )
        }
    }

    legend_df <- panel_df[!duplicated(panel_df$line_label), c("line_label", "method_id"), drop = FALSE]
    legend_df <- legend_df[order(legend_df$method_id), , drop = FALSE]

    legend(
        "topright",
        legend = legend_df$line_label,
        col = method_colors[as.character(legend_df$method_id)],
        pch = method_pch[as.character(legend_df$method_id)],
        lty = 1,
        lwd = 2,
        bty = "n",
        cex = 0.95
    )
}

mtext(sprintf("Setting %d: log(Relative Brier Score) Profiles", setting_id), outer = TRUE, cex = 1.3, font = 2)
dev.off()

cat("Plot saved:\n")
cat("- ", plot_file, "\n", sep = "")
cat("Data saved:\n")
cat("- ", data_out, "\n", sep = "")
cat("Missing-file log:\n")
cat("- ", missing_out, "\n", sep = "")
if (nonpositive_rbs_count > 0) {
    cat(sprintf("Warning: %d non-positive RBS values encountered; log(RBS) set to NA for those points.\n", nonpositive_rbs_count))
}
cat("MRTS k by method:\n")
for (method_id in method_ids) {
    cat(sprintf("- method %d: k=%d\n", method_id, as.integer(mrts_k_by_method[as.character(method_id)])))
}
cat(sprintf("Datasets used: %d\n", length(dataset_ids)))
