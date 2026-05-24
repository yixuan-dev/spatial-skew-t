rm(list = ls())

# =====================================================================
# Relative-Brier plot (method 2..5, baseline vs +MRTS), side by side.
#
# Reads pre-computed Brier scores from
#   <simstudy>/scores<setting>_<K>mrts.RData
# instead of re-running BrierScore on the fit objects. This avoids the
# data-alignment hazard that previously appeared when simdata.RData
# was regenerated after the fits/scores were cached.
#
#   * CLI override:        Rscript plot_relative_brier.R --setting=5 --k=15
#   * Env-var override:    SIMSTUDY_MRTS_PLOT_K=10 Rscript ...
#
# Required caches per (setting, K):
#   scores<setting>_0mrts.RData    (baseline, gives left panel + Gaussian ref)
#   scores<setting>_<K>mrts.RData  (MRTS K cache, gives right panel)
# =====================================================================
DEFAULT_MRTS_K     <- 20L
DEFAULT_SETTING_ID <- 4L
DEFAULT_SCORES_DIR <- "."   # relative to <simstudy>/

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
    script_path <- normalizePath(sub("^--file=", "", script_arg[1]),
                                 winslash = "/", mustWork = FALSE)
    script_dir <- dirname(script_path)
    if (dir.exists(script_dir)) setwd(script_dir)
}

# ---- CLI flag parsing (--key=value) ---------------------------------
cli_flags <- list()
for (arg in commandArgs(trailingOnly = TRUE)) {
    m <- regmatches(arg, regexec("^--([^=]+)=(.*)$", arg))[[1]]
    if (length(m) == 3L) cli_flags[[m[2]]] <- m[3]
}
get_cli <- function(key) {
    v <- cli_flags[[key]]
    if (is.null(v)) "" else trimws(v)
}

simstudy_dir <- normalizePath(file.path("..", ".."),
                              winslash = "/", mustWork = TRUE)
setwd(simstudy_dir)

parse_pos_int <- function(raw, src) {
    val <- suppressWarnings(as.integer(raw))
    if (!is.finite(val) || is.na(val) || val < 1L || val != as.numeric(raw)) {
        stop(sprintf("%s must be a single positive integer (got '%s')",
                     src, raw), call. = FALSE)
    }
    val
}

# ---- resolve setting / k / scores_dir -------------------------------
setting_raw <- get_cli("setting")
if (!nzchar(setting_raw)) {
    setting_raw <- Sys.getenv("SIMSTUDY_MRTS_PLOT_SETTING",
                              unset = as.character(DEFAULT_SETTING_ID))
}
setting_id <- parse_pos_int(setting_raw, "setting")

k_raw <- get_cli("k")
k_src <- "--k"
if (!nzchar(k_raw)) {
    k_raw <- trimws(Sys.getenv("SIMSTUDY_MRTS_PLOT_K", unset = ""))
    k_src <- "SIMSTUDY_MRTS_PLOT_K"
}
mrts_k_active <- if (nzchar(k_raw)) parse_pos_int(k_raw, k_src) else DEFAULT_MRTS_K

scores_dir <- get_cli("scores-dir")
if (!nzchar(scores_dir)) scores_dir <- DEFAULT_SCORES_DIR
if (!dir.exists(scores_dir)) {
    stop(sprintf("Scores dir not found: %s", scores_dir), call. = FALSE)
}

cat(sprintf("Setting = %d   MRTS k = %d   scores_dir = %s\n",
            setting_id, mrts_k_active, scores_dir))

# ---- load scores cache ----------------------------------------------
load_scores_cache <- function(setting_id, mrts_k, scores_dir) {
    fname <- sprintf("scores%d_%dmrts.RData", setting_id, mrts_k)
    fpath <- file.path(scores_dir, fname)
    if (!file.exists(fpath)) {
        stop(sprintf("Scores cache not found: %s", fpath), call. = FALSE)
    }
    e <- new.env(parent = emptyenv())
    load(fpath, envir = e)
    needed <- c("brier.score", "probs", "methods", "datasets")
    miss <- setdiff(needed, ls(e))
    if (length(miss) > 0L) {
        stop(sprintf("Cache %s missing objects: %s",
                     fpath, paste(miss, collapse = ", ")), call. = FALSE)
    }
    list(
        file        = fpath,
        brier.score = e$brier.score,
        probs       = e$probs,
        methods     = e$methods,
        datasets    = e$datasets,
        mrts_k      = if (exists("mrts_k", envir = e)) e$mrts_k else mrts_k,
        setting     = if (exists("setting", envir = e)) e$setting else setting_id
    )
}

cache_base <- load_scores_cache(setting_id, 0L, scores_dir)
cache_mrts <- load_scores_cache(setting_id, mrts_k_active, scores_dir)
cat(sprintf("Loaded: %s\n", cache_base$file))
cat(sprintf("Loaded: %s\n", cache_mrts$file))

# ---- pick target quantiles (subset of cache probs, fuzzy match) -----
target_probs <- c(seq(0.90, 0.99, by = 0.01), 0.995)
match_probs <- function(targets, cache_probs, tol = 1e-8) {
    vapply(targets, function(t) {
        idx <- which(abs(cache_probs - t) < tol)
        if (length(idx) == 0L) NA_integer_ else idx[1]
    }, integer(1))
}
prob_idx      <- match_probs(target_probs, cache_base$probs)
prob_idx_mrts <- match_probs(target_probs, cache_mrts$probs)
if (any(is.na(prob_idx))) {
    miss <- target_probs[is.na(prob_idx)]
    stop(sprintf("Cache %s missing probs: %s",
                 cache_base$file, paste(miss, collapse = ", ")),
         call. = FALSE)
}
if (any(is.na(prob_idx_mrts))) {
    miss <- target_probs[is.na(prob_idx_mrts)]
    stop(sprintf("Cache %s missing probs: %s",
                 cache_mrts$file, paste(miss, collapse = ", ")),
         call. = FALSE)
}

# ---- methods 2..5 (and method 1 as Gaussian reference) --------------
method_ids <- 2:5
method_label_map <- c(
    "2" = "Skew-t, K=1",
    "3" = "t, K=1, q(0.80)",
    "4" = "Skew-t, K=5",
    "5" = "t, K=5, q(0.80)"
)

verify_methods <- function(cache, needed) {
    miss <- setdiff(needed, cache$methods)
    if (length(miss) > 0L) {
        stop(sprintf("Cache %s missing methods: %s",
                     cache$file, paste(miss, collapse = ", ")), call. = FALSE)
    }
}
verify_methods(cache_base, c(1L, method_ids))
verify_methods(cache_mrts, method_ids)

m_idx_base <- match(c(1L, method_ids), cache_base$methods)
m_idx_mrts <- match(method_ids,         cache_mrts$methods)

# ---- mean over datasets, then ratio vs Gaussian baseline ------------
# cache$brier.score dim = [probs, datasets, methods]
mean_bs_base <- apply(cache_base$brier.score[prob_idx, , m_idx_base, drop = FALSE],
                      c(1, 3), mean, na.rm = TRUE)
mean_bs_mrts <- apply(cache_mrts$brier.score[prob_idx_mrts, , m_idx_mrts, drop = FALSE],
                      c(1, 3), mean, na.rm = TRUE)
# mean_bs_base columns: 1=method 1, 2..5 = method_ids
mean_ref <- mean_bs_base[, 1]
mean_baseline_models <- mean_bs_base[, -1, drop = FALSE]   # methods 2..5 (no MRTS)
mean_mrts_models     <- mean_bs_mrts                       # methods 2..5 (+ MRTS)

relative_baseline <- sweep(mean_baseline_models, 1, mean_ref, FUN = "/")
relative_mrts     <- sweep(mean_mrts_models,     1, mean_ref, FUN = "/")
relative_baseline[!is.finite(relative_baseline)] <- NA_real_
relative_mrts[!is.finite(relative_mrts)]         <- NA_real_

log_baseline <- suppressWarnings(log(relative_baseline))
log_mrts     <- suppressWarnings(log(relative_mrts))
log_baseline[!is.finite(log_baseline)] <- NA_real_
log_mrts[!is.finite(log_mrts)]         <- NA_real_

nonpos_count <- sum(
    (is.finite(relative_baseline) & relative_baseline <= 0) |
    (is.finite(relative_mrts) & relative_mrts <= 0),
    na.rm = TRUE
)

# ---- assemble long-format plot data ---------------------------------
panel_mrts_label <- sprintf("method 2-5 + mrts (K=%d)", mrts_k_active)
make_panel_rows <- function(rel, log_rel, panel_name, mrts_k) {
    out <- list()
    for (j in seq_along(method_ids)) {
        m <- method_ids[j]
        method_label <- unname(method_label_map[as.character(m)])
        label <- if (is.na(mrts_k)) method_label else
            sprintf("%s + mrts(k=%d)", method_label, mrts_k)
        out[[j]] <- data.frame(
            quantile           = target_probs,
            relative_brier     = rel[, j],
            log_relative_brier = log_rel[, j],
            panel              = panel_name,
            method_id          = m,
            mrts_k             = mrts_k,
            line_label         = label,
            stringsAsFactors   = FALSE
        )
    }
    do.call(rbind, out)
}
plot_df <- rbind(
    make_panel_rows(relative_baseline, log_baseline,
                    "method 2-5", NA_integer_),
    make_panel_rows(relative_mrts, log_mrts,
                    panel_mrts_label, mrts_k_active)
)
plot_df$panel <- factor(plot_df$panel,
                        levels = c("method 2-5", panel_mrts_label))

# ---- outputs --------------------------------------------------------
output_dir <- file.path("comparison_mrts", "plot")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

data_out <- file.path(output_dir,
    sprintf("setting%d_relative_brier_plot_data-K%d.csv",
            setting_id, mrts_k_active))
write.csv(plot_df, data_out, row.names = FALSE)

plot_file <- file.path(output_dir,
    sprintf("setting%d_relative_brier_method2to5_vs_mrts-K%d.png",
            setting_id, mrts_k_active))

# Shared palette with code/analysis/simstudy/plots.R: indexed by method id 1..8.
mlty <- c(1, 1, 3, 3, 5, 6, 2, 4)
mpch <- c(21, 22, 23, 24, 25, 4, 8, 9)
mcol <- c("gray30", "firebrick4", "dodgerblue4", "firebrick1",
          "dodgerblue1", "darkgreen", "purple4", "darkorange2")
mbg  <- c("gray70", "firebrick2", "dodgerblue2", "firebrick1",
          "dodgerblue1", "lightgreen", "plum", "moccasin")

png(filename = plot_file, width = 2600, height = 1300, res = 220)
par(mfrow = c(1, 2), mar = c(5.2, 6.2, 3.8, 1.5), oma = c(0, 0, 2, 0))

y_all <- plot_df$log_relative_brier[is.finite(plot_df$log_relative_brier)]
if (length(y_all) == 0) y_all <- c(0)
y_lim <- range(c(y_all, 0), na.rm = TRUE)
y_pad <- 0.05 * diff(y_lim)
if (!is.finite(y_pad) || y_pad == 0) y_pad <- 0.05
y_lim <- c(y_lim[1] - y_pad, y_lim[2] + y_pad)

axis_at <- target_probs
axis_show <- c(seq(0.90, 0.99, by = 0.01), 0.995)
axis_labels <- ifelse(axis_at %in% axis_show,
                      format(axis_at, trim = TRUE, scientific = FALSE), "")

for (panel_name in levels(plot_df$panel)) {
    panel_df <- plot_df[plot_df$panel == panel_name, , drop = FALSE]
    panel_df <- panel_df[order(panel_df$method_id, panel_df$quantile), ,
                         drop = FALSE]

    all_methods <- sort(unique(panel_df$method_id))
    first_method <- all_methods[1]
    first_df <- panel_df[panel_df$method_id == first_method, , drop = FALSE]

    plot(
        x = first_df$quantile,
        y = first_df$log_relative_brier,
        type = "o",
        pch  = mpch[first_method],
        lty  = mlty[first_method],
        col  = mcol[first_method],
        bg   = mbg[first_method],
        lwd  = 2,
        cex  = 1.1,
        ylim = y_lim,
        xaxt = "n",
        xlab = "Threshold quantile",
        ylab = "log(Relative Brier score) (vs method 1)",
        main = panel_name
    )

    axis(1, at = axis_at, labels = axis_labels, cex.axis = 0.9, las = 2)
    abline(h = 0, lty = 2, col = "gray60")

    if (length(all_methods) > 1) {
        for (mid in all_methods[-1]) {
            method_df <- panel_df[panel_df$method_id == mid, , drop = FALSE]
            lines(x = method_df$quantile,
                  y = method_df$log_relative_brier,
                  lty = mlty[mid], col = mcol[mid], lwd = 2)
            points(x = method_df$quantile,
                   y = method_df$log_relative_brier,
                   pch = mpch[mid], col = mcol[mid], bg = mbg[mid], cex = 1.1)
        }
    }

    legend_df <- panel_df[!duplicated(panel_df$line_label),
                          c("line_label", "method_id"), drop = FALSE]
    legend_df <- legend_df[order(legend_df$method_id), , drop = FALSE]

    legend("topleft",
           legend = legend_df$line_label,
           col    = mcol[legend_df$method_id],
           pt.bg  = mbg[legend_df$method_id],
           pch    = mpch[legend_df$method_id],
           lty    = mlty[legend_df$method_id],
           lwd    = 2,
           bty    = "n",
           cex    = 0.85)
}

mtext(sprintf("Setting %d (MRTS K=%d): log(Relative Brier Score) Profiles",
              setting_id, mrts_k_active),
      outer = TRUE, cex = 1.3, font = 2)
dev.off()

cat("Plot saved:\n  - ", plot_file, "\n", sep = "")
cat("Data saved:\n  - ", data_out, "\n", sep = "")
cat(sprintf("Datasets used (baseline cache): %d\n", length(cache_base$datasets)))
cat(sprintf("Datasets used (MRTS cache):     %d\n", length(cache_mrts$datasets)))
if (nonpos_count > 0) {
    cat(sprintf("Warning: %d non-positive relative-Brier values; log set to NA.\n",
                nonpos_count))
}
