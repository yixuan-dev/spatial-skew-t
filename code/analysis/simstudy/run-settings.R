#########################################################################
# Simulation replication script
# Data setting: configurable (default: 5 => Skew-t, K = 5, lambda = 3)
# Run selected analysis methods for one or more dataset indices.
# Optionally append MRTS basis covariates to methods 1 through 5.
#########################################################################

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
    script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = FALSE)
    script_dir <- dirname(script_path)
    if (dir.exists(script_dir)) {
        setwd(script_dir)
    }
}

source("./ar2_load.R", chdir = TRUE)
source("./mrts_cov_helpers.R", chdir = TRUE)

# -----------------------------
# User controls
# -----------------------------
default_setting <- 5
iters <- 20000
burn <- 10000
update <- 1000
thin <- 1

# Parse command-line arguments
# Usage:
#   Rscript run-settings.R [--setting=<id>|--setting <id>] [datasets] [workers] [ms_threads] [methods] [mrts_k]
# Notes:
#   - --setting (if provided) must be the first argument after script name.
#   - methods remains numeric in 1..6.
#   - mrts_k is optional and, when provided, adds MRTS-augmented variants of methods 1..5.
# Examples:
#   methods="1:6"
#   methods="(2,3,6)"
#   mrts_k="15"
#   mrts_k="c(5,10,15)"

parse_dataset_spec <- function(dataset_str) {
    dataset_ids <- parse_index_expr(dataset_str, "datasets")
    if (any(dataset_ids < 1 | dataset_ids > 50)) {
        stop("datasets must be integers in 1..50", call. = FALSE)
    }
    dataset_ids
}

parse_methods_spec <- function(methods_str) {
    method_ids <- parse_index_expr(methods_str, "methods")
    if (any(method_ids < 1 | method_ids > 6)) {
        stop("methods must be integers in 1..6", call. = FALSE)
    }
    sort(unique(as.integer(method_ids)))
}

parse_mrts_k_spec <- function(mrts_str) {
    expr <- trimws(mrts_str)
    if (!nzchar(expr)) {
        return(integer(0))
    }
    parse_positive_int_expr(expr, "mrts_k")
}

parse_setting_spec <- function(setting_str) {
    setting_ids <- parse_index_expr(setting_str, "setting")
    if (length(setting_ids) != 1) {
        stop("setting must evaluate to a single integer", call. = FALSE)
    }

    setting_id <- setting_ids[[1]]
    if (setting_id < 1) {
        stop("setting must be a positive integer", call. = FALSE)
    }

    max_setting <- NA_integer_
    if (exists("y", inherits = TRUE)) {
        y_dim <- dim(y)
        if (!is.null(y_dim) && length(y_dim) >= 4 && !is.na(y_dim[4])) {
            max_setting <- as.integer(y_dim[4])
        }
    }

    if (!is.na(max_setting) && setting_id > max_setting) {
        stop(sprintf("setting must be in 1..%d", max_setting), call. = FALSE)
    }

    setting_id
}

extract_setting_arg <- function(args) {
    if (length(args) == 0) {
        return(list(args = character(0), setting = NULL))
    }

    has_setting_anywhere <- any(grepl("^--setting=", args)) || any(args == "--setting")
    setting_value <- NULL
    remaining_args <- args

    if (grepl("^--setting=", args[1])) {
        setting_value <- sub("^--setting=", "", args[1])
        remaining_args <- if (length(args) > 1) args[-1] else character(0)
    } else if (args[1] == "--setting") {
        if (length(args) < 2) {
            stop("Missing value after --setting", call. = FALSE)
        }
        setting_value <- args[2]
        remaining_args <- if (length(args) > 2) args[-c(1, 2)] else character(0)
    } else if (has_setting_anywhere) {
        stop("--setting must be the first argument after script name", call. = FALSE)
    }

    if (any(grepl("^--setting=", remaining_args)) || any(remaining_args == "--setting")) {
        stop("Specify --setting only once and place it first", call. = FALSE)
    }

    list(args = remaining_args, setting = setting_value)
}

get_method_catalog <- function() {
    get_simstudy_method_catalog(include_maxstable = TRUE)
}

build_run_plan <- function(method_ids, mrts_k_values) {
    catalog <- get_method_catalog()
    baseline <- catalog[catalog$method_id %in% method_ids, , drop = FALSE]
    baseline$mrts_k <- NA_integer_
    baseline$output_tag <- NA_character_

    plan_cols <- c(
        "method_id", "method_key", "label", "runner", "method", "skew",
        "thresh_all", "thresh_quant", "nknots", "use_exponential_fallback",
        "mrts_k", "output_tag"
    )
    mrts_rows <- baseline[0, plan_cols, drop = FALSE]
    eligible_ids <- intersect(method_ids, 1:5)
    if (length(mrts_k_values) > 0 && length(eligible_ids) > 0) {
        eligible <- catalog[catalog$method_id %in% eligible_ids, , drop = FALSE]
        combos <- expand.grid(
            method_id = eligible$method_id,
            mrts_k = mrts_k_values,
            KEEP.OUT.ATTRS = FALSE,
            stringsAsFactors = FALSE
        )
        mrts_rows <- merge(combos, eligible, by = "method_id", sort = FALSE)
        mrts_rows$output_tag <- vapply(mrts_rows$mrts_k, format_mrts_tag, FUN.VALUE = character(1))
        mrts_rows$method_key <- paste0(mrts_rows$method_id, "+", mrts_rows$output_tag)
        mrts_rows$label <- paste0(mrts_rows$label, " + MRTS(", mrts_rows$mrts_k, ")")
        mrts_rows <- mrts_rows[, plan_cols, drop = FALSE]
    }

    if (length(mrts_k_values) > 0 && length(eligible_ids) == 0) {
        cat("MRTS requested, but selected methods do not include methods 1 through 5. No MRTS-augmented tasks added.\n")
    }

    plan <- rbind(
        baseline[, plan_cols, drop = FALSE],
        mrts_rows
    )
    rownames(plan) <- NULL
    plan$plan_id <- seq_len(nrow(plan))
    plan
}

run_method_mcmc <- function(plan_id, dataset_id) {
    spec <- method_plan[method_plan$plan_id == plan_id, , drop = FALSE]
    if (nrow(spec) != 1) {
        stop(sprintf("Expected exactly one plan row for plan_id=%s", plan_id), call. = FALSE)
    }

    obs <- c(rep(TRUE, nrow(y) - ntest), rep(FALSE, ntest))
    y.d <- y[, , dataset_id, setting]
    y.o <- y.d[obs, ]
    x.o <- x[obs, , , drop = FALSE]
    s.o <- s[obs, , drop = FALSE]
    x.p <- x[!obs, , , drop = FALSE]
    s.p <- s[!obs, , drop = FALSE]

    mrts_meta <- NULL
    x.train.use <- x.o
    x.pred.use <- x.p

    if (!is.na(spec$mrts_k[1])) {
        mrts_meta <- build_mrts_covariates(
            S_train = s.o,
            S_pred = s.p,
            nt = dim(x.o)[2],
            k = spec$mrts_k[1]
        )
        extended_x <- append_mrts_covariates(
            X_train = x.o,
            X_pred = x.p,
            mrts_cov = mrts_meta
        )
        x.train.use <- extended_x$train
        x.pred.use <- extended_x$pred
    }

    set.seed(get_simstudy_seed(
        setting_id = setting,
        method_id = spec$method_id[1],
        dataset_id = dataset_id,
        mrts_k = spec$mrts_k[1]
    ))

    outputfile <- build_simstudy_result_file(
        results_dir = results_dir,
        setting_id = setting,
        method_id = spec$method_id[1],
        dataset_id = dataset_id,
        mrts_k = spec$mrts_k[1]
    )

    cat(sprintf("[Dataset %d][Method %s] start\n", dataset_id, spec$method_key[1]))
    if (!is.null(mrts_meta)) {
        cat(sprintf(
            "[Dataset %d][Method %s] MRTS source=%s cols=%d/%d/%d\n",
            dataset_id,
            spec$method_key[1],
            mrts_meta$source,
            mrts_meta$original_cols,
            mrts_meta$kept_cols,
            mrts_meta$dropped_cols
        ))
    }

    fit.1 <- NULL
    if (isTRUE(spec$use_exponential_fallback[1])) {
        fit.1 <- tryCatch(
            mcmc(
                y = y.o,
                x = x.train.use,
                s = s.o,
                s.pred = s.p,
                x.pred = x.pred.use,
                method = spec$method[1],
                skew = isTRUE(spec$skew[1]),
                thresh.all = spec$thresh_all[1],
                thresh.quant = isTRUE(spec$thresh_quant[1]),
                nknots = spec$nknots[1],
                iterplot = FALSE,
                iters = iters,
                burn = burn,
                update = update,
                min.s = c(0, 0),
                max.s = c(10, 10),
                temporalw = FALSE,
                temporaltau = FALSE,
                temporalz = FALSE,
                rho.upper = 15,
                nu.upper = 10
            ),
            error = function(e) {
                mcmc(
                    y = y.o,
                    x = x.train.use,
                    s = s.o,
                    s.pred = s.p,
                    x.pred = x.pred.use,
                    method = spec$method[1],
                    skew = isTRUE(spec$skew[1]),
                    thresh.all = spec$thresh_all[1],
                    thresh.quant = isTRUE(spec$thresh_quant[1]),
                    nknots = spec$nknots[1],
                    iterplot = FALSE,
                    iters = iters,
                    burn = burn,
                    update = update,
                    min.s = c(0, 0),
                    max.s = c(10, 10),
                    temporalw = FALSE,
                    temporaltau = FALSE,
                    temporalz = FALSE,
                    rho.upper = 15,
                    nu.upper = 10,
                    cov.model = "exponential"
                )
            }
        )
    } else {
        fit.1 <- mcmc(
            y = y.o,
            x = x.train.use,
            s = s.o,
            s.pred = s.p,
            x.pred = x.pred.use,
            method = spec$method[1],
            skew = isTRUE(spec$skew[1]),
            thresh.all = spec$thresh_all[1],
            thresh.quant = isTRUE(spec$thresh_quant[1]),
            nknots = spec$nknots[1],
            iterplot = FALSE,
            iters = iters,
            burn = burn,
            update = update,
            min.s = c(0, 0),
            max.s = c(10, 10),
            temporalw = FALSE,
            temporaltau = FALSE,
            temporalz = FALSE,
            rho.upper = 15,
            nu.upper = 10
        )
    }

    if (is.na(spec$mrts_k[1])) {
        save(fit.1, file = outputfile)
    } else {
        analysis_spec <- spec
        save(fit.1, analysis_spec, mrts_meta, file = outputfile)
    }

    rm(fit.1)
    if (exists("analysis_spec")) {
        rm(analysis_spec)
    }
    if (exists("mrts_meta")) {
        rm(mrts_meta)
    }
    gc()
    cat(sprintf("[Dataset %d][Method %s] done -> %s\n", dataset_id, spec$method_key[1], outputfile))
}

run_method_maxstable <- function(dataset_id) {
    spec <- method_plan[method_plan$method_id == 6 & is.na(method_plan$mrts_k), , drop = FALSE]
    if (nrow(spec) != 1) {
        stop("Expected exactly one baseline max-stable specification.", call. = FALSE)
    }

    obs <- c(rep(TRUE, nrow(y) - ntest), rep(FALSE, ntest))
    y.d <- y[, , dataset_id, setting]
    y.o <- y.d[obs, ]
    x.o <- x[obs, , , drop = FALSE]
    s.o <- s[obs, , drop = FALSE]
    x.p <- x[!obs, , , drop = FALSE]
    s.p <- s[!obs, , drop = FALSE]

    set.seed(get_simstudy_seed(setting, 6L, dataset_id))
    outputfile <- build_simstudy_result_file(
        results_dir = results_dir,
        setting_id = setting,
        method_id = 6L,
        dataset_id = dataset_id
    )

    knots.x <- seq(1, 9, length = 12)
    knots <- expand.grid(knots.x, knots.x)
    y.ms <- t(y.o)
    thresh <- quantile(y.ms, probs = 0.80, na.rm = TRUE)

    cat(sprintf("[Dataset %d][Method 6] Max-stable, T=q(0.80) start\n", dataset_id))
    fit.1 <- maxstable(
        y = y.ms,
        x = x.o,
        s = s.o,
        sp = s.p,
        xp = x.p,
        thresh = thresh,
        knots = knots,
        iters = iters,
        burn = burn,
        update = update,
        threads = ms_threads,
        thin = thin
    )

    save(fit.1, file = outputfile)
    rm(fit.1)
    gc()
    cat(sprintf("[Dataset %d][Method 6] done -> %s\n", dataset_id, outputfile))
}

args_raw <- commandArgs(trailingOnly = TRUE)
parsed_args <- extract_setting_arg(args_raw)
args <- parsed_args$args
setting_flag_value <- parsed_args$setting

datasets_spec <- if (length(args) >= 1) args[1] else "1"
workers <- if (length(args) >= 2) as.integer(args[2]) else 1
ms_threads <- if (length(args) >= 3) as.integer(args[3]) else 2
methods_spec <- if (length(args) >= 4) args[4] else "1:6"
mrts_k_spec <- if (length(args) >= 5) args[5] else ""

setting <- if (!is.null(setting_flag_value)) {
    parse_setting_spec(setting_flag_value)
} else {
    parse_setting_spec(as.character(default_setting))
}

if (is.na(workers) || workers < 1) {
    stop("workers must be a positive integer", call. = FALSE)
}
if (is.na(ms_threads) || ms_threads < 1) {
    stop("ms_threads must be a positive integer", call. = FALSE)
}

dataset_ids <- parse_dataset_spec(datasets_spec)
method_ids <- parse_methods_spec(methods_spec)
mrts_k_values <- parse_mrts_k_spec(mrts_k_spec)

results_dir <- "results"
if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
}

method_plan <- build_run_plan(method_ids = method_ids, mrts_k_values = mrts_k_values)
write.csv(method_plan, file.path(results_dir, sprintf("run-plan-setting-%d.csv", setting)), row.names = FALSE)

cat(sprintf(
    "Run setting=%d, datasets=%s, iters=%d, burn=%d, workers=%d, ms_threads=%d\n",
    setting, paste(dataset_ids, collapse = ","), iters, burn, workers, ms_threads
))
cat(sprintf("Methods to run: %s\n", paste(method_plan$method_key, collapse = ", ")))
if (length(mrts_k_values) > 0) {
    cat(sprintf("Requested MRTS k values: %s\n", paste(mrts_k_values, collapse = ", ")))
}

mcmc_plan <- method_plan[method_plan$runner == "mcmc", , drop = FALSE]
run_maxstable <- any(method_plan$runner == "maxstable")
mcmc_tasks <- if (nrow(mcmc_plan) > 0) {
    expand.grid(
        dataset = dataset_ids,
        plan_id = mcmc_plan$plan_id,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
} else {
    data.frame(dataset = integer(0), plan_id = integer(0))
}

workers_use <- min(workers, nrow(mcmc_tasks))

if (workers_use > 1) {
    cat(sprintf(
        "Running MCMC methods in parallel with %d workers across %d tasks...\n",
        workers_use,
        nrow(mcmc_tasks)
    ))
    cl <- parallel::makeCluster(workers_use, type = "PSOCK")

    wd <- getwd()
    parallel::clusterExport(cl, varlist = c("wd"), envir = environment())
    parallel::clusterEvalQ(cl, {
        setwd(wd)
        source("./ar2_load.R", chdir = TRUE)
        source("./mrts_cov_helpers.R", chdir = TRUE)
        NULL
    })
    parallel::clusterExport(
        cl,
        varlist = c(
            "setting", "iters", "burn", "update", "thin", "results_dir",
            "method_plan", "mcmc_tasks", "run_method_mcmc"
        ),
        envir = environment()
    )

    elapsed <- parallel::parLapply(cl, seq_len(nrow(mcmc_tasks)), function(i) {
        d <- mcmc_tasks$dataset[i]
        pid <- mcmc_tasks$plan_id[i]
        tic <- proc.time()
        run_method_mcmc(pid, d)
        unname((proc.time() - tic)[3])
    })
    for (i in seq_len(nrow(mcmc_tasks))) {
        method_key <- method_plan$method_key[match(mcmc_tasks$plan_id[i], method_plan$plan_id)]
        cat(sprintf(
            "elapsed (dataset %d, method %s): %.2f sec\n",
            mcmc_tasks$dataset[i],
            method_key,
            elapsed[[i]]
        ))
    }

    parallel::stopCluster(cl)
} else {
    for (i in seq_len(nrow(mcmc_tasks))) {
        d <- mcmc_tasks$dataset[i]
        pid <- mcmc_tasks$plan_id[i]
        method_key <- method_plan$method_key[match(pid, method_plan$plan_id)]
        tic <- proc.time()
        run_method_mcmc(pid, d)
        toc <- proc.time()
        cat(sprintf("elapsed (dataset %d, method %s): %.2f sec\n", d, method_key, (toc - tic)[3]))
    }
}

if (run_maxstable) {
    for (d in dataset_ids) {
        tic <- proc.time()
        run_method_maxstable(d)
        toc <- proc.time()
        cat(sprintf("elapsed (dataset %d, method 6): %.2f sec\n", d, (toc - tic)[3]))
    }
}

cat("All specified methods and datasets finished.\n")
