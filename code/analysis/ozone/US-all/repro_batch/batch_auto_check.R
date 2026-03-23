rm(list = ls())

# Repro batch auto runner + health check
# Usage examples:
#   Rscript repro_batch/batch_auto_check.R
#   Rscript repro_batch/batch_auto_check.R --settings=33,34,35,36,38,39,40,41,43,44,45,46
#   Rscript repro_batch/batch_auto_check.R --settings=51,52,53,54,55,56 --label=F3
#   Rscript repro_batch/batch_auto_check.R --batch=F2
#   Rscript repro_batch/batch_auto_check.R --batch=AR2_SMOKE --label=AR2-SMOKE
#   Rscript repro_batch/batch_auto_check.R --dry-run # just check which results are present without running any scripts
#   Rscript repro_batch/batch_auto_check.R --force
#   Rscript repro_batch/batch_auto_check.R --full-done
#   Rscript repro_batch/batch_auto_check.R --full-done --dry-run
#   Rscript repro_batch/batch_auto_check.R --full-done --dry-run --stop-on-batch-failure
#   Rscript repro_batch/batch_auto_check.R --full-done --dry-run --continue-on-error
# Notes:
# - The --full-done mode runs all batches sequentially (F1 to F6) and creates checkpoints after each batch. By default, it continues even if a batch has failures, but you can change this behavior with --stop-on-batch-failure or --continue-on-error.
#
# & "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" ".\repro_batch\batch_auto_check.R" --batch=F6 --label=F6 --stop-on-batch-failure

main <- function() {
    default_f1_settings <- c(1, 2, 3, 4, 5, 7, 8, 9, 11, 12, 13, 15, 16, 17)
    ar2_smoke_settings <- c(101, 102, 103)

    batch_map <- list(
        F1 = default_f1_settings,
        F2 = c(33, 34, 35, 36, 38, 39, 40, 41, 43, 44, 45, 46),
        F3 = c(51, 52, 53, 54, 55, 56),
        F4 = c(57, 58, 59, 60, 61, 62),
        F5 = c(63, 64, 65, 66, 67, 68),
        F6 = c(69, 70, 71, 72, 73, 74),
        F7 = 101:108,
        F8 = 109:116,
        F9 = 117:124,
        AR2_SMOKE = ar2_smoke_settings,
        A1 = ar2_smoke_settings
    )

    args <- commandArgs(trailingOnly = TRUE)
    dry_run <- "--dry-run" %in% args
    force_run <- "--force" %in% args
    full_done <- "--full-done" %in% args
    stop_on_batch_failure <- "--stop-on-batch-failure" %in% args
    continue_on_error_flag <- "--continue-on-error" %in% args

    if (stop_on_batch_failure && continue_on_error_flag) {
        stop("Please choose one policy: --stop-on-batch-failure OR --continue-on-error, not both.")
    }

    # default behavior keeps historical behavior: continue even if a batch has failures
    continue_on_error <- if (stop_on_batch_failure) FALSE else TRUE

    get_arg_value <- function(args, key) {
        prefix <- paste0("--", key, "=")
        hit <- grep(paste0("^", prefix), args, value = TRUE)
        if (length(hit) > 0) {
            return(sub(prefix, "", hit[1]))
        }
        idx <- which(args == paste0("--", key))
        if (length(idx) > 0 && idx[1] < length(args)) {
            return(args[idx[1] + 1])
        }
        return(NA_character_)
    }

    parse_settings <- function(x) {
        x <- gsub("\\s+", "", x)
        if (!nzchar(x)) {
            stop("--settings provided but empty.")
        }
        vals <- suppressWarnings(as.integer(strsplit(x, ",", fixed = TRUE)[[1]]))
        if (any(is.na(vals))) {
            stop("--settings must be a comma-separated integer list, e.g. --settings=33,34,35")
        }
        vals <- unique(vals)
        vals <- vals[order(vals)]
        vals
    }

    settings_arg <- get_arg_value(args, "settings")
    batch_arg <- toupper(get_arg_value(args, "batch"))
    label_arg <- get_arg_value(args, "label")

    if (full_done && (!is.na(settings_arg) || !is.na(batch_arg))) {
        stop("--full-done cannot be combined with --settings or --batch.")
    }

    if (!is.na(settings_arg) && !is.na(batch_arg)) {
        stop("Please provide either --settings or --batch, not both.")
    }

    if (full_done) {
        run_settings <- NULL
        run_label <- ifelse(is.na(label_arg), "full-done", label_arg)
    } else if (!is.na(settings_arg)) {
        run_settings <- parse_settings(settings_arg)
        run_label <- ifelse(is.na(label_arg), "custom", label_arg)
    } else if (!is.na(batch_arg)) {
        if (!batch_arg %in% names(batch_map)) {
            stop("Unknown --batch value. Use one of: ", paste(names(batch_map), collapse = ", "))
        }
        run_settings <- batch_map[[batch_arg]]
        run_label <- ifelse(is.na(label_arg), batch_arg, label_arg)
    } else {
        run_settings <- default_f1_settings
        run_label <- ifelse(is.na(label_arg), "F1", label_arg)
    }

    script_path_arg <- grep("^--file=", commandArgs(), value = TRUE)
    if (length(script_path_arg) > 0) {
        script_path <- sub("^--file=", "", script_path_arg[1])
        script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
    } else {
        script_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
    }

    project_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
    results_dir <- file.path(project_dir, "results")
    output_dir <- file.path(script_dir, "output")

    if (!dir.exists(results_dir)) {
        dir.create(results_dir, recursive = TRUE)
    }
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }

    cat("Run label:", run_label, "\n")
    if (full_done) {
        cat("Mode: full-done (F1 -> F6 sequential)\n")
        cat("Batch failure policy:", ifelse(continue_on_error, "continue-on-error", "stop-on-batch-failure"), "\n")
    } else {
        cat("Settings:", paste(run_settings, collapse = ", "), "\n")
    }
    cat("dry_run =", dry_run, "; force_run =", force_run, "\n\n")

    run_one_setting <- function(setting, project_dir, results_dir, dry_run = FALSE, force_run = FALSE) {
        script_file <- file.path(project_dir, sprintf("us-all-%d.R", setting))
        result_file <- file.path(results_dir, sprintf("us-all-%d.RData", setting))

        start_time <- Sys.time()
        elapsed_sec <- NA_real_
        run_status <- NA_character_
        error_message <- ""

        if (!file.exists(script_file)) {
            run_status <- "script_missing"
        } else if (dry_run) {
            run_status <- ifelse(file.exists(result_file), "dry_present", "dry_missing")
        } else if (file.exists(result_file) && !force_run) {
            run_status <- "already_present"
        } else {
            err <- NULL
            tic <- proc.time()[3]

            tryCatch(
                {
                    script_env <- new.env(parent = globalenv())
                    source(script_file, chdir = TRUE, local = script_env)
                },
                error = function(e) {
                    err <<- conditionMessage(e)
                }
            )

            elapsed_sec <- proc.time()[3] - tic

            if (is.null(err)) {
                if (file.exists(result_file)) {
                    run_status <- "ok"
                } else {
                    run_status <- "missing_result"
                    error_message <- "Script finished but result file not found."
                }
            } else {
                error_message <- err
                if (file.exists(result_file)) {
                    run_status <- "error_result_exists"
                } else {
                    run_status <- "error"
                }
            }
        }

        result_exists <- file.exists(result_file)

        cat(sprintf(
            "[%s] setting=%d | status=%s | result=%s\n",
            format(Sys.time(), "%H:%M:%S"),
            setting,
            run_status,
            ifelse(result_exists, "FOUND", "MISSING")
        ))

        if (!identical(error_message, "")) {
            cat("  error:", error_message, "\n")
        }

        data.frame(
            setting = setting,
            script_file = normalizePath(script_file, winslash = "/", mustWork = FALSE),
            result_file = normalizePath(result_file, winslash = "/", mustWork = FALSE),
            run_status = run_status,
            result_exists = result_exists,
            elapsed_sec = elapsed_sec,
            checked_at = format(start_time, "%Y-%m-%d %H:%M:%S"),
            error_message = error_message,
            stringsAsFactors = FALSE
        )
    }

    failure_statuses <- c(
        "script_missing",
        "missing_result",
        "error",
        "error_result_exists",
        "dry_missing"
    )

    write_checkpoint <- function(status_df, run_label, output_dir, failure_statuses, summary_title = "Batch Auto-check Summary") {
        failed_settings <- status_df$setting[status_df$run_status %in% failure_statuses]

        stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
        status_csv_time <- file.path(output_dir, sprintf("%s-run-status-%s.csv", run_label, stamp))
        status_csv_latest <- file.path(output_dir, sprintf("%s-run-status-latest.csv", run_label))
        fail_txt_time <- file.path(output_dir, sprintf("%s-failures-%s.txt", run_label, stamp))
        fail_txt_latest <- file.path(output_dir, sprintf("%s-failures-latest.txt", run_label))

        write.csv(status_df, status_csv_time, row.names = FALSE)
        write.csv(status_df, status_csv_latest, row.names = FALSE)

        if (length(failed_settings) > 0) {
            fail_lines <- c(
                sprintf("Generated at: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                sprintf("Run label: %s", run_label),
                sprintf("Failed settings (%d):", length(failed_settings)),
                paste(failed_settings, collapse = ", "),
                "",
                "Rerun helpers:",
                paste(sprintf("source(\"us-all-%d.R\")", failed_settings), collapse = "\n")
            )
        } else {
            fail_lines <- c(
                sprintf("Generated at: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                sprintf("Run label: %s", run_label),
                "No failed settings in this run."
            )
        }

        writeLines(fail_lines, con = fail_txt_time)
        writeLines(fail_lines, con = fail_txt_latest)

        cat("\n=====", summary_title, "=====\n")
        print(table(status_df$run_status, useNA = "ifany"))

        cat("\nSaved:\n")
        cat("-", normalizePath(status_csv_time, winslash = "/", mustWork = FALSE), "\n")
        cat("-", normalizePath(status_csv_latest, winslash = "/", mustWork = FALSE), "\n")
        cat("-", normalizePath(fail_txt_time, winslash = "/", mustWork = FALSE), "\n")
        cat("-", normalizePath(fail_txt_latest, winslash = "/", mustWork = FALSE), "\n")

        if (length(failed_settings) > 0) {
            cat("\nFailed settings:", paste(failed_settings, collapse = ", "), "\n")
        } else {
            cat("\nAll settings passed file-existence checks.\n")
        }

        list(
            failed_settings = failed_settings,
            status_csv_time = status_csv_time,
            status_csv_latest = status_csv_latest,
            fail_txt_time = fail_txt_time,
            fail_txt_latest = fail_txt_latest
        )
    }

    run_settings_block <- function(settings_vec, block_label, project_dir, results_dir, dry_run = FALSE, force_run = FALSE) {
        cat("\n--- Running", block_label, "---\n")
        cat("Settings:", paste(settings_vec, collapse = ", "), "\n")

        status_list <- lapply(
            settings_vec,
            run_one_setting,
            project_dir = project_dir,
            results_dir = results_dir,
            dry_run = dry_run,
            force_run = force_run
        )

        status_df <- do.call(rbind, status_list)
        status_df$batch <- block_label
        status_df
    }

    if (full_done) {
        full_batch_order <- c("F1", "F2", "F3", "F4", "F5", "F6")
        all_status <- list()
        completed_batches <- character(0)
        stop_triggered <- FALSE

        for (b in full_batch_order) {
            batch_status <- run_settings_block(
                settings_vec = batch_map[[b]],
                block_label = b,
                project_dir = project_dir,
                results_dir = results_dir,
                dry_run = dry_run,
                force_run = force_run
            )

            all_status[[b]] <- batch_status

            # per-batch checkpoint
            checkpoint_info <- write_checkpoint(
                status_df = batch_status,
                run_label = b,
                output_dir = output_dir,
                failure_statuses = failure_statuses,
                summary_title = paste("Checkpoint", b)
            )
            invisible(checkpoint_info)

            completed_batches <- c(completed_batches, b)
            if (!continue_on_error && length(checkpoint_info$failed_settings) > 0) {
                cat("\nStopping early because", b, "has failures and --stop-on-batch-failure is active.\n")
                stop_triggered <- TRUE
                break
            }
        }

        full_status_df <- do.call(rbind, all_status)
        full_done_info <- write_checkpoint(
            status_df = full_status_df,
            run_label = run_label,
            output_dir = output_dir,
            failure_statuses = failure_statuses,
            summary_title = "Full-done Auto-check Summary"
        )
        invisible(full_done_info)

        if (stop_triggered) {
            remaining_batches <- setdiff(full_batch_order, completed_batches)
            cat("Skipped batches:", paste(remaining_batches, collapse = ", "), "\n")
        }
    } else {
        status_df <- run_settings_block(
            settings_vec = run_settings,
            block_label = run_label,
            project_dir = project_dir,
            results_dir = results_dir,
            dry_run = dry_run,
            force_run = force_run
        )

        checkpoint_info <- write_checkpoint(
            status_df = status_df,
            run_label = run_label,
            output_dir = output_dir,
            failure_statuses = failure_statuses,
            summary_title = "Batch Auto-check Summary"
        )
        invisible(checkpoint_info)
    }
}

main()
