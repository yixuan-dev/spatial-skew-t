# E1: MRTS constant-column A/B under Morris prop MCMC
#
# Variant A (trim):  prop_trim_constant = TRUE  (current default)
# Variant B (keep):  prop_trim_constant = FALSE (retain mrts constant column)
#
# Phase A: simdata setting 1, Gaussian
# Phase B: simdata setting 4, Skew-t K=1
#
# Constant injection: Y <- Y + a on all observed cells; true beta0 = 10 + a
#
# Env overrides: E1_ITERS, E1_BURN, E1_PROP_K, E1_INJECTION (comma-separated),
#               E1_DATASETS (e.g. "1:10" or "1,3,5"; default "1:10")

suppressMessages({
  library(fields)
  library(emulator)
})

prop_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
simdata_candidates <- c(
  file.path(prop_dir, "../../analysis/simstudy/simdata.RData"),
  file.path(prop_dir, "../analysis/simstudy/simdata.RData")
)
simdata_path <- simdata_candidates[file.exists(simdata_candidates)][1]
if (is.na(simdata_path)) {
  stop("simdata.RData not found.", call. = FALSE)
}

parse_numeric_vec <- function(x, default) {
  if (is.null(x) || !nzchar(x)) {
    return(default)
  }
  as.numeric(strsplit(x, ",", fixed = TRUE)[[1]])
}

parse_dataset_ids <- function(x, nsets) {
  x <- trimws(x)
  if (grepl(":", x, fixed = TRUE)) {
    parts <- as.integer(strsplit(x, ":", fixed = TRUE)[[1]])
    if (length(parts) != 2L) {
      stop("E1_DATASETS range must look like start:end", call. = FALSE)
    }
    return(seq.int(parts[1], parts[2]))
  }
  as.integer(strsplit(x, ",", fixed = TRUE)[[1]])
}

load(simdata_path)
dataset_ids <- parse_dataset_ids(Sys.getenv("E1_DATASETS", unset = "1:10"), nsets = nsets)
if (any(dataset_ids < 1L | dataset_ids > nsets)) {
  stop(sprintf("E1_DATASETS out of range 1:%d", nsets), call. = FALSE)
}
source(file.path(prop_dir, "load_ar2.R"))
source(file.path(prop_dir, "prop_utils.R"))
source(file.path(prop_dir, "prop_basis.R"))
source(file.path(prop_dir, "prop_covariance.R"))
source(file.path(prop_dir, "prop_imputation.R"))
source(file.path(prop_dir, "prop_modules.R"))
source(file.path(prop_dir, "mcmc_prop.R"))

iters <- as.integer(Sys.getenv("E1_ITERS", unset = "8000"))
burn <- as.integer(Sys.getenv("E1_BURN", unset = "3000"))
prop_k <- as.integer(Sys.getenv("E1_PROP_K", unset = "20"))
injections <- parse_numeric_vec(Sys.getenv("E1_INJECTION", unset = ""), c(0, 2, 5))
beta_true_base <- c(10, 0, 0)

e1_seed <- function(setting_id, method_id, dataset_id, prop_k, injection_a) {
  as.integer(method_id) * 100000L +
    as.integer(prop_k) * 1000L +
    as.integer(setting_id) * 100L +
    as.integer(dataset_id) +
    as.integer(round(injection_a * 10))
}

ess_from_chain <- function(x) {
  x <- as.numeric(x)
  n <- length(x)
  if (n < 3L) {
    return(NA_real_)
  }
  acf_obj <- stats::acf(x, plot = FALSE, lag.max = min(n - 1L, 100L))
  acf_vals <- as.numeric(acf_obj$acf[-1, 1, 1])
  pos_acf <- acf_vals[acf_vals > 0]
  if (length(pos_acf) == 0L) {
    return(n)
  }
  n / (1 + 2 * sum(pos_acf))
}

prepare_obs_data <- function(setting_id, dataset_id) {
  obs <- c(rep(TRUE, nrow(y) - ntest), rep(FALSE, ntest))
  y.d <- y[, , dataset_id, setting_id]
  list(
    y = y.d[obs, , drop = FALSE],
    x = x[obs, , , drop = FALSE],
    s = s[obs, , drop = FALSE]
  )
}

run_e1_cell <- function(phase, setting_id, method_id, spec, injection_a, trim_constant, dataset_id) {
  dat <- prepare_obs_data(setting_id, dataset_id)
  y_inj <- dat$y + injection_a
  seed <- e1_seed(setting_id, method_id, dataset_id, prop_k, injection_a)
  set.seed(seed)

  tic <- proc.time()
  fit <- mcmc(
    y = y_inj,
    x = dat$x,
    s = dat$s,
    s.pred = NULL,
    x.pred = NULL,
    method = spec$method,
    skew = isTRUE(spec$skew),
    thresh.all = spec$thresh_all,
    thresh.quant = isTRUE(spec$thresh_quant),
    nknots = spec$nknots,
    iterplot = FALSE,
    iters = iters,
    burn = burn,
    update = max(500L, iters %/% 10L),
    thin = 1L,
    min.s = c(0, 0),
    max.s = c(10, 10),
    temporalw = FALSE,
    temporaltau = FALSE,
    temporalz = FALSE,
    prop_k = prop_k,
    prop_cov_update_every = 1L,
    prop_trim_constant = isTRUE(trim_constant)
  )
  wall_sec <- unname((proc.time() - tic)[3])

  beta0_chain <- fit$beta[, 1]
  beta0_mean <- mean(beta0_chain)
  beta0_sd <- stats::sd(beta0_chain)
  beta0_true <- beta_true_base[1] + injection_a

  data.frame(
    replicate = dataset_id,
    phase = phase,
    setting = setting_id,
    method_id = method_id,
    variant = if (trim_constant) "A_trim" else "B_keep",
    trim_constant = isTRUE(trim_constant),
    injection_a = injection_a,
    beta0_true = beta0_true,
    beta0_mean = beta0_mean,
    beta0_sd = beta0_sd,
    beta0_bias = beta0_mean - beta0_true,
    beta0_ess = ess_from_chain(beta0_chain),
    sigma_xi_mean = mean(fit$prop$sigma_xi),
    rank_kept = fit$prop$basis_rank_kept[1],
    kept_cols = fit$prop$basis_kept_cols[1],
    wall_sec = wall_sec,
    seed = seed,
    stringsAsFactors = FALSE
  )
}

judge_cell <- function(row, phase) {
  bias_tol <- if (phase == "A_gaussian") 0.5 else 1.0
  bias_ok <- abs(row$beta0_bias) < bias_tol
  if (bias_ok) "PASS" else "WARN"
}

summarize_injection <- function(df_sub) {
  if (nrow(df_sub) < 2L) {
    return(list(slope = NA_real_, intercept = NA_real_, injection_ok = NA))
  }
  fit_lm <- stats::lm(beta0_mean ~ injection_a, data = df_sub)
  coefs <- stats::coef(fit_lm)
  slope <- unname(coefs["injection_a"])
  intercept <- unname(coefs["(Intercept)"])
  injection_ok <- is.finite(slope) && abs(slope - 1) < 0.25 && abs(intercept - beta_true_base[1]) < 1.0
  list(slope = slope, intercept = intercept, injection_ok = injection_ok)
}

phase_specs <- list(
  A_gaussian = list(
    setting_id = 1L,
    method_id = 1L,
    spec = list(
      method = "gaussian",
      skew = FALSE,
      thresh_all = 0,
      thresh_quant = TRUE,
      nknots = 1L
    )
  ),
  B_skewt_k1 = list(
    setting_id = 4L,
    method_id = 2L,
    spec = list(
      method = "t",
      skew = TRUE,
      thresh_all = 0,
      thresh_quant = TRUE,
      nknots = 1L
    )
  )
)

rows <- list()
for (dataset_id in dataset_ids) {
  cat(sprintf("\n[E1] ===== replicate / dataset %d =====\n", dataset_id))
  for (phase_name in names(phase_specs)) {
    cfg <- phase_specs[[phase_name]]
    for (trim_constant in c(TRUE, FALSE)) {
      for (injection_a in injections) {
        cat(sprintf(
          "[E1] rep=%d %s variant=%s a=%.1f ...\n",
          dataset_id,
          phase_name,
          if (trim_constant) "A_trim" else "B_keep",
          injection_a
        ))
        rows[[length(rows) + 1L]] <- run_e1_cell(
          phase = phase_name,
          setting_id = cfg$setting_id,
          method_id = cfg$method_id,
          spec = cfg$spec,
          injection_a = injection_a,
          trim_constant = trim_constant,
          dataset_id = dataset_id
        )
      }
    }
  }
}

results <- do.call(rbind, rows)
results$bias_flag <- mapply(judge_cell, split(results, seq_len(nrow(results))), results$phase,
  SIMPLIFY = TRUE, USE.NAMES = FALSE)

variance_flags <- rep(NA_character_, nrow(results))
for (dataset_id in dataset_ids) {
  for (phase_name in names(phase_specs)) {
    for (injection_a in injections) {
      idx <- which(
        results$replicate == dataset_id &
          results$phase == phase_name &
          results$injection_a == injection_a
      )
      sd_trim <- results$beta0_sd[results$replicate == dataset_id &
        results$phase == phase_name &
        results$injection_a == injection_a &
        results$variant == "A_trim"]
      sd_keep <- results$beta0_sd[results$replicate == dataset_id &
        results$phase == phase_name &
        results$injection_a == injection_a &
        results$variant == "B_keep"]
      if (length(sd_trim) == 1L && length(sd_keep) == 1L && sd_keep > 2 * sd_trim) {
        variance_flags[idx[results$variant[idx] == "B_keep"]] <- "WARN_variance"
      }
    }
  }
}
results$variance_flag <- variance_flags

injection_summary <- list()
for (dataset_id in dataset_ids) {
  for (phase_name in names(phase_specs)) {
    for (variant in c("A_trim", "B_keep")) {
      sub <- results[results$replicate == dataset_id &
        results$phase == phase_name & results$variant == variant, , drop = FALSE]
      inj <- summarize_injection(sub)
      injection_summary[[paste(dataset_id, phase_name, variant, sep = "|")]] <- inj
    }
  }
}

absorption_flags <- rep(NA_character_, nrow(results))
for (dataset_id in dataset_ids) {
  for (phase_name in names(phase_specs)) {
    for (variant in c("A_trim", "B_keep")) {
      sub <- results[results$replicate == dataset_id &
        results$phase == phase_name & results$variant == variant, , drop = FALSE]
      if (nrow(sub) < 2L) next
      ord <- order(sub$injection_a)
      sub <- sub[ord, , drop = FALSE]
      d_bias <- diff(sub$beta0_bias)
      d_sigma <- diff(sub$sigma_xi_mean)
      if (any(abs(d_bias) > 1) && all(abs(d_sigma) < 0.05)) {
        idx <- which(
          results$replicate == dataset_id &
            results$phase == phase_name & results$variant == variant
        )
        absorption_flags[idx] <- "WARN_absorption"
      }
    }
  }
}
results$absorption_flag <- absorption_flags

out_dir <- file.path(prop_dir, "output", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_csv <- file.path(out_dir, "e1_constant_column_ab.csv")
out_rep_csv <- file.path(out_dir, "e1_constant_column_ab_replicates.csv")
out_summary_csv <- file.path(out_dir, "e1_constant_column_ab_summary.csv")
write.csv(results, out_rep_csv, row.names = FALSE)
if (length(dataset_ids) == 1L) {
  write.csv(results, out_csv, row.names = FALSE)
} else {
  write.csv(results, out_csv, row.names = FALSE)
}

agg_keys <- c("phase", "variant", "injection_a")
summary_rows <- lapply(split(results, interaction(results[agg_keys], drop = TRUE)), function(sub) {
  data.frame(
    phase = sub$phase[1],
    variant = sub$variant[1],
    injection_a = sub$injection_a[1],
    n_reps = nrow(sub),
    beta0_bias_mean = mean(sub$beta0_bias),
    beta0_bias_sd = stats::sd(sub$beta0_bias),
    beta0_sd_mean = mean(sub$beta0_sd),
    beta0_sd_sd = stats::sd(sub$beta0_sd),
    sigma_xi_mean_avg = mean(sub$sigma_xi_mean),
    bias_pass_rate = mean(sub$bias_flag == "PASS"),
    variance_warn_rate = mean(sub$variance_flag == "WARN_variance", na.rm = TRUE),
    absorption_warn_rate = mean(sub$absorption_flag == "WARN_absorption", na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})
summary_df <- do.call(rbind, summary_rows)
rownames(summary_df) <- NULL
write.csv(summary_df, out_summary_csv, row.names = FALSE)

cat("\n=== E1 replicate summary (aggregated) ===\n")
print(summary_df)

cat("\n=== E1 results (first replicate only) ===\n")
first_rep <- results[results$replicate == dataset_ids[1], , drop = FALSE]
print(first_rep[, c(
  "replicate", "phase", "variant", "injection_a", "beta0_mean", "beta0_sd",
  "beta0_bias", "sigma_xi_mean", "rank_kept", "bias_flag", "variance_flag"
)])

cat("\n=== Injection slopes per replicate (beta0_mean ~ a; target slope=1) ===\n")
for (nm in names(injection_summary)) {
  inj <- injection_summary[[nm]]
  cat(sprintf(
    "%s: slope=%.3f intercept=%.3f injection_ok=%s\n",
    nm, inj$slope, inj$intercept,
    if (isTRUE(inj$injection_ok)) "PASS" else if (is.na(inj$injection_ok)) "NA" else "WARN"
  ))
}

slope_df <- do.call(rbind, lapply(names(injection_summary), function(nm) {
  parts <- strsplit(nm, "|", fixed = TRUE)[[1]]
  inj <- injection_summary[[nm]]
  data.frame(
    replicate = as.integer(parts[1]),
    phase = parts[2],
    variant = parts[3],
    slope = inj$slope,
    intercept = inj$intercept,
    injection_ok = isTRUE(inj$injection_ok),
    stringsAsFactors = FALSE
  )
}))
slope_summary <- aggregate(
  cbind(slope, intercept) ~ phase + variant,
  data = slope_df,
  FUN = mean,
  na.rm = TRUE
)
ok_rates <- aggregate(
  injection_ok ~ phase + variant,
  data = slope_df,
  FUN = mean,
  na.rm = TRUE
)
names(ok_rates)[3] <- "injection_pass_rate"
slope_summary <- merge(slope_summary, ok_rates, by = c("phase", "variant"), sort = FALSE)
cat("\n=== Injection slope summary across replicates ===\n")
print(slope_summary)

recommend <- "inconclusive"
b_better_bias <- 0L
b_ok_variance <- 0L
n_cells <- 0L
for (dataset_id in dataset_ids) {
  for (phase_name in names(phase_specs)) {
    for (injection_a in injections) {
      n_cells <- n_cells + 1L
      a_row <- results[results$replicate == dataset_id &
        results$phase == phase_name & results$variant == "A_trim" &
        results$injection_a == injection_a, , drop = FALSE]
      b_row <- results[results$replicate == dataset_id &
        results$phase == phase_name & results$variant == "B_keep" &
        results$injection_a == injection_a, , drop = FALSE]
      if (nrow(a_row) != 1L || nrow(b_row) != 1L) next
      if (abs(b_row$beta0_bias) + 0.05 < abs(a_row$beta0_bias)) {
        b_better_bias <- b_better_bias + 1L
      }
      if (is.na(b_row$variance_flag) || b_row$variance_flag != "WARN_variance") {
        b_ok_variance <- b_ok_variance + 1L
      }
    }
  }
}
if (b_better_bias >= n_cells %/% 2 && b_ok_variance >= (2L * n_cells) %/% 3) {
  recommend <- "prefer_keep_constant"
} else if (all(results$bias_flag == "PASS", na.rm = TRUE) &&
    mean(results$variance_flag == "WARN_variance", na.rm = TRUE) < 0.5) {
  recommend <- "keep_trim_default"
} else {
  recommend <- "needs_constraint_or_case_by_case"
}

cat(sprintf("\n=== E1 recommendation (%d replicates): %s ===\n", length(dataset_ids), recommend))
cat(sprintf("B better bias in %d / %d cells; B variance OK in %d / %d cells\n",
            b_better_bias, n_cells, b_ok_variance, n_cells))
cat(sprintf("Wrote %s\n", out_csv))
cat(sprintf("Wrote %s\n", out_rep_csv))
cat(sprintf("Wrote %s\n", out_summary_csv))
