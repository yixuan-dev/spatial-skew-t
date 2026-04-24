script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0L) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = FALSE)
  script_dir <- dirname(script_path)
  if (dir.exists(script_dir)) {
    setwd(script_dir)
  }
}

source("prop_utils.R")
source("prop_basis.R")
source("prop_covariance.R")
source("prop_modules.R")

expect_error_contains <- function(expr, pattern) {
  matched <- FALSE
  err_msg <- NULL
  tryCatch(
    force(expr),
    error = function(e) {
      err_msg <<- conditionMessage(e)
      matched <<- grepl(pattern, err_msg, fixed = TRUE)
    }
  )

  if (!matched) {
    stop(sprintf("Expected error containing '%s', got '%s'.", pattern, err_msg), call. = FALSE)
  }
}

run_test <- function(label, expr) {
  cat(sprintf("[TEST] %s\n", label))
  force(expr)
  cat(sprintf("[PASS] %s\n", label))
}

set.seed(20260424)

s_obs <- matrix(runif(24), ncol = 2)
s_pred <- matrix(runif(10), ncol = 2)
n_obs <- nrow(s_obs)
n_pred <- nrow(s_pred)
nt <- 4L
p <- 3L

y <- matrix(rnorm(n_obs * nt), nrow = n_obs, ncol = nt)
x_beta <- matrix(rnorm(n_obs * nt), nrow = n_obs, ncol = nt)
zg <- matrix(rnorm(n_obs * nt), nrow = n_obs, ncol = nt)
taug <- matrix(rexp(n_obs * nt, rate = 1) + 0.5, nrow = n_obs, ncol = nt)
x_pred <- array(rnorm(n_pred * nt * p), dim = c(n_pred, nt, p))
beta <- rnorm(p)
lambda <- 0.35

run_test("build_basis_matrix returns aligned basis blocks", {
  basis <- build_basis_matrix(s_obs = s_obs, s_pred = s_pred, k = 4L)
  stopifnot(
    is.matrix(basis$F_obs),
    is.matrix(basis$F_pred),
    nrow(basis$F_obs) == n_obs,
    nrow(basis$F_pred) == n_pred,
    ncol(basis$F_obs) == ncol(basis$F_pred),
    basis$rank_kept <= basis$requested_k
  )
})

run_test("build_basis_matrix rejects incompatible prediction coordinates", {
  bad_pred <- matrix(runif(15), ncol = 3)
  expect_error_contains(
    build_basis_matrix(s_obs = s_obs, s_pred = bad_pred, k = 4L),
    "s_pred must have 2 columns"
  )
})

run_test("build_basis_matrix rejects non-positive k", {
  expect_error_contains(
    build_basis_matrix(s_obs = s_obs, k = 0L),
    "k must be >="
  )
})

run_test("update_residuals returns raw and standardized residuals", {
  residual_block <- update_residuals(
    y = y,
    x_beta = x_beta,
    lambda = lambda,
    zg = zg,
    taug = taug
  )
  expected_raw <- y - x_beta - lambda * zg
  expected_std <- sqrt(taug) * expected_raw
  stopifnot(
    isTRUE(all.equal(residual_block$raw, expected_raw, tolerance = 1e-10)),
    isTRUE(all.equal(residual_block$std, expected_std, tolerance = 1e-10))
  )
})

run_test("update_residuals rejects non-positive taug", {
  bad_taug <- taug
  bad_taug[1, 1] <- 0
  expect_error_contains(
    update_residuals(y = y, x_beta = x_beta, taug = bad_taug),
    "taug must contain only positive values"
  )
})

run_test("update_M_closed_form returns symmetric M_hat and compatible rotation", {
  basis <- build_basis_matrix(s_obs = s_obs, s_pred = s_pred, k = 4L)
  residual_block <- update_residuals(y = y, x_beta = x_beta, lambda = lambda, zg = zg, taug = taug)
  m_update <- update_M_closed_form(
    F_obs = basis$F_obs,
    residuals_std = residual_block$std,
    F_pred = basis$F_pred,
    sigma_floor = 1e-6
  )
  stopifnot(
    is.matrix(m_update$rotation),
    is.matrix(m_update$M_hat),
    nrow(m_update$M_hat) == ncol(basis$F_obs),
    ncol(m_update$M_hat) == ncol(basis$F_obs),
    isTRUE(all.equal(m_update$M_hat, t(m_update$M_hat), tolerance = 1e-10)),
    m_update$effective_rank <= ncol(basis$F_obs)
  )
})

run_test("update_M_closed_form rejects row-mismatched residuals", {
  basis <- build_basis_matrix(s_obs = s_obs, s_pred = s_pred, k = 4L)
  expect_error_contains(
    update_M_closed_form(
      F_obs = basis$F_obs,
      residuals_std = matrix(rnorm((n_obs - 1L) * nt), nrow = n_obs - 1L),
      F_pred = basis$F_pred
    ),
    "residuals_std must have"
  )
})

run_test("build_R_from_G returns a unit-diagonal correlation matrix", {
  basis <- build_basis_matrix(s_obs = s_obs, s_pred = s_pred, k = 4L)
  residual_block <- update_residuals(y = y, x_beta = x_beta, lambda = lambda, zg = zg, taug = taug)
  m_update <- update_M_closed_form(F_obs = basis$F_obs, residuals_std = residual_block$std, F_pred = basis$F_pred)
  cov_state <- build_R_from_G(F_obs = basis$F_obs, M_update = m_update, F_pred = basis$F_pred)
  stopifnot(
    is.matrix(cov_state$r_obs),
    nrow(cov_state$r_obs) == n_obs,
    ncol(cov_state$r_obs) == n_obs,
    max(abs(diag(cov_state$r_obs) - 1)) < 1e-6
  )
})

run_test("build_R_from_G rejects incompatible eigenvalue length", {
  basis <- build_basis_matrix(s_obs = s_obs, s_pred = s_pred, k = 4L)
  residual_block <- update_residuals(y = y, x_beta = x_beta, taug = taug)
  m_update <- update_M_closed_form(F_obs = basis$F_obs, residuals_std = residual_block$std, F_pred = basis$F_pred)
  m_update$lowrank_eigs <- m_update$lowrank_eigs[-1]
  expect_error_contains(
    build_R_from_G(F_obs = basis$F_obs, M_update = m_update, F_pred = basis$F_pred),
    "M_update$lowrank_eigs must have length"
  )
})

run_test("quadform_logdet returns finite likelihood components", {
  basis <- build_basis_matrix(s_obs = s_obs, s_pred = s_pred, k = 4L)
  residual_block <- update_residuals(y = y, x_beta = x_beta, lambda = lambda, zg = zg, taug = taug)
  m_update <- update_M_closed_form(F_obs = basis$F_obs, residuals_std = residual_block$std, F_pred = basis$F_pred)
  cov_state <- build_R_from_G(F_obs = basis$F_obs, M_update = m_update, F_pred = basis$F_pred)
  ll_obj <- quadform_logdet(cov_state = cov_state, residuals_std = residual_block$std)
  stopifnot(all(is.finite(unlist(ll_obj))))
})

run_test("quadform_logdet rejects row-mismatched residuals", {
  basis <- build_basis_matrix(s_obs = s_obs, s_pred = s_pred, k = 4L)
  residual_block <- update_residuals(y = y, x_beta = x_beta, taug = taug)
  m_update <- update_M_closed_form(F_obs = basis$F_obs, residuals_std = residual_block$std, F_pred = basis$F_pred)
  cov_state <- build_R_from_G(F_obs = basis$F_obs, M_update = m_update, F_pred = basis$F_pred)
  expect_error_contains(
    quadform_logdet(
      cov_state = cov_state,
      residuals_std = matrix(rnorm((n_obs - 2L) * nt), nrow = n_obs - 2L)
    ),
    "residuals_std must have"
  )
})

run_test("predict_at_new_sites returns deterministic output when draw = FALSE", {
  basis <- build_basis_matrix(s_obs = s_obs, s_pred = s_pred, k = 4L)
  residual_block <- update_residuals(y = y, x_beta = x_beta, lambda = lambda, zg = zg, taug = taug)
  m_update <- update_M_closed_form(F_obs = basis$F_obs, residuals_std = residual_block$std, F_pred = basis$F_pred)
  cov_state <- build_R_from_G(F_obs = basis$F_obs, M_update = m_update, F_pred = basis$F_pred)
  z_pred <- matrix(rnorm(n_pred * nt), nrow = n_pred, ncol = nt)
  sigma_pred <- matrix(rexp(n_pred * nt, rate = 1) + 0.5, nrow = n_pred, ncol = nt)
  pred_a <- predict_at_new_sites(
    cov_state = cov_state,
    residuals_raw = residual_block$raw,
    taug = taug,
    x_pred = x_pred,
    beta = beta,
    lambda = lambda,
    z_pred = z_pred,
    sigma_pred = sigma_pred,
    draw = FALSE
  )
  pred_b <- predict_at_new_sites(
    cov_state = cov_state,
    residuals_raw = residual_block$raw,
    taug = taug,
    x_pred = x_pred,
    beta = beta,
    lambda = lambda,
    z_pred = z_pred,
    sigma_pred = sigma_pred,
    draw = FALSE
  )
  stopifnot(
    is.matrix(pred_a),
    nrow(pred_a) == n_pred,
    ncol(pred_a) == nt,
    isTRUE(all.equal(pred_a, pred_b, tolerance = 1e-10))
  )
})

run_test("predict_at_new_sites rejects non-positive sigma_pred", {
  basis <- build_basis_matrix(s_obs = s_obs, s_pred = s_pred, k = 4L)
  residual_block <- update_residuals(y = y, x_beta = x_beta, taug = taug)
  m_update <- update_M_closed_form(F_obs = basis$F_obs, residuals_std = residual_block$std, F_pred = basis$F_pred)
  cov_state <- build_R_from_G(F_obs = basis$F_obs, M_update = m_update, F_pred = basis$F_pred)
  bad_sigma_pred <- matrix(1, n_pred, nt)
  bad_sigma_pred[1, 1] <- 0
  expect_error_contains(
    predict_at_new_sites(
      cov_state = cov_state,
      residuals_raw = residual_block$raw,
      taug = taug,
      x_pred = x_pred,
      beta = beta,
      sigma_pred = bad_sigma_pred
    ),
    "sigma_pred must contain only positive values"
  )
})

cat("All prop module tests passed.\n")
