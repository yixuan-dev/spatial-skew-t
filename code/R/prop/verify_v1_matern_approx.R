# V1: MRTS approximation quality for Matern spatial covariance
#
# Core question: how well does R_MRTS(K) approximate R_Matern(rho, nu)?
#
# Model:  v_t ~ iid N(0, R_Matern).  Feed to TH -> R_hat(K).
# Metrics (all relative or per-observation to be scale-free):
#   M1  ||R̂ - R_M||_F  /  ||R_M||_F              overall approximation
#   M2  |lambda_1(R̂) - lambda_1(R_M)| / lambda_1  dominant eigenvalue
#   M3  |log|R̂| - log|R_M|| / n                  per-obs log-det error
#   M4  mean_t |v_t'R̂^{-1}v_t - v_t'R_M^{-1}v_t| / n  per-obs QF error
#
# Baseline: R = I (independence), gives upper bound on M1-M4.
# Varies: K in {5,10,20,40}, nu in {0.5,1.5,2.5}, rho in {0.05,0.15,0.40}
# Fixed:  n=100 sites in [0,1]^2, T=500 oracle samples, same seed throughout.

suppressMessages(library(autoFRK))
prop_dir <- "d:/Github/spatial-skew-t/code/R/prop"
source(file.path(prop_dir, "prop_utils.R"))
source(file.path(prop_dir, "prop_basis.R"))
source(file.path(prop_dir, "prop_covariance.R"))
source(file.path(prop_dir, "prop_modules.R"))

set.seed(20260601)

## ---- Matern correlation function ------------------------------------------
matern_corr <- function(d, rho, nu) {
  d <- pmax(d, 0)
  if (nu == 0.5) return(exp(-d / rho))
  if (nu == 1.5) { x <- sqrt(3) * d / rho; return((1 + x) * exp(-x)) }
  if (nu == 2.5) { x <- sqrt(5) * d / rho; return((1 + x + x^2/3) * exp(-x)) }
  x <- pmax(sqrt(2*nu) * d / rho, 1e-10)
  (2^(1-nu)/gamma(nu)) * x^nu * besselK(x, nu)
}

build_matern <- function(loc, rho, nu) {
  d <- as.matrix(dist(loc))
  R <- matern_corr(d, rho, nu)
  diag(R) <- 1
  prop_symmetrize(R)
}

## ---- fixed design ----------------------------------------------------------
n        <- 100L
TT       <- 500L
K_vals   <- c(5L, 10L, 20L, 40L)
nu_vals  <- c(0.5, 1.5, 2.5)
rho_vals <- c(0.05, 0.15, 0.40)

loc <- cbind(runif(n), runif(n))   # fixed across all (K, nu, rho)

## ---- precompute MRTS bases (one per K) ------------------------------------
bases <- lapply(K_vals, function(K) {
  build_basis_matrix(s_obs = loc, k = K)
})
names(bases) <- paste0("K", K_vals)

## ---- metrics ---------------------------------------------------------------
compute_metrics <- function(R_true, cov_state, v_oracle) {
  R_hat   <- cov_state$r_obs
  n_loc   <- nrow(R_true)

  # M1: Frobenius relative error
  m1 <- norm(R_hat - R_true, "F") / norm(R_true, "F")

  # M2: leading eigenvalue relative error
  e_true <- eigen(R_true, only.values = TRUE, symmetric = TRUE)$values[1]
  e_hat  <- eigen(R_hat,  only.values = TRUE, symmetric = TRUE)$values[1]
  m2 <- abs(e_hat - e_true) / e_true

  # M3: |log|R̂| - log|R_M|| / n
  logdet_hat  <- cov_state$logdet_r_obs
  logdet_true <- sum(log(pmax(
    eigen(R_true, only.values = TRUE, symmetric = TRUE)$values, 1e-15)))
  m3 <- abs(logdet_hat - logdet_true) / n_loc

  # M4: mean_t |v'R̂^{-1}v - v'R_true^{-1}v| / n
  prec_true <- chol2inv(prop_stable_chol(R_true))
  qf_true   <- colSums(v_oracle * (prec_true %*% v_oracle))
  rinv_hat  <- prop_apply_rinv_obs(cov_state, v_oracle)
  qf_hat    <- colSums(v_oracle * rinv_hat)
  m4 <- mean(abs(qf_hat - qf_true)) / n_loc

  c(M1 = m1, M2 = m2, M3 = m3, M4 = m4)
}

## ---- independence baseline -------------------------------------------------
indep_metrics <- function(R_true, v_oracle) {
  n_loc     <- nrow(R_true)
  R_hat     <- diag(n_loc)
  m1 <- norm(R_hat - R_true, "F") / norm(R_true, "F")
  e_true <- eigen(R_true, only.values = TRUE, symmetric = TRUE)$values[1]
  m2 <- abs(1 - e_true) / e_true
  logdet_true <- sum(log(pmax(
    eigen(R_true, only.values = TRUE, symmetric = TRUE)$values, 1e-15)))
  m3 <- abs(0 - logdet_true) / n_loc      # log|I| = 0
  prec_true <- chol2inv(prop_stable_chol(R_true))
  qf_true   <- colSums(v_oracle * (prec_true %*% v_oracle))
  qf_indep  <- colSums(v_oracle^2)         # R=I => R^{-1}=I
  m4 <- mean(abs(qf_indep - qf_true)) / n_loc
  c(M1 = m1, M2 = m2, M3 = m3, M4 = m4)
}

## ---- main loop ------------------------------------------------------------
results <- vector("list", length(nu_vals) * length(rho_vals) *
                    (length(K_vals) + 1L))
idx <- 1L

for (nu in nu_vals) {
  for (rho in rho_vals) {
    R_matern <- build_matern(loc, rho, nu)
    L_chol   <- t(prop_stable_chol(R_matern))    # lower Cholesky for sampling
    v_oracle <- L_chol %*% matrix(rnorm(n * TT), n, TT)

    # Independence baseline
    bm <- indep_metrics(R_matern, v_oracle)
    results[[idx]] <- data.frame(nu = nu, rho = rho, K = 0L,
                                 label = "indep",
                                 M1 = bm["M1"], M2 = bm["M2"],
                                 M3 = bm["M3"], M4 = bm["M4"],
                                 row.names = NULL)
    idx <- idx + 1L

    for (ki in seq_along(K_vals)) {
      K  <- K_vals[ki]
      bs <- bases[[ki]]
      Fo <- bs$F_obs
      K0 <- bs$rank_kept

      m_upd <- update_M_closed_form(
        F_obs         = Fo,
        residuals_std = v_oracle,
        sigma_eps_sq  = 0,
        sigma_floor   = 1e-6
      )
      cs <- build_R_from_G(Fo, m_upd)
      mt <- compute_metrics(R_matern, cs, v_oracle)

      results[[idx]] <- data.frame(nu = nu, rho = rho, K = K,
                                   label = paste0("K=", K),
                                   M1 = mt["M1"], M2 = mt["M2"],
                                   M3 = mt["M3"], M4 = mt["M4"],
                                   row.names = NULL)
      idx <- idx + 1L
    }
  }
}
results <- do.call(rbind, results)

## ---- print results --------------------------------------------------------
cat("=== V1: MRTS approximation quality for Matern covariance ===\n")
cat(sprintf("n=%d  T=%d  locations: uniform [0,1]^2\n\n", n, TT))

for (nu in nu_vals) {
  cat(sprintf("--- nu = %.1f ---\n", nu))
  cat(sprintf("%-8s %-6s | %-8s %-8s %-8s %-8s\n",
              "rho", "K", "M1(Frob)", "M2(eig1)", "M3(ldet)", "M4(QF)"))
  cat(strrep("-", 58), "\n")
  sub <- results[results$nu == nu, ]
  for (rho in rho_vals) {
    rows <- sub[sub$rho == rho, ]
    rows <- rows[order(rows$K), ]
    for (i in seq_len(nrow(rows))) {
      r <- rows[i, ]
      lbl <- if (r$K == 0) "indep" else paste0("K=", r$K)
      cat(sprintf("%-8.2f %-6s | %-8.4f %-8.4f %-8.4f %-8.4f\n",
                  rho, lbl, r$M1, r$M2, r$M3, r$M4))
    }
    cat("\n")
  }
}

## ---- summary: K needed for M1 < 0.10 per (nu, rho) -----------------------
cat("=== K needed for M1 (Frobenius) < 0.10 ===\n")
cat(sprintf("%-6s | %-8s %-8s %-8s\n", "nu", "rho=0.05", "rho=0.15", "rho=0.40"))
cat(strrep("-", 36), "\n")
for (nu in nu_vals) {
  row_str <- sprintf("%-6.1f |", nu)
  for (rho in rho_vals) {
    sub <- results[results$nu == nu & results$rho == rho & results$K > 0, ]
    sub <- sub[order(sub$K), ]
    k_thresh <- sub$K[which(sub$M1 < 0.10)[1]]
    row_str <- paste0(row_str, sprintf(" %-8s",
                      ifelse(is.na(k_thresh), ">40", paste0("K=", k_thresh))))
  }
  cat(row_str, "\n")
}
