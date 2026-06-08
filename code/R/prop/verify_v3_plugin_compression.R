# V3: Does plug-in M compress posterior uncertainty beyond V2?
#
# Gaussian model: Y_t = X beta + v_t,  v_t ~ iid N(0, R_Matern),  sigma=1 known.
#
# Comparison:
#   Fixed R_hat : R_hat computed once from T_pre=500 pilot data (V2 scenario)
#   Adaptive    : Gibbs chain  beta -> R_hat(Y - X beta) -> beta
#                 R_hat updated from current residuals at EACH iteration
#
# Key metric per rho:
#   SD_adaptive(beta_j) / SD_fixed(beta_j)
#   > 1 : adaptive widens posterior (plug-in adds noise)
#   ~ 1 : plug-in M does not change posterior width (V2 dominates)
#   < 1 : plug-in M compresses uncertainty (concern)
#
# Analytical SD under fixed R_hat already available from V2 formula.

suppressMessages(library(autoFRK))
prop_dir <- "d:/Github/spatial-skew-t/code/R/prop"
for (f in c("prop_utils.R","prop_basis.R","prop_covariance.R","prop_modules.R"))
  source(file.path(prop_dir, f))

set.seed(20260601)

n      <- 50L; T_rep <- 30L; K <- 20L; nu <- 0.5
n_chain <- 30L; M <- 3000L; burn <- 500L
rhos   <- c(0.05, 0.15, 0.40)
p      <- 2L

loc    <- cbind(runif(n), runif(n))
X      <- cbind(1, loc[, 1])
beta_true <- c(2.0, 1.5)

matern_exp <- function(d, rho) exp(-d / rho)
d_mat <- as.matrix(dist(loc))

basis <- build_basis_matrix(s_obs = loc, k = K)
Fo    <- basis$F_obs

## analytical SD under fixed R_hat (uses V2 sandwich formula) ---------------
analytic_sd_fixed <- function(R_hat, R_M) {
  prec_hat <- chol2inv(prop_stable_chol(R_hat))
  XtWX     <- crossprod(X, prec_hat %*% X)
  XtWXinv  <- solve(XtWX)
  WX       <- prec_hat %*% X
  sandwich <- XtWXinv %*% crossprod(WX, R_M %*% WX) %*% XtWXinv
  M_W      <- prec_hat - WX %*% XtWXinv %*% t(WX)
  tr_M_R   <- sum(M_W * R_M)
  # nominal SE under known sigma=1
  se_nom   <- sqrt(diag(XtWXinv) / T_rep)
  # actual SD (sandwich under R_M)
  sd_true  <- sqrt(diag(sandwich) / T_rep)
  list(se_nom = se_nom, sd_true = sd_true)
}

## one Gibbs chain -----------------------------------------------------------
run_adaptive_chain <- function(Y_mat, Fo) {
  Ybar  <- rowMeans(Y_mat)
  beta  <- coef(lm(Ybar ~ X - 1))
  beta_chain <- matrix(NA, M, p)

  for (m in seq_len(M)) {
    # residuals under current beta
    resid <- Y_mat - as.vector(X %*% beta)

    # update R_hat via TH
    m_upd <- update_M_closed_form(Fo, resid, sigma_eps_sq = 0,
                                  sigma_floor = 1e-6)
    cs <- build_R_from_G(Fo, m_upd, F_pred = NULL)
    prec_curr <- chol2inv(prop_stable_chol(cs$r_obs))

    # GLS posterior for beta | R_hat, sigma=1 known
    XtWX <- crossprod(X, prec_curr %*% X)
    XtWy <- crossprod(X, prec_curr %*% Ybar)
    beta_h <- solve(XtWX, XtWy)
    V_beta <- solve(XtWX) / T_rep
    L_V    <- t(chol(prop_symmetrize(V_beta)))
    beta   <- beta_h + L_V %*% rnorm(p)
    beta_chain[m, ] <- beta
  }
  beta_chain[-(1:burn), ]
}

## main ----------------------------------------------------------------------
results <- vector("list", length(rhos))

for (ri in seq_along(rhos)) {
  rho  <- rhos[ri]
  R_M  <- matern_exp(d_mat, rho); diag(R_M) <- 1
  L_M  <- t(prop_stable_chol(R_M))

  # fixed R_hat from T_pre=500 pilot
  v_pilot  <- L_M %*% matrix(rnorm(n * 500L), n, 500L)
  m_pilot  <- update_M_closed_form(Fo, v_pilot, sigma_eps_sq = 0, sigma_floor = 1e-6)
  cs_fixed <- build_R_from_G(Fo, m_pilot)
  R_hat    <- cs_fixed$r_obs

  asd <- analytic_sd_fixed(R_hat, R_M)  # analytical SD under fixed R_hat

  # adaptive chains: n_chain independent datasets
  sd_adaptive <- matrix(NA, n_chain, p)
  for (nc in seq_len(n_chain)) {
    Y_mat <- as.vector(X %*% beta_true) + L_M %*% matrix(rnorm(n * T_rep), n, T_rep)
    ch    <- run_adaptive_chain(Y_mat, Fo)
    sd_adaptive[nc, ] <- apply(ch, 2, sd)
  }
  mean_sd_adapt <- colMeans(sd_adaptive)

  results[[ri]] <- list(
    rho             = rho,
    sd_fixed_nom    = asd$se_nom,       # nominal SD under fixed R_hat
    sd_fixed_true   = asd$sd_true,      # actual (sandwich) SD under fixed R_hat
    sd_adaptive     = mean_sd_adapt,    # empirical SD from adaptive chain
    compress_ratio  = mean_sd_adapt / asd$sd_true  # key metric
  )
  cat(sprintf("rho=%.2f done\n", rho))
}

## output -------------------------------------------------------------------
cat("\n=== V3: Plug-in M posterior compression ===\n")
cat(sprintf("n=%d  T=%d  K=%d  nu=%.1f  chains=%d  M=%d\n\n",
            n, T_rep, K, nu, n_chain, M))
cat(sprintf("%-6s %-6s | %-14s %-14s | %-14s | %-14s\n",
            "rho","param",
            "SD_fixed(true)","SD_adaptive",
            "compress_ratio","interpretation"))
cat(strrep("-", 80), "\n")
for (ri in seq_along(rhos)) {
  r <- results[[ri]]
  for (j in 1:p) {
    pname <- if (j==1) "beta_1" else "beta_2"
    ratio <- r$compress_ratio[j]
    interp <- if (ratio > 1.05) "wider (noise added)"
              else if (ratio < 0.95) "narrower (compressed)"
              else "unchanged"
    cat(sprintf("%-6.2f %-6s | %-14.5f %-14.5f | %-14.4f | %s\n",
                r$rho, pname,
                r$sd_fixed_true[j], r$sd_adaptive[j],
                ratio, interp))
  }
  cat("\n")
}
cat("=== VERDICT ===\n")
for (ri in seq_along(rhos)) {
  r <- results[[ri]]
  max_dev <- max(abs(r$compress_ratio - 1))
  cat(sprintf("rho=%.2f  max|compress_ratio-1|=%.3f  [%s]\n",
              r$rho, max_dev,
              ifelse(max_dev < 0.10, "PASS: V3 negligible", "WARN: V3 effect present")))
}
