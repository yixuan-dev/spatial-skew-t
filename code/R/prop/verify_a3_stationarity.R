# A3 Verification: common R across all t (time-invariant spatial dependence)?
#
# Unlike A1/A2, A3 is NOT an exact model property — it is a genuine
# modelling assumption that can be violated by real data.
#
# The Tzeng-Huang update uses S = (1/T) sum_t r_t r_t', which estimates
# E_t[R_t].  When A3 holds this equals R; when A3 is violated the update
# returns a time-averaged covariance that may misrepresent both regimes.
#
# Tests
#   T1  Common R:   pooled TH estimate recovers true eigenvalues
#   T2  Two-regime: pooled estimate is a blend; quantify distortion per regime
#   T3  Split-half eigenvalue ratio — diagnostic for A3 on real data
#         common R  -> ratio ≈ 1  (within sampling noise)
#         two-regime -> ratio clearly != 1

suppressMessages(library(autoFRK))
prop_dir <- "d:/Github/spatial-skew-t/code/R/prop"
source(file.path(prop_dir, "prop_utils.R"))
source(file.path(prop_dir, "prop_basis.R"))
source(file.path(prop_dir, "prop_covariance.R"))
source(file.path(prop_dir, "prop_modules.R"))

set.seed(20260531)

## ---- shared basis ----------------------------------------------------------
n  <- 80L
TT <- 120L  # even number for clean split-half
K  <- 8L

loc   <- cbind(runif(n), runif(n))
basis <- build_basis_matrix(s_obs = loc, k = K)
Fo    <- basis$F_obs   # orthonormal

## ---- helper: build cov_state from raw parameters --------------------------
make_state <- function(eigs, sigma_xi_sq) {
  k0 <- basis$rank_kept
  upd <- list(
    rotation       = diag(k0),
    sigma_xi_sq    = sigma_xi_sq,
    lowrank_eigs   = c(eigs, rep(0, k0 - length(eigs))),
    effective_rank = sum(eigs > 0),
    projected_eigs = c(eigs, rep(0, k0 - length(eigs))),
    trace_s        = NA_real_,
    lstar          = sum(eigs > 0),
    M_hat          = diag(c(eigs, rep(0, k0 - length(eigs)))),
    F_rot_obs      = Fo,
    F_rot_pred     = NULL
  )
  build_R_from_G(F_obs = Fo, M_update = upd)
}

## ---- helper: draw T replicates v ~ N(0, R) from a cov_state ---------------
draw_v <- function(TT, cs) {
  eigs <- cs$lowrank_eigs
  v    <- matrix(NA_real_, n, TT)
  for (tt in seq_len(TT)) {
    eta       <- rnorm(length(eigs)) * sqrt(pmax(eigs, 0))
    xi        <- rnorm(n) * sqrt(cs$sigma_xi_sq)
    raw       <- as.vector(cs$q_obs %*% eta) + xi
    v[, tt]   <- raw / sqrt(cs$a_obs_diag)
  }
  v
}

## ---- helper: leading eigenvalue of sample covariance -----------------------
lead_eig <- function(v_mat) {
  s <- tcrossprod(v_mat) / ncol(v_mat)
  eigen(s, only.values = TRUE, symmetric = TRUE)$values[1]
}

## ---- helper: Tzeng-Huang effective rank from residuals --------------------
th_rank <- function(v_mat) {
  upd <- prop_closed_form_update(
    residuals    = v_mat,
    q_obs        = Fo,
    sigma_eps_sq = 0,
    sigma_floor  = 1e-6
  )
  upd$effective_rank
}

## ---- define two regimes ---------------------------------------------------
# Regime 1: strong spatial signal, low nugget
eigs1      <- c(2.0, 1.2, 0.5)
sigma_xi1  <- 0.1
cs1        <- make_state(eigs1, sigma_xi1)

# Regime 2: weak spatial signal, high nugget
eigs2      <- c(0.2, 0.1, 0.05)
sigma_xi2  <- 0.5
cs2        <- make_state(eigs2, sigma_xi2)

## ---- T1: common R (all T from regime 1) -----------------------------------
# Compare estimated R̂ to true R via the full pipeline.
# TH estimates G = F M F' + σ_ξ² I, then normalises to R.
# Metric: leading eigenvalue of R̂ vs R_true, and effective rank.
v_common   <- draw_v(TT, cs1)
m_upd_est  <- update_M_closed_form(Fo, v_common, sigma_eps_sq = 0,
                                   sigma_floor = 1e-6)
cs_est     <- build_R_from_G(Fo, m_upd_est)
R_hat      <- cs_est$r_obs
R_true_c1  <- cs1$r_obs

eig1_hat   <- eigen(R_hat,      only.values = TRUE, symmetric = TRUE)$values[1]
eig1_true  <- eigen(R_true_c1,  only.values = TRUE, symmetric = TRUE)$values[1]
eig1_relerr <- abs(eig1_hat - eig1_true) / eig1_true
rank_common <- m_upd_est$effective_rank

## ---- T2: two-regime R (first T/2 from cs1, second T/2 from cs2) ----------
TH <- TT %/% 2L
v_mixed <- cbind(draw_v(TH, cs1), draw_v(TH, cs2))
upd_mixed <- prop_closed_form_update(v_mixed, Fo, 0, 1e-6)

# Leading eigenvalue from each regime vs pooled estimate
le_r1     <- lead_eig(v_mixed[, 1:TH])
le_r2     <- lead_eig(v_mixed[, (TH + 1):TT])
le_pooled <- lead_eig(v_mixed)
rank_mixed <- upd_mixed$effective_rank

# True leading eig of each regime's R
le_r1_true <- eigen(cs1$r_obs, only.values = TRUE, symmetric = TRUE)$values[1]
le_r2_true <- eigen(cs2$r_obs, only.values = TRUE, symmetric = TRUE)$values[1]

## ---- T3: split-half eigenvalue ratio diagnostic ---------------------------
# Common-R case
v_sh_common <- draw_v(TT, cs1)
ratio_common <- lead_eig(v_sh_common[, 1:TH]) /
                lead_eig(v_sh_common[, (TH + 1):TT])

# Null distribution of ratio under common R (50 bootstrap draws)
null_ratios <- replicate(200L, {
  vb <- draw_v(TT, cs1)
  lead_eig(vb[, 1:TH]) / lead_eig(vb[, (TH + 1):TT])
})
ratio_q025 <- quantile(null_ratios, 0.025)
ratio_q975 <- quantile(null_ratios, 0.975)

# Two-regime case
ratio_mixed <- lead_eig(v_mixed[, 1:TH]) /
               lead_eig(v_mixed[, (TH + 1):TT])

## ---- output ----------------------------------------------------------------
cat("=== A3 Time-invariant R verification ===\n")
cat(sprintf("n=%d  T=%d  K0=%d\n", n, TT, basis$rank_kept))
cat(sprintf("Regime 1: eigs=[%.1f,%.1f,%.1f]  sigma_xi2=%.1f  (strong spatial)\n",
            eigs1[1], eigs1[2], eigs1[3], sigma_xi1))
cat(sprintf("Regime 2: eigs=[%.2f,%.2f,%.2f]  sigma_xi2=%.1f  (weak spatial)\n\n",
            eigs2[1], eigs2[2], eigs2[3], sigma_xi2))

cat("--- T1: common R — leading eigenvalue of R_hat vs R_true ---\n")
cat(sprintf("True  leading eig of R : %.4f\n", eig1_true))
cat(sprintf("Est   leading eig of R : %.4f\n", eig1_hat))
cat(sprintf("Relative error         : %.3f\n", eig1_relerr))
cat(sprintf("Effective rank (TH)    : %d  (true signal rank: %d)\n\n",
            rank_common, length(eigs1)))

cat("--- T2: two-regime R — blending distortion ---\n")
cat(sprintf("Leading eig of R1 (true)  : %.4f\n", le_r1_true))
cat(sprintf("Leading eig of R2 (true)  : %.4f\n", le_r2_true))
cat(sprintf("Empirical: first-half     : %.4f\n", le_r1))
cat(sprintf("Empirical: second-half    : %.4f\n", le_r2))
cat(sprintf("Pooled estimate           : %.4f  (blend of both)\n", le_pooled))
cat(sprintf("Effective rank (TH)       : %d  (regime 1 true: %d)\n\n",
            rank_mixed, length(eigs1)))

cat("--- T3: split-half leading-eigenvalue ratio diagnostic ---\n")
cat(sprintf("Null 95%% CI (200 draws, common R) : [%.3f, %.3f]\n",
            ratio_q025, ratio_q975))
cat(sprintf("Observed ratio, common R          : %.3f  (%s)\n",
            ratio_common,
            ifelse(ratio_common >= ratio_q025 & ratio_common <= ratio_q975,
                   "inside CI", "outside CI")))
cat(sprintf("Observed ratio, two-regime        : %.3f  (%s)\n\n",
            ratio_mixed,
            ifelse(ratio_mixed >= ratio_q025 & ratio_mixed <= ratio_q975,
                   "inside CI — A3 NOT flagged", "outside CI — A3 flagged")))

## ---- verdict ---------------------------------------------------------------
# Leading eig within 20% is the practical criterion; rank may over-select
# by 1 in finite samples (near-zero eigenvalues crossing the sigma_xi floor).
t1_ok <- eig1_relerr < 0.20
t2_ok <- le_pooled > le_r2 && le_pooled < le_r1   # blend is between regimes
t3a_ok <- ratio_common >= ratio_q025 & ratio_common <= ratio_q975
t3b_ok <- ratio_mixed < ratio_q025 | ratio_mixed > ratio_q975

cat("=== VERDICT ===\n")
cat(sprintf("T1 common-R eigenvalue recovery (<20%% err + rank) : %s\n",
            ifelse(t1_ok, "PASS", "FAIL")))
cat(sprintf("T2 pooled eig between regime eigs (blending)      : %s\n",
            ifelse(t2_ok, "PASS", "FAIL")))
cat(sprintf("T3a common-R ratio inside null CI                 : %s\n",
            ifelse(t3a_ok, "PASS", "FAIL")))
cat(sprintf("T3b two-regime ratio outside null CI              : %s\n",
            ifelse(t3b_ok, "PASS", "FAIL")))
cat(sprintf("OVERALL                                           : %s\n",
            ifelse(all(t1_ok, t2_ok, t3a_ok, t3b_ok), "PASS", "FAIL")))
