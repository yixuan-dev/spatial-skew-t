# Numerical cross-check: prop closed-form update vs autoFRK cMLEimat (s = 0)
# Verifies that prop_closed_form_update / update_M_closed_form reproduce the
# Tzeng-Huang estimator as implemented in autoFRK, on basis-invariant quantities:
#   - noise variance  v (sigma_xi^2)
#   - implied low-rank covariance  G = F M F' + v I   (n x n)
#   - retained eigenvalue spectrum
#
# M itself differs by a basis transform (autoFRK keeps M in the raw mrts basis,
# prop orthonormalises the basis so F'F = I), so we compare G, not M directly.

suppressMessages({
  library(autoFRK)
})

prop_dir <- "d:/Github/spatial-skew-t/code/R/prop"
source(file.path(prop_dir, "prop_utils.R"))
source(file.path(prop_dir, "prop_basis.R"))
source(file.path(prop_dir, "prop_covariance.R"))
source(file.path(prop_dir, "prop_modules.R"))

set.seed(20260531)

## ---- synthetic data -------------------------------------------------------
n  <- 60L      # spatial locations
TT <- 40L      # replicates (T > K for identifiability)
K  <- 8L       # basis dimension

loc <- cbind(runif(n), runif(n))           # 2-D locations in unit square

# Common orthonormalised basis fed to BOTH estimators, so the comparison
# isolates the estimator (F'F = I removes the (F'F)^{-1/2} difference and any
# basis-construction divergence such as prop trimming the constant column).
basis  <- build_basis_matrix(s_obs = loc, k = K)
Fo     <- basis$F_obs                        # n x K0, orthonormal (F'F = I)
K0     <- ncol(Fo)

# low-rank signal (rank 3) + iid noise, so L* > 0 is non-trivial
rank_true <- 3L
eta <- matrix(rnorm(K0 * TT), K0, TT)
eta[(rank_true + 1L):K0, ] <- 0
signal <- Fo %*% eta
noise  <- matrix(rnorm(n * TT, sd = 0.5), n, TT)
data   <- signal + noise                    # n x T

## ---- autoFRK estimator (s = 0, identity D, no missing) --------------------
auto <- cMLEimat(Fo, data, s = 0, wSave = TRUE)
# cMLEimat returns out$s = INPUT noise floor, out$v = ESTIMATED noise; the
# total noise variance in the fitted model is s + v (indeMLE later renames
# v -> s when the input s is 0).
v_auto <- auto$s + auto$v                    # estimated total noise variance
M_auto <- auto$M
G_auto <- Fo %*% M_auto %*% t(Fo) + diag(v_auto, n)

## ---- prop estimator -------------------------------------------------------
upd    <- update_M_closed_form(F_obs = Fo, residuals_std = data,
                               sigma_eps_sq = 0, sigma_floor = 0)
v_prop <- upd$sigma_xi_sq
M_prop <- upd$M_hat
# implied G in prop's (rotated, orthonormal) basis
Fr     <- upd$F_rot_obs
G_prop <- Fr %*% diag(upd$lowrank_eigs, length(upd$lowrank_eigs)) %*% t(Fr) +
          diag(v_prop, n)

## ---- comparisons ----------------------------------------------------------
cat("=== Tzeng-Huang estimator: prop vs autoFRK (s = 0) ===\n\n")
cat(sprintf("n = %d,  T = %d,  K0 = %d (common orthonormal basis),  rank = %d\n\n",
            n, TT, K0, rank_true))

cat(sprintf("noise variance v:   autoFRK = %.10f\n", v_auto))
cat(sprintf("                    prop    = %.10f\n", v_prop))
cat(sprintf("                    |diff|  = %.3e\n\n", abs(v_auto - v_prop)))

cat(sprintf("retained eigs autoFRK : %s\n",
            paste(sprintf("%.4f", sort(eigen(M_auto, only.values = TRUE)$values,
                  decreasing = TRUE)), collapse = ", ")))
cat(sprintf("retained eigs prop    : %s\n\n",
            paste(sprintf("%.4f", sort(upd$lowrank_eigs, decreasing = TRUE)),
                  collapse = ", ")))

m_absdiff <- max(abs(M_auto - M_prop))
m_reldiff <- m_absdiff / max(abs(M_auto))
cat(sprintf("M (same basis, %d x %d):\n", K0, K0))
cat(sprintf("   max |M_auto - M_prop|        = %.3e\n", m_absdiff))
cat(sprintf("   max |M_auto - M_prop| / |M|  = %.3e\n\n", m_reldiff))

g_absdiff <- max(abs(G_auto - G_prop))
g_reldiff <- g_absdiff / max(abs(G_auto))
cat(sprintf("implied covariance G = F M F' + v I  (%d x %d):\n", n, n))
cat(sprintf("   max |G_auto - G_prop|        = %.3e\n", g_absdiff))
cat(sprintf("   max |G_auto - G_prop| / |G|  = %.3e\n\n", g_reldiff))

ok_v <- abs(v_auto - v_prop) < 1e-7
ok_M <- m_reldiff < 1e-7
ok_G <- g_reldiff < 1e-7
cat("=== VERDICT ===\n")
cat(sprintf("noise variance match : %s\n", ifelse(ok_v, "PASS", "FAIL")))
cat(sprintf("M (same basis) match : %s\n", ifelse(ok_M, "PASS", "FAIL")))
cat(sprintf("implied G match      : %s\n", ifelse(ok_G, "PASS", "FAIL")))
cat(sprintf("OVERALL              : %s\n",
            ifelse(ok_v && ok_M && ok_G, "PASS", "FAIL")))
