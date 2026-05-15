# ============================================================================
# Regression test for Fix 2 of the y_t update under the AR(2) prior.
#
# Under an AR(2) prior, the value y_t appears in three factors of the joint:
#   - self        : p(y_t     | y_{t-1}, y_{t-2})
#   - one ahead   : p(y_{t+1} | y_t,     y_{t-1})   -- y_t is lag-1
#   - two ahead   : p(y_{t+2} | y_{t+1}, y_t)       -- y_t is lag-2   <-- MISSING
# The current updateTauTS_AR2 / updateZTS_AR2 / updateKnotsTS_AR2 functions
# include only the first two factors in the MH ratio. Detailed balance then
# holds for a modified target that drops the third factor, so the chain's
# stationary distribution is NOT the AR(2) prior. In particular, the
# stationary lag-2 autocorrelation should equal gamma_2 = phi1*gamma1 + phi2;
# under the buggy update it does not.
#
# Test design: a prior-only y_t MH update (no data likelihood) on a single
# 1D chain. Run the BUGGY and PATCHED versions on the same RNG sequence and
# measure the empirical lag-1 and lag-2 autocorrelations averaged over many
# iterations. Assert:
#   - patched lag-2 ACF is close to gamma_2 (within tol)
#   - buggy   lag-2 ACF is FAR from gamma_2 (outside tol)
#   - patched lag-1 ACF is also close to gamma_1 (sanity check)
# ============================================================================

script_path <- (function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
        return(normalizePath(sub("^--file=", "", file_arg[1])))
    }
    return(normalizePath(sys.frames()[[1]]$ofile))
})()
ar2_dir <- normalizePath(file.path(dirname(script_path), ".."))
source(file.path(ar2_dir, "auxfunctions_ar2.R"))

# ---- Prior-only y_t MH update. include_lag2 toggles the t+2 correction.
updateY_AR2_prior <- function(y, phi1, phi2, mh, moments, include_lag2) {
    nt <- length(y)
    for (t in seq_len(nt)) {
        can <- rnorm(1, y[t], mh)
        lag1 <- if (t >= 2) y[t - 1] else NULL
        lag2 <- if (t >= 3) y[t - 2] else NULL

        # self
        R <- ar2_conditional_log_density(can,  t, lag1, lag2, phi1, phi2, moments) -
             ar2_conditional_log_density(y[t], t, lag1, lag2, phi1, phi2, moments)

        # t+1 (y_t is lag-1)
        if (t < nt) {
            next.lag2 <- if (t >= 2) y[t - 1] else NULL
            R <- R +
                 ar2_conditional_log_density(
                     y[t + 1], t + 1, can,  next.lag2, phi1, phi2, moments) -
                 ar2_conditional_log_density(
                     y[t + 1], t + 1, y[t], next.lag2, phi1, phi2, moments)
        }

        # t+2 (y_t is lag-2). Present only when Fix 2 is applied.
        if (include_lag2 && (t + 2) <= nt) {
            R <- R +
                 ar2_conditional_log_density(
                     y[t + 2], t + 2, y[t + 1], can,  phi1, phi2, moments) -
                 ar2_conditional_log_density(
                     y[t + 2], t + 2, y[t + 1], y[t], phi1, phi2, moments)
        }

        if (!is.na(R) && log(runif(1)) < R) y[t] <- can
    }
    y
}

# ---- Run a chain and compute average lag-h ACF across time and iterations.
run_chain_and_acf <- function(nt, phi1, phi2, mh, n.iter, n.burn,
                              include_lag2, seed) {
    moments <- get_ar2_standard_moments(phi1, phi2)
    set.seed(seed)
    y <- rnorm(nt)
    Y <- matrix(NA, n.iter, nt)
    for (iter in seq_len(n.iter)) {
        y <- updateY_AR2_prior(y, phi1, phi2, mh, moments, include_lag2)
        Y[iter, ] <- y
    }
    Y.post <- Y[(n.burn + 1):n.iter, , drop = FALSE]

    # Pool across iterations: marginal stats and pairwise covariances.
    means <- colMeans(Y.post)
    vars  <- apply(Y.post, 2, var)
    cov.h <- function(h) {
        idx <- seq_len(nt - h)
        sapply(idx, function(t) cov(Y.post[, t], Y.post[, t + h]))
    }
    acf1 <- mean(cov.h(1)) / mean(vars)
    acf2 <- mean(cov.h(2)) / mean(vars)
    list(mean = mean(means), var = mean(vars), acf1 = acf1, acf2 = acf2)
}

# ---- Truth: pick (phi1, phi2) so gamma_1 and gamma_2 are well separated
# so any distortion in the lag-2 structure is easy to detect.
phi1.t <- 0.7
phi2.t <- -0.4
moments.t <- get_ar2_standard_moments(phi1.t, phi2.t)

gamma0 <- moments.t$gamma0
gamma1 <- moments.t$gamma1
gamma2 <- moments.t$gamma2

cat(sprintf("Truth: phi1=%.2f, phi2=%.2f\n", phi1.t, phi2.t))
cat(sprintf("       gamma_0=%.4f, gamma_1=%.4f, gamma_2=%.4f\n",
            gamma0, gamma1, gamma2))
cat(sprintf("       innov_var=%.4f\n\n", moments.t$innov_var))

# ---- Common MCMC settings
nt     <- 80
mh     <- 0.5
n.iter <- 8000
n.burn <- 2000
seed   <- 4242

cat("Running prior-only chains (no data likelihood)...\n\n")

out.buggy <- run_chain_and_acf(nt, phi1.t, phi2.t, mh, n.iter, n.burn,
                                include_lag2 = FALSE, seed = seed)
out.patch <- run_chain_and_acf(nt, phi1.t, phi2.t, mh, n.iter, n.burn,
                                include_lag2 = TRUE,  seed = seed)

cat(sprintf("%-7s | %8s | %8s | %8s | %8s\n",
            "chain", "mean", "var", "acf1", "acf2"))
cat(strrep("-", 60), "\n", sep = "")
cat(sprintf("%-7s | %+8.4f | %8.4f | %+8.4f | %+8.4f\n",
            "target",  0,        gamma0,        gamma1,        gamma2))
cat(sprintf("%-7s | %+8.4f | %8.4f | %+8.4f | %+8.4f\n",
            "buggy",   out.buggy$mean, out.buggy$var,
            out.buggy$acf1, out.buggy$acf2))
cat(sprintf("%-7s | %+8.4f | %8.4f | %+8.4f | %+8.4f\n",
            "patched", out.patch$mean, out.patch$var,
            out.patch$acf1, out.patch$acf2))

# ---- Pass criteria
# 1. Patched lag-1 and lag-2 ACFs match the theoretical gamma values within tol
# 2. Buggy lag-2 ACF deviates from gamma_2 by more than tol
# 3. Patched chain is strictly closer to gamma_2 than buggy chain
tol.match <- 0.05
tol.bug   <- 0.05

err.patch.1 <- abs(out.patch$acf1 - gamma1)
err.patch.2 <- abs(out.patch$acf2 - gamma2)
err.bug.2   <- abs(out.buggy$acf2 - gamma2)

cat(sprintf("\nlag-1 ACF error  (patched) : %.4f  (tol %.3f)\n",
            err.patch.1, tol.match))
cat(sprintf("lag-2 ACF error  (patched) : %.4f  (tol %.3f)\n",
            err.patch.2, tol.match))
cat(sprintf("lag-2 ACF error  (buggy)   : %.4f  (must exceed %.3f)\n",
            err.bug.2, tol.bug))

pass <- (err.patch.1 < tol.match) &&
        (err.patch.2 < tol.match) &&
        (err.bug.2   > tol.bug)   &&
        (err.patch.2 < err.bug.2)

if (pass) {
    cat("\nPASS: lag-2 MH correction restores the stationary AR(2) prior.\n")
    quit(status = 0)
} else {
    cat("\nFAIL: either the patched chain does not match the target,\n")
    cat("       or the buggy chain does not visibly deviate.\n")
    quit(status = 1)
}
