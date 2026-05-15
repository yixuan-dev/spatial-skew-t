# ============================================================================
# Production-level regression test for Fix 2.
#
# The companion file test_fix2_lag2_mh.R verifies the math of the lag-2 MH
# correction via an inline mini-sampler -- but that test never executes the
# patched code in updateZTS_AR2 / updateTauTS_AR2 / updateKnotsTS_AR2. This
# test closes that gap by driving the real updateZTS_AR2 with a zero-data
# configuration so that the AR(2) prior alone governs the stationary
# distribution of z.star, then checking that the empirical lag-1 and lag-2
# autocorrelations of z.star match gamma_1 and gamma_2.
#
# Zero-data trick:
#   y       = 0     => observed residuals are zero
#   x.beta  = 0
#   lambda  = 0     => the data likelihood term -0.5 * quad.form(prec, ...) is
#                      identically 0 for both current and candidate proposals,
#                      and the |z| skew likelihood contribution vanishes too.
# Hence updateZTS_AR2's MH ratio reduces to the AR(2) prior contribution.
#
# Phi-update isolation:
#   updateZTS_AR2 ends with a call to updatePhiAR2TS. We stub that out with a
#   no-op that returns phi unchanged, so phi stays fixed at the truth and the
#   z.star equilibrium is governed only by the y_t update (i.e. Fix 2).
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

suppressPackageStartupMessages({
    library(emulator)        # provides quad.form
})
source(file.path(ar2_dir, "auxfunctions.R"))      # transform, hn.cop, ...
source(file.path(ar2_dir, "auxfunctions_ar2.R"))  # ar2_conditional_log_density
source(file.path(ar2_dir, "update_params_ar2.R")) # updateZTS_AR2, updatePhiAR2TS

# ---- Stub: prevent the internal updatePhiAR2TS call from drifting phi.
# updateZTS_AR2 calls updatePhiAR2TS(...) at the end of every invocation.
# We replace it with a no-op so phi stays exactly at the truth value passed
# in, isolating Fix 2's behaviour from Fix 3's.
updatePhiAR2TS <- function(data, phi1, phi2, day.mar,
                           att.phi1, acc.phi1, mh.phi1,
                           att.phi2, acc.phi2, mh.phi2,
                           prior.mean = 0, prior.sd = 0.5) {
    list(phi1 = phi1, phi2 = phi2,
         att.phi1 = att.phi1, acc.phi1 = acc.phi1,
         att.phi2 = att.phi2, acc.phi2 = acc.phi2)
}

# ---- Helper: build the zero-data input bundle for updateZTS_AR2.
make_inputs <- function(ns, nt, nknots, phi1, phi2, init.zstar = NULL, mh = 0.5) {
    # tau fixed at 1 => sig = 1; the half-normal copula transform is the
    # identity between standard half-normal and standard normal under |.|
    tau  <- matrix(1, nknots, nt)
    taug <- matrix(1, ns,     nt)

    # initial z (HN scale). If init.zstar is supplied, transform it back so
    # the chain starts on the stationary AR(2). Otherwise use a fixed value.
    if (is.null(init.zstar)) {
        z <- matrix(0.5, nknots, nt)
    } else {
        # z.star -> z via the inverse copula (sig = 1):
        z <- hn.invcop(x = init.zstar, sig = matrix(1, nknots, nt))
        z[z < 1e-6] <- 1e-6
    }
    g  <- matrix(1, ns, nt)
    zg <- matrix(0, ns, nt)
    for (t in seq_len(nt)) zg[, t] <- z[g[, t], t]

    list(
        z      = z,            zg     = zg,
        y      = matrix(0, ns, nt),
        lambda = 0,
        x.beta = matrix(0, ns, nt),
        phi1   = phi1,         phi2   = phi2,
        tau    = tau,          taug   = taug,
        g      = g,            prec   = diag(ns),
        acc    = matrix(0, nknots, nt),
        att    = matrix(0, nknots, nt),
        mh     = matrix(mh, nknots, nt),
        acc.phi1 = 0, att.phi1 = 0, mh.phi1 = 0.1,
        acc.phi2 = 0, att.phi2 = 0, mh.phi2 = 0.1
    )
}

# ---- One chain: returns the n.iter x nt matrix of z.star samples.
run_zstar_chain <- function(ns, nt, nknots, phi1, phi2, mh, n.iter, n.burn,
                            seed) {
    set.seed(seed)
    moments <- get_ar2_standard_moments(phi1, phi2)

    # Initialise z.star from the stationary AR(2) to skip transient.
    init.zstar <- simulate_ar2_standard(nt = nt, n = nknots,
                                         phi1 = phi1, phi2 = phi2)
    st <- make_inputs(ns, nt, nknots, phi1, phi2,
                      init.zstar = init.zstar, mh = mh)

    Zstar <- matrix(NA_real_, n.iter, nt)
    sig.mat <- matrix(1, nknots, nt)   # tau = 1 => sig = 1 throughout

    for (iter in seq_len(n.iter)) {
        out <- updateZTS_AR2(
            z = st$z,           zg = st$zg,
            y = st$y,           lambda = st$lambda, x.beta = st$x.beta,
            phi1 = phi1,        phi2 = phi2,
            tau = st$tau,       taug = st$taug,
            g = st$g,           prec = st$prec,
            acc = st$acc,       att = st$att, mh = st$mh,
            acc.phi1 = st$acc.phi1, att.phi1 = st$att.phi1,
            mh.phi1 = st$mh.phi1,
            acc.phi2 = st$acc.phi2, att.phi2 = st$att.phi2,
            mh.phi2 = st$mh.phi2
        )
        st$z  <- out$z
        st$zg <- out$zg
        st$acc <- out$acc;       st$att <- out$att
        st$acc.phi1 <- out$acc.phi1; st$att.phi1 <- out$att.phi1
        st$acc.phi2 <- out$acc.phi2; st$att.phi2 <- out$att.phi2

        # Record z.star = hn.cop(z, sig = 1) for the (sole) knot.
        zstar.mat <- hn.cop(x = st$z, sig = sig.mat)
        Zstar[iter, ] <- zstar.mat[1, ]
    }

    Zstar[(n.burn + 1):n.iter, , drop = FALSE]
}

# ---- Pool ACFs across time and iterations.
acf_summary <- function(Z) {
    nt <- ncol(Z)
    means <- colMeans(Z)
    vars  <- apply(Z, 2, var)
    cov.h <- function(h) {
        idx <- seq_len(nt - h)
        sapply(idx, function(t) cov(Z[, t], Z[, t + h]))
    }
    list(
        mean = mean(means), var = mean(vars),
        acf1 = mean(cov.h(1)) / mean(vars),
        acf2 = mean(cov.h(2)) / mean(vars)
    )
}

# ===========================================================================
# Truth: same (phi1, phi2) as the math test, with well-separated gamma_1, gamma_2
# ===========================================================================
phi1.t <- 0.7
phi2.t <- -0.4
mom <- get_ar2_standard_moments(phi1.t, phi2.t)

cat(sprintf("Truth: phi1=%.2f, phi2=%.2f  =>  gamma_1=%.4f, gamma_2=%.4f\n\n",
            phi1.t, phi2.t, mom$gamma1, mom$gamma2))

# Smaller MCMC than the math test: updateZTS_AR2 is heavier (full nks*nt loop
# per iteration with the copula transform), so we use fewer iterations and
# more time points to keep the wall clock reasonable.
ns     <- 5
nt     <- 60
nknots <- 1
mh     <- 0.5
n.iter <- 3000
n.burn <- 500

cat(sprintf("Driving updateZTS_AR2 with zero data (ns=%d, nt=%d, nknots=%d)...\n\n",
            ns, nt, nknots))

t0 <- Sys.time()
Z <- run_zstar_chain(ns, nt, nknots, phi1.t, phi2.t, mh, n.iter, n.burn,
                     seed = 4242)
t1 <- Sys.time()
out <- acf_summary(Z)

cat(sprintf("%-7s | %8s | %8s | %8s | %8s\n",
            "chain", "mean", "var", "acf1", "acf2"))
cat(strrep("-", 60), "\n", sep = "")
cat(sprintf("%-7s | %+8.4f | %8.4f | %+8.4f | %+8.4f\n",
            "target", 0, 1, mom$gamma1, mom$gamma2))
cat(sprintf("%-7s | %+8.4f | %8.4f | %+8.4f | %+8.4f\n",
            "real",   out$mean, out$var, out$acf1, out$acf2))

cat(sprintf("\nWall: %.1fs\n", as.numeric(difftime(t1, t0, units = "secs"))))

# Pass criteria (uncalibrated thresholds; will be tightened in a later
# Monte Carlo pass).
tol.match <- 0.05

err1 <- abs(out$acf1 - mom$gamma1)
err2 <- abs(out$acf2 - mom$gamma2)
cat(sprintf("\nlag-1 ACF error (real updateZTS_AR2): %.4f  (tol %.3f)\n",
            err1, tol.match))
cat(sprintf("lag-2 ACF error (real updateZTS_AR2): %.4f  (tol %.3f)\n",
            err2, tol.match))

pass <- err1 < tol.match && err2 < tol.match
if (pass) {
    cat("\nPASS: updateZTS_AR2 reproduces the stationary AR(2) prior.\n")
    quit(status = 0)
} else {
    cat("\nFAIL: production updateZTS_AR2 does not reproduce the AR(2) prior.\n")
    cat("       If this fails on patched code, Fix 2 has been regressed in\n")
    cat("       update_params_ar2.R; re-add the t+2 correction block.\n")
    quit(status = 1)
}
