# ============================================================================
# Regression test: the phi-block and the latent-blocks of the AR(2) sampler
# must target the SAME joint density p(Y | phi) (compatible Gibbs blocks).
#
# History
# -------
# updatePhiAR2TS originally evaluated only the Box-Jenkins CONDITIONAL
# likelihood (transitions t >= 3, via ar2_transition_loglik), inherited from
# the AR(1) code where skipping the initial term is an exact identity: for
# p = 1 the initial factor is N(Y_1; 0, gamma_0) with gamma_0 = 1 pinned,
# hence phi-free. For AR(2) the initial block contains
#     p(Y_2 | Y_1, phi) = N(gamma_1 Y_1, 1 - gamma_1^2),
#     gamma_1 = phi1 / (1 - phi2),
# which IS phi-dependent. The latent-block updaters (updateZTS_AR2 etc.)
# include this term (the Fix 1 t = 2 branch of get_ar2_conditional_params),
# so the two blocks targeted different distributions -- an incompatible
# (pseudo-)Gibbs sampler whose composite chain is invariant for neither.
#
# The fix (2026-07-15) replaced ar2_transition_loglik in updatePhiAR2TS with
# ar2_stationary_loglik, which adds the phi-dependent initial-block factor.
# A replicated Geweke prior-recovery experiment (T = 4/6, up to 6 chains of
# 40k-120k sweeps) showed the conditional-only version deviates from the
# truncated prior in E[phi2] (z = 1.8-2.1, same sign in two independent
# designs) while the stationary version passes all checks (|z| <= 1).
# Full derivation and measurements:
#     tex/ar2_phi_exact_likelihood/ar2_phi_exact_likelihood.tex
#
# Layers
# ------
#   (1) STRUCTURAL: the body of updatePhiAR2TS calls ar2_stationary_loglik
#       and contains no live call to ar2_transition_loglik. Strong guard
#       against re-introduction of the conditional-only likelihood.
#
#   (2) EXACT IDENTITY (theorem-level, no MCMC): for random stationary
#       (phi1, phi2) and random data,
#           ar2_stationary_loglik - ar2_transition_loglik
#             == sum_k ar2_conditional_log_density(y2[k], t = 2, lag1 = y1[k])
#       to 1e-10. This pins the added term to the EXACT density the
#       latent-blocks use at t = 2 -- block compatibility, term by term.
#
#   (3) GEWEKE SMOKE GATE: short replicated prior-recovery run through the
#       production updatePhiAR2TS (likelihood off, latent prior-only Y-block
#       mirroring updateZTS_AR2's t / t+1 / t+2 structure). With the fix the
#       phi-marginal must recover the N(0, 0.5^2)^2 prior truncated to the
#       stationarity triangle. Gate |z| < 4 with fixed seeds: catches gross
#       breakage deterministically; the subtle conditional-only bug needs
#       the heavier offline design above (documented, not re-run here).
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
src_path <- file.path(ar2_dir, "update_params_ar2.R")

source(file.path(ar2_dir, "auxfunctions_ar2.R"))
source(src_path)

extract_fn_lines <- function(src.lines, fn.name) {
    start <- grep(paste0("^", fn.name, "\\s*<-\\s*function"), src.lines)
    if (length(start) == 0) {
        stop(sprintf("could not locate function %s", fn.name))
    }
    start <- start[1]
    depth <- 0
    end <- NA_integer_
    seen.open <- FALSE
    count_char <- function(ch, s) nchar(gsub(paste0("[^", ch, "]"), "", s))

    for (i in seq(start, length(src.lines))) {
        ln <- sub("#.*$", "", src.lines[i])
        opens <- count_char("\\{", ln)
        closes <- count_char("\\}", ln)
        if (opens > 0) seen.open <- TRUE
        depth <- depth + opens - closes
        if (seen.open && depth == 0) {
            end <- i
            break
        }
    }
    if (is.na(end)) {
        stop(sprintf("could not find matching closing brace for %s", fn.name))
    }
    src.lines[start:end]
}

# ---------------------------------------------------------------------------
# (1) Structural assertion
# ---------------------------------------------------------------------------

src <- readLines(src_path)
phi.body <- sub("#.*$", "", extract_fn_lines(src, "updatePhiAR2TS"))

has.stationary <- any(grepl("ar2_stationary_loglik\\s*\\(", phi.body))
has.transition <- any(grepl("ar2_transition_loglik\\s*\\(", phi.body))
struct.ok <- has.stationary && !has.transition

if (!struct.ok) {
    stop(sprintf(
        paste0("STRUCTURAL FAIL: updatePhiAR2TS must call ",
               "ar2_stationary_loglik (found: %s) and must not call ",
               "ar2_transition_loglik (found: %s). The conditional-only ",
               "likelihood makes the phi-block target a different ",
               "distribution than the latent-blocks; see header."),
        has.stationary, has.transition
    ))
}
cat("  [1/3] structural: PASS\n")

# ---------------------------------------------------------------------------
# (2) Exact identity: added term == latent-block t = 2 density
# ---------------------------------------------------------------------------

set.seed(314)
max.err <- 0
for (r in 1:200) {
    repeat {
        phi1 <- runif(1, -1.9, 1.9)
        phi2 <- runif(1, -0.95, 0.95)
        if (check_ar2_stability(phi1, phi2)) break
    }
    K <- sample(1:5, 1)
    Tt <- sample(4:12, 1)
    Y <- matrix(rnorm(K * Tt), K, Tt)
    m <- get_ar2_standard_moments(phi1, phi2)

    cur <- Y[, 3:Tt, drop = FALSE]
    l1 <- Y[, 2:(Tt - 1), drop = FALSE]
    l2 <- Y[, 1:(Tt - 2), drop = FALSE]

    lhs <- ar2_stationary_loglik(cur, l1, l2, Y[, 1], Y[, 2],
                                 phi1 = phi1, phi2 = phi2, moments = m) -
        ar2_transition_loglik(cur, l1, l2,
                              phi1 = phi1, phi2 = phi2, moments = m)
    rhs <- sum(vapply(1:K, function(k) {
        ar2_conditional_log_density(Y[k, 2], t = 2, lag1 = Y[k, 1],
                                    phi1 = phi1, phi2 = phi2, moments = m)
    }, numeric(1)))
    max.err <- max(max.err, abs(lhs - rhs))
}
if (max.err > 1e-10) {
    stop(sprintf(
        paste0("IDENTITY FAIL: the initial-block term added by ",
               "ar2_stationary_loglik differs from the latent-block t = 2 ",
               "density by %.3e; the two Gibbs blocks no longer share one ",
               "joint density."), max.err
    ))
}
cat(sprintf("  [2/3] exact identity (200 random draws): PASS (max err %.2e)\n",
            max.err))

# ---------------------------------------------------------------------------
# (3) Geweke smoke gate (prior-recovery through the production phi-block)
# ---------------------------------------------------------------------------

Tt <- 4
K <- 3

update_Y_prior <- function(Y, phi1, phi2, mh = 1.0) {
    m <- get_ar2_standard_moments(phi1, phi2)
    for (t in 1:Tt) {
        for (k in 1:K) {
            can <- rnorm(1, Y[k, t], mh)
            lag1 <- if (t >= 2) Y[k, t - 1] else NULL
            lag2 <- if (t >= 3) Y[k, t - 2] else NULL
            R <- ar2_conditional_log_density(can, t, lag1, lag2, phi1, phi2, m) -
                ar2_conditional_log_density(Y[k, t], t, lag1, lag2, phi1, phi2, m)
            if (t < Tt) {
                nl2 <- if (t >= 2) Y[k, t - 1] else NULL
                R <- R +
                    ar2_conditional_log_density(Y[k, t + 1], t + 1, can, nl2, phi1, phi2, m) -
                    ar2_conditional_log_density(Y[k, t + 1], t + 1, Y[k, t], nl2, phi1, phi2, m)
            }
            if ((t + 2) <= Tt) {
                R <- R +
                    ar2_conditional_log_density(Y[k, t + 2], t + 2, Y[k, t + 1], can, phi1, phi2, m) -
                    ar2_conditional_log_density(Y[k, t + 2], t + 2, Y[k, t + 1], Y[k, t], phi1, phi2, m)
            }
            if (!is.na(R) && log(runif(1)) < R) Y[k, t] <- can
        }
    }
    Y
}

run_chain <- function(seed, iters = 30000, burn = 3000, mh.phi = 0.5) {
    set.seed(seed)
    phi1 <- 0
    phi2 <- 0
    Y <- simulate_ar2_standard(Tt, K, 0, 0)
    keep <- matrix(NA_real_, iters, 2)
    for (i in 1:iters) {
        Y <- update_Y_prior(Y, phi1, phi2)
        up <- updatePhiAR2TS(
            data = Y, phi1 = phi1, phi2 = phi2, day.mar = 2,
            att.phi1 = 0, acc.phi1 = 0, mh.phi1 = mh.phi,
            att.phi2 = 0, acc.phi2 = 0, mh.phi2 = mh.phi
        )
        phi1 <- up$phi1
        phi2 <- up$phi2
        keep[i, ] <- c(phi1, phi2)
    }
    colMeans(keep[(burn + 1):iters, ])
}

# truncated-prior reference moments by rejection sampling
set.seed(99)
ref <- matrix(rnorm(1.2e6, 0, 0.5), ncol = 2)
ref <- ref[mapply(check_ar2_stability, ref[, 1], ref[, 2]), , drop = FALSE]
ref.mean <- colMeans(ref)

nrep <- 4
chains <- t(sapply(1:nrep, function(r) run_chain(seed = 500 + r)))
mu <- colMeans(chains)
se <- apply(chains, 2, sd) / sqrt(nrep)
z <- (mu - ref.mean) / se

cat(sprintf(
    "  [3/3] geweke gate: E[phi1] = %.4f (ref %.4f, z = %.2f), E[phi2] = %.4f (ref %.4f, z = %.2f)\n",
    mu[1], ref.mean[1], z[1], mu[2], ref.mean[2], z[2]
))
if (any(abs(z) >= 4)) {
    stop(sprintf(
        paste0("GEWEKE FAIL: prior-recovery through updatePhiAR2TS deviates ",
               "from the truncated prior (z = %.2f, %.2f; gate |z| < 4). ",
               "The phi-block no longer targets the same joint as the ",
               "latent-blocks."), z[1], z[2]
    ))
}

cat("\nPASS: phi-block and latent-blocks share one stationary AR(2) joint density.\n")
