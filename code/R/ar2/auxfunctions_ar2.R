# Check if AR(2) parameters are stationary
check_ar2_stability <- function(phi1, phi2, tol = 1e-08) {
    cond1 <- abs(phi2) < (1 - tol)
    cond2 <- (phi1 + phi2) < (1 - tol)
    cond3 <- (phi2 - phi1) < (1 - tol)
    cond1 && cond2 && cond3
}

get_ar2_standard_moments <- function(phi1, phi2, tol = 1e-06, param_name = NULL) {
    if (!check_ar2_stability(phi1, phi2, tol = tol / 10)) {
        if (!is.null(param_name)) {
            stop(paste0(
                "AR(2) parameters on ", param_name,
                " are not stationary; please adjust phi1.", param_name,
                " and phi2.", param_name
            ))
        } else {
            stop("AR(2) parameters are not stationary; please adjust phi1 and phi2")
        }
    }

    denom <- 1 - phi2
    if (abs(denom) < tol) {
        if (!is.null(param_name)) {
            stop(paste0("phi2.", param_name, " is too close to 1; AR(2) cannot remain stable"))
        } else {
            stop("phi2 is too close to 1; AR(2) cannot remain stable")
        }
    }

    gamma0 <- 1
    gamma1 <- phi1 * gamma0 / denom
    gamma2 <- phi1 * gamma1 + phi2 * gamma0
    innov_var <- gamma0 - phi1 * gamma1 - phi2 * gamma2
    if (innov_var <= tol) {
        if (!is.null(param_name)) {
            stop(paste0(
                "AR(2) innovation variance on ", param_name,
                " is non-positive; please check phi1.", param_name,
                " and phi2.", param_name
            ))
        } else {
            stop("AR(2) innovation variance is non-positive; please check phi1 and phi2")
        }
    }

    list(
        gamma0 = gamma0,
        gamma1 = gamma1,
        gamma2 = gamma2,
        innov_var = innov_var,
        innov_sd = sqrt(innov_var)
    )
}

get_ar2_conditional_params <- function(t, lag1 = NULL, lag2 = NULL,
                                       phi1, phi2, moments = NULL) {
    if (is.null(moments)) {
        moments <- get_ar2_standard_moments(phi1, phi2)
    }

    if (t <= 0) {
        stop("AR(2) time index must be positive")
    }

    if (t == 1) {
        return(list(mean = 0, sd = sqrt(moments$gamma0)))
    }

    if (t == 2) {
        if (is.null(lag1)) {
            stop("lag1 is required for AR(2) t = 2 conditional density")
        }
        # Fix 1: stationary AR(2) conditional from the bivariate joint
        #   (Y_1, Y_2) ~ N_2(0, [[1, gamma_1], [gamma_1, 1]]) ,
        # which yields  Y_2 | Y_1 = lag1  ~  N(gamma_1 * lag1, 1 - gamma_1^2).
        # The previous form N(phi_1 * lag1, sigma_eps^2) was the AR(2)
        # recursion extrapolated backward with a fictitious Y_0 = 0; the
        # two coincide only when phi_2 = 0. Using the stationary form
        # preserves Var(Y_t) = 1 for all t and matches the joint draw
        # in simulate_ar2_standard / draw_ar2_initial.
        g1    <- moments$gamma1
        sd.t2 <- sqrt(moments$gamma0 - g1^2)
        return(list(mean = g1 * lag1, sd = sd.t2))
    }

    if (is.null(lag1) || is.null(lag2)) {
        stop("lag1 and lag2 are required for AR(2) t > 2 conditional density")
    }

    list(mean = phi1 * lag1 + phi2 * lag2, sd = moments$innov_sd)
}

ar2_conditional_log_density <- function(x, t, lag1 = NULL, lag2 = NULL,
                                        phi1, phi2, moments = NULL) {
    params <- get_ar2_conditional_params(
        t = t, lag1 = lag1, lag2 = lag2,
        phi1 = phi1, phi2 = phi2, moments = moments
    )

    dnorm(x, params$mean, params$sd, log = TRUE)
}

ar2_transition_loglik <- function(current, lag1, lag2, phi1, phi2, moments = NULL) {
    if (is.null(moments)) {
        moments <- get_ar2_standard_moments(phi1, phi2)
    }

    sum(dnorm(current, phi1 * lag1 + phi2 * lag2, moments$innov_sd, log = TRUE))
}

# Function to simulate stationary AR(2) time series
draw_ar2_initial <- function(n, gamma0, gamma1) {
    Sigma <- matrix(c(gamma0, gamma1, gamma1, gamma0), nrow = 2, ncol = 2) # covariance matrix
    chol.Sigma <- tryCatch(
        chol(Sigma),
        error = function(e) {
            eig.inv(Sigma, inv = FALSE, logdet = FALSE)$sd.mtx
        }
    ) # Cholesky decomposition

    z <- matrix(rnorm(2 * n), nrow = 2, ncol = n) # 2 x n matrix
    draws <- t(chol.Sigma %*% z) # n x 2 matrix
    return(draws)
}

simulate_ar2_standard <- function(nt, n, phi1, phi2, tol = 1e-06, param_name = NULL) {
    moments <- get_ar2_standard_moments(phi1, phi2, tol = tol, param_name = param_name)
    gamma0 <- moments$gamma0
    gamma1 <- moments$gamma1
    innov.sd <- moments$innov_sd

    ar2 <- matrix(0, nrow = n, ncol = nt)

    if (nt == 1) {
        ar2[, 1] <- rnorm(n, 0, sqrt(gamma0))
        return(ar2)
    }

    init.states <- draw_ar2_initial(n, gamma0, gamma1)
    ar2[, 1] <- init.states[, 1]
    ar2[, 2] <- init.states[, 2]

    if (nt >= 3) {
        for (t in 3:nt) {
            ar2[, t] <- phi1 * ar2[, t - 1] + phi2 * ar2[, t - 2] + rnorm(n, 0, innov.sd)
        }
    }

    return(ar2)
}

makeTauAR2 <- function(nt, nknots, tau.alpha, tau.beta, phi1, phi2) {
    tau.star <- simulate_ar2_standard(
        nt = nt, n = nknots, phi1 = phi1, phi2 = phi2,
        param_name = "tau"
    )
    tau <- gamma.invcop(x = tau.star, tau.alpha, tau.beta)
    return(tau)
}

makeZAR2 <- function(nt, nknots, tau, phi1, phi2) {
    z.star <- simulate_ar2_standard(
        nt = nt, n = nknots, phi1 = phi1, phi2 = phi2,
        param_name = "z"
    )
    sd <- 1 / sqrt(tau)
    z <- hn.invcop(x = z.star, sig = sd)
    return(z)
}

makeKnotsAR2 <- function(nt, nknots, s, phi1, phi2) {
    knots.star <- knots <- array(NA, dim = c(nknots, 2, nt))

    knots.star[, 1, ] <- simulate_ar2_standard(
        nt = nt, n = nknots, phi1 = phi1, phi2 = phi2,
        param_name = "w"
    )
    knots.star[, 2, ] <- simulate_ar2_standard(
        nt = nt, n = nknots, phi1 = phi1, phi2 = phi2,
        param_name = "w"
    )

    min.s1 <- min(s[, 1])
    max.s1 <- max(s[, 1])
    min.s2 <- min(s[, 2])
    max.s2 <- max(s[, 2])

    knots[, 1, ] <- transform$inv.probit(knots.star[, 1, ], lower = min.s1, upper = max.s1)
    knots[, 2, ] <- transform$inv.probit(knots.star[, 2, ], lower = min.s2, upper = max.s2)

    return(knots)
}

# ---- Generalized AR(p) simulation ----

# Check if AR(p) parameters define a stationary process.
# Uses the companion matrix: all eigenvalues must lie strictly inside the unit circle.
check_arp_stability <- function(phi, tol = 1e-8) {
    p <- length(phi)
    if (p == 0) return(TRUE)
    if (p == 1) return(abs(phi[1]) < (1 - tol))

    F_mat <- matrix(0, p, p)
    F_mat[1, ] <- phi
    F_mat[2:p, 1:(p - 1)] <- diag(p - 1)

    max(Mod(eigen(F_mat, only.values = TRUE)$values)) < (1 - tol)
}

# Compute the stationary covariance of the AR(p) initial state vector and
# the corresponding innovation variance, normalised so that gamma(0) = 1.
#
# Method: solve the discrete Lyapunov equation
#   Sigma = F Sigma F' + Q   (unit innovation variance Q[1,1] = 1)
# via direct linear solve on the vectorised form, then rescale.
#
# Returns: list(Sigma        = p×p Toeplitz initial-state covariance,
#               sigma2_innov = innovation variance for gamma(0)=1)
get_arp_stationary_params <- function(phi, tol = 1e-8) {
    p <- length(phi)

    F_mat <- matrix(0, p, p)
    F_mat[1, ] <- phi
    if (p > 1) F_mat[2:p, 1:(p - 1)] <- diag(p - 1)

    Q <- matrix(0, p, p)
    Q[1, 1] <- 1  # unit innovation; will normalise below

    # (I - F ⊗ F) vec(Sigma) = vec(Q)
    A_lyap <- diag(p * p) - kronecker(F_mat, F_mat)
    Sigma_raw <- matrix(solve(A_lyap, as.vector(Q)), p, p)

    gamma0_raw <- Sigma_raw[1, 1]
    if (gamma0_raw <= tol) {
        stop("AR(p) stationary variance is non-positive; check phi coefficients")
    }

    Sigma <- Sigma_raw / gamma0_raw          # normalised: Sigma[1,1] = 1
    sigma2_innov <- 1 / gamma0_raw           # innovation variance for gamma(0) = 1

    if (sigma2_innov <= tol) {
        stop("AR(p) innovation variance is non-positive; check phi coefficients")
    }

    list(Sigma = Sigma, sigma2_innov = sigma2_innov)
}

# Simulate n independent stationary AR(p) paths of length nt.
# Returns an n × nt matrix with marginal variance = 1.
simulate_arp_standard <- function(nt, n, phi, tol = 1e-6, param_name = NULL) {
    phi <- as.numeric(phi)
    p <- length(phi)

    if (!check_arp_stability(phi, tol = tol / 10)) {
        nm <- if (!is.null(param_name)) paste0(".", param_name) else ""
        stop(paste0(
            "AR(", p, ") parameters on phi", nm,
            " are not stationary; please check the phi coefficients"
        ))
    }

    out <- matrix(0.0, nrow = n, ncol = nt)
    if (nt == 0) return(out)

    # AR(1): closed-form innovation SD avoids the Lyapunov solve
    if (p == 1) {
        sigma_innov <- sqrt(1 - phi[1]^2)
        out[, 1] <- rnorm(n)
        for (t in seq_len(nt - 1) + 1) {
            out[, t] <- phi[1] * out[, t - 1] + sigma_innov * rnorm(n)
        }
        return(out)
    }

    # AR(p), p >= 2: draw initial p time steps jointly from stationary distribution
    sp <- get_arp_stationary_params(phi, tol = tol)
    sigma_innov <- sqrt(sp$sigma2_innov)

    chol_Sigma <- tryCatch(
        chol(sp$Sigma),
        error = function(e) eig.inv(sp$Sigma, inv = FALSE, logdet = FALSE)$sd.mtx
    )
    # Draw n realisations of the p-vector; result is n × p
    init_draws <- t(chol_Sigma %*% matrix(rnorm(p * n), p, n))
    out[, 1:min(p, nt)] <- init_draws[, 1:min(p, nt), drop = FALSE]

    if (nt > p) {
        for (t in (p + 1):nt) {
            ar_mean <- rep(0.0, n)
            for (j in seq_len(p)) ar_mean <- ar_mean + phi[j] * out[, t - j]
            out[, t] <- ar_mean + sigma_innov * rnorm(n)
        }
    }

    return(out)
}

# Exact stationary Gaussian ARFIMA(0, d, 0) / fractional-Gaussian-noise sampler.
# Sibling of simulate_arp_standard: returns an n x nt matrix whose one-point
# marginal is EXACTLY N(0, 1) (variance 1, Gaussian), so the copula transforms
# downstream (gamma.invcop / hn.invcop / inv.probit) stay valid. Unlike the
# finite-order AR path, this carries genuine long memory (hyperbolic ACF,
# sum_h |rho(h)| = Inf), the misspecified DGP that AR(1)/AR(2) analysis models
# cannot represent.
#
# Method: Davies-Harte / Wood-Chan exact circulant embedding. The target
# autocorrelation obeys rho(h)/rho(h-1) = (h-1+d)/(h-d) with rho(0)=1 (Hosking
# 1981); since gamma(0)=1 the autocovariance equals the autocorrelation.
simulate_arfima_standard <- function(nt, n, d, tol = 1e-8, param_name = NULL) {
    if (!is.finite(d) || d <= -0.5 || d >= 0.5) {
        stop(sprintf("ARFIMA(0,d,0) requires d in (-0.5, 0.5); got %s", format(d)))
    }
    N <- as.integer(nt)
    out <- matrix(0.0, nrow = n, ncol = max(N, 0L))
    if (N <= 0L) return(out)
    if (N == 1L) { out[] <- rnorm(n); return(out) }

    # Autocovariances gamma(0..N-1) (== autocorrelations, since gamma(0)=1),
    # built by ratio recursion to avoid gamma-function overflow.
    g <- numeric(N)
    g[1L] <- 1
    for (h in seq_len(N - 1L)) g[h + 1L] <- g[h] * ((h - 1L + d) / (h - d))

    # Minimal circulant embedding, size m = 2(N-1); first column is the
    # symmetric extension (gamma_0..gamma_{N-1}, gamma_{N-2}..gamma_1).
    m <- 2L * (N - 1L)
    cseq <- c(g, if (N >= 3L) g[(N - 1L):2L] else numeric(0))
    lambda <- Re(stats::fft(cseq))                 # circulant eigenvalues (real, >= 0)
    neg_tol <- tol * max(abs(lambda))
    if (any(lambda < -neg_tol)) {
        nm <- if (!is.null(param_name)) paste0(" for '", param_name, "'") else ""
        stop(sprintf(
            "Davies-Harte embedding%s has negative eigenvalues (min %.3e, d=%.4f); not PSD-embeddable at nt=%d.",
            nm, min(lambda), d, N))
    }
    lambda[lambda < 0] <- 0
    sl <- sqrt(lambda)
    half <- m %/% 2L                               # = N - 1

    for (i in seq_len(n)) {
        # Hermitian-symmetric complex weights with E|Z_k|^2 = 1, so that
        # Re(fft(sqrt(lambda) * Z)) / sqrt(m) is real with covariance gamma.
        Z <- complex(length.out = m)
        Z[1L] <- rnorm(1L)                         # k = 0 (real)
        Z[half + 1L] <- rnorm(1L)                  # k = m/2 (real)
        if (half >= 2L) {
            kR <- 2L:half                          # 0-based k = 1 .. half-1
            Z[kR] <- complex(real = rnorm(half - 1L), imaginary = rnorm(half - 1L)) / sqrt(2)
            Z[m - kR + 2L] <- Conj(Z[kR])          # enforce Hermitian symmetry
        }
        y <- Re(stats::fft(sl * Z)) / sqrt(m)
        out[i, ] <- y[1L:N]
    }
    out
}

# Dispatch a latent-process spec to the matching unit-variance N(0,1) sampler:
# finite-order AR(p) via simulate_arp_standard, long memory via
# simulate_arfima_standard. Both return an n x nt matrix with marginal N(0,1).
draw_latent_standard <- function(spec, nt, n, param_name = NULL) {
    if (!is.null(spec$type) && identical(spec$type, "arfima")) {
        return(simulate_arfima_standard(nt = nt, n = n, d = spec$d, param_name = param_name))
    }
    simulate_arp_standard(nt = nt, n = n, phi = spec$coeffs, param_name = param_name)
}

# Yule-Walker pseudo-true coefficients: the AR(1) and AR(2) coefficients an
# analysis model converges to (in the M-closest projection sense) when the true
# process is ARFIMA(0, d, 0). Uses the first two theoretical autocorrelations,
# rho(1) = d/(1-d) and rho(2) = d(1+d)/((1-d)(2-d)). The AR(2) projection has a
# POSITIVE phi2 = (rho2 - rho1^2)/(1 - rho1^2) for all d in (0, 0.5) -- the
# signature that AR(1) is misspecified and AR(2) captures the lag-2 partial
# autocorrelation. Returns d, Hurst exponent, rho1/rho2, and both projections.
arfima_yw_projection <- function(d) {
    if (!is.finite(d) || d <= -0.5 || d >= 0.5) {
        stop(sprintf("arfima_yw_projection: d must be in (-0.5, 0.5); got %s", format(d)))
    }
    rho1 <- d / (1 - d)
    rho2 <- d * (1 + d) / ((1 - d) * (2 - d))
    ar2_phi1 <- rho1 * (1 - rho2) / (1 - rho1^2)
    ar2_phi2 <- (rho2 - rho1^2) / (1 - rho1^2)
    list(
        d = d, hurst = d + 0.5,
        rho1 = rho1, rho2 = rho2,
        ar1 = rho1,
        ar2 = c(phi1 = ar2_phi1, phi2 = ar2_phi2)
    )
}

# Theoretical ARFIMA(0, d, 0) autocorrelation rho(0..H) via the ratio recursion
# rho(h)/rho(h-1) = (h-1+d)/(h-d), rho(0)=1. Used by the implied-ACF diagnostic
# to overlay the true long-memory ACF against the AR(1)/AR(2) geometric fits.
arfima_acf <- function(d, H) {
    r <- numeric(H + 1L)
    r[1L] <- 1
    for (h in seq_len(H)) r[h + 1L] <- r[h] * ((h - 1L + d) / (h - d))
    r
}

makeTauARP <- function(nt, nknots, tau.alpha, tau.beta, spec) {
    tau.star <- draw_latent_standard(spec, nt = nt, n = nknots, param_name = "tau")
    tau <- gamma.invcop(x = tau.star, tau.alpha, tau.beta)
    return(tau)
}

makeZARP <- function(nt, nknots, tau, spec) {
    z.star <- draw_latent_standard(spec, nt = nt, n = nknots, param_name = "z")
    sd <- 1 / sqrt(tau)
    z <- hn.invcop(x = z.star, sig = sd)
    return(z)
}

makeKnotsARP <- function(nt, nknots, s, spec) {
    knots.star <- knots <- array(NA, dim = c(nknots, 2, nt))

    knots.star[, 1, ] <- draw_latent_standard(spec, nt = nt, n = nknots, param_name = "w")
    knots.star[, 2, ] <- draw_latent_standard(spec, nt = nt, n = nknots, param_name = "w")

    min.s1 <- min(s[, 1]); max.s1 <- max(s[, 1])
    min.s2 <- min(s[, 2]); max.s2 <- max(s[, 2])

    knots[, 1, ] <- transform$inv.probit(knots.star[, 1, ], lower = min.s1, upper = max.s1)
    knots[, 2, ] <- transform$inv.probit(knots.star[, 2, ], lower = min.s2, upper = max.s2)

    return(knots)
}

parse_phi_spec <- function(phi, param_name) {
    # Long-memory spec: list(type = "arfima", d = <d>). Intercepted before the
    # numeric-vector path so the fractional differencing parameter is not
    # coerced through unlist(). order = Inf makes the order >= 1 dispatch in
    # rpotspatTS_arp select the make*ARP path, where draw_latent_standard then
    # routes on $type.
    if (is.list(phi) && identical(phi$type, "arfima")) {
        d <- phi$d
        if (is.null(d) || !is.numeric(d) || length(d) != 1 ||
            !is.finite(d) || d <= -0.5 || d >= 0.5) {
            stop(sprintf(
                "phi.%s: arfima spec requires a single numeric d in (-0.5, 0.5)",
                param_name))
        }
        return(list(order = Inf, type = "arfima", d = as.numeric(d), coeffs = numeric(0)))
    }

    coeffs <- if (is.list(phi)) {
        unlist(phi, recursive = TRUE, use.names = FALSE)
    } else {
        phi
    }

    if (length(coeffs) == 0 || is.null(coeffs)) {
        return(list(order = 0, coeffs = numeric(0)))
    }

    if (!is.numeric(coeffs)) {
        stop(paste0(
            "phi.", param_name,
            " must be numeric. Preferred input is a numeric vector c(phi1, ..., phip) for AR(p); list is kept for backward compatibility."
        ))
    }

    coeffs <- as.numeric(coeffs)

    if (any(!is.finite(coeffs))) {
        stop(paste0("phi.", param_name, " contains non-finite values"))
    }

    if (length(coeffs) == 1 && coeffs[1] == 0) {
        return(list(order = 0, coeffs = coeffs))
    }

    list(order = length(coeffs), coeffs = coeffs)
}

# Geometric-anisotropy ("deformed") exponential covariance.
#   C(s,s') = gamma * exp(-||A (s-s')|| / rho),  diag = 1
# where A = R(theta) %*% diag(c(1, ratio)) deforms 2D coordinates so that
# the kernel becomes isotropic exponential in the deformed space.
#
# Args:
#   s     : ns x 2 matrix of coordinates
#   gamma : marginal correlation scale (off-diagonal multiplier)
#   rho   : range parameter (deformed-space scale)
#   theta : rotation angle in radians for the principal axis
#   ratio : aspect ratio for the second axis (1 = isotropic)
CorFxDef <- function(s, gamma, rho, theta = 0, ratio = 1) {
    ns <- nrow(s)
    if (rho < 1e-6) {
        return(diag(1, nrow = ns))
    }
    R <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
    A <- R %*% diag(c(1, ratio))
    s.def <- s %*% t(A)
    d.def <- as.matrix(stats::dist(s.def))
    cor <- gamma * exp(-d.def / rho)
    diag(cor) <- 1
    return(cor)
}

ensure_rpotspatTS_arp_dependencies <- function() {
    required_symbols <- c(
        "CorFx", "mem", "makeKnotsTS", "gamma.invcop", "hn.invcop", "transform", "g.Rcpp"
    )

    missing_symbols <- required_symbols[!vapply(
        required_symbols,
        exists,
        logical(1),
        mode = "any",
        inherits = TRUE
    )]

    if (length(missing_symbols) > 0) {
        load_all_candidates <- c(
            "code/00_core/load_all.R",
            "00_core/load_all.R",
            "load_all.R"
        )
        load_all_file <- load_all_candidates[file.exists(load_all_candidates)][1]

        if (!is.na(load_all_file) && length(load_all_file) == 1) {
            source(load_all_file)
        }
    }

    missing_symbols <- required_symbols[!vapply(
        required_symbols,
        exists,
        logical(1),
        mode = "any",
        inherits = TRUE
    )]

    if (length(missing_symbols) > 0) {
        stop(paste0(
            "rpotspatTS_arp() is missing required objects: ",
            paste(missing_symbols, collapse = ", "),
            ". Please run source('code/00_core/load_all.R') before calling this function."
        ))
    }

    invisible(TRUE)
}

rpotspatTS_arp <- function(nt, x, s, beta,
                           gamma, nu, rho,
                           lambda, tau.alpha, tau.beta,
                           nknots, dist,
                           phi.z = 0, phi.w = 0, phi.tau = 0,
                           cov.type = c("matern", "deformed"),
                           theta = 0, ratio = 1) {
    ensure_rpotspatTS_arp_dependencies()
    cov.type <- match.arg(cov.type)

    p <- dim(x)[3]
    ns <- nrow(s)

    y <- matrix(NA, ns, nt)
    tau <- matrix(NA, nknots, nt)
    z <- matrix(NA, nknots, nt)
    g <- matrix(NA, ns, nt)
    tau.alpha <- tau.alpha / 2 # reparameterizaiton
    tau.beta <- tau.beta / 2 # reparameterization

    d <- as.matrix(dist(s))

    if (lambda == 0) {
        skew <- FALSE
    } else {
        skew <- TRUE
    }

    C <- if (cov.type == "deformed") {
        CorFxDef(s = s, gamma = gamma, rho = rho, theta = theta, ratio = ratio)
    } else {
        CorFx(d = d, gamma = gamma, rho = rho, nu = nu)
    }
    chol.C <- chol(C)
    t.chol.C <- t(chol.C)

    phi.tau.spec <- parse_phi_spec(phi.tau, "tau")
    phi.z.spec <- parse_phi_spec(phi.z, "z")
    phi.w.spec <- parse_phi_spec(phi.w, "w")

    if (dist == "t") {
        if (phi.tau.spec$order >= 1) {
            tau <- makeTauARP(
                nt = nt, nknots = nknots,
                tau.alpha = tau.alpha, tau.beta = tau.beta,
                spec = phi.tau.spec
            )
        } else {
            tau <- matrix(rgamma(nknots * nt, tau.alpha, tau.beta), nknots, nt)
        }
    } else if (dist == "gaussian") {
        tau <- matrix(0.25, nknots, nt)
    }
    sd <- 1 / sqrt(tau)

    if (skew) {
        if (phi.z.spec$order >= 1) {
            z <- makeZARP(nt = nt, nknots = nknots, tau = tau, spec = phi.z.spec)
        } else {
            z <- abs(matrix(rnorm(nknots * nt, 0, sd), nknots, nt))
        }
    } else {
        z <- matrix(0, nrow = nknots, ncol = nt)
    }

    if (phi.w.spec$order >= 1) {
        knots <- makeKnotsARP(nt = nt, nknots = nknots, s = s, spec = phi.w.spec)
    } else {
        knots <- makeKnotsTS(nt = nt, nknots = nknots, s = s, phi = 0)
    }

    for (t in 1:nt) {
        knots.t <- matrix(knots[, , t], nknots, 2)
        g <- mem(s, knots.t)
        zg.t <- z[g, t]
        taug.t <- sqrt(tau[g, t])
        sdg <- 1 / taug.t

        if (p == 1) {
            x.beta <- matrix(x[, t, ], ns, 1) * beta
        } else {
            x.beta <- x[, t, ] %*% beta
        }

        mu <- x.beta + lambda * zg.t
        y.t <- mu + t.chol.C %*% matrix(rnorm(ns, 0, sdg), ns, 1)
        y[, t] <- y.t
    }

    results <- list(y = y, tau = tau, z = z, knots = knots)
    return(results)
}


################################################################
# Arguments:
#   preds(iters, np, nt): mcmc predictions at validation
#                         locations
#   probs(nprobs): sample quantiles for scoring
#   validate(np, nt): validation data
#
# Returns:
#   score(nprobs): a single quantile score per quantile
################################################################
Quant_Score <- function(preds, probs, validate) {
    # preds: (iters × np × nt) MCMC 樣本
    # probs: c(0.9, 0.95, 0.98, 0.99) 要評估的分位數
    # validate: (np × nt) 真實驗證資料

    np <- nrow(validate) # number of prediction sites
    nt <- ncol(validate) # number of prediction days
    nprobs <- length(probs) # number of quantiles to find quantile score

    # we get the predicted quantile over all sites and times
    pred.quants <- quantile(preds, probs = probs, na.rm = T)

    scores.sites <- array(NA, dim = c(nprobs, np, nt)) # QS per site per time

    for (q in 1:nprobs) {
        diff <- pred.quants[q] - validate # np x nt matrix of differences
        i <- diff >= 0 # diff >= 0 means qhat is larger
        scores.sites[q, , ] <- 2 * (i - probs[q]) * diff
    }

    scores <- apply(scores.sites, 1, mean, na.rm = T) # average over times
    return(scores)
}

################################################################
# Arguments:
#   preds(iters, np, nt): mcmc predictions at validation
#                         locations
#   probs(nprobs): sample quantiles for scoring
#   validate(np, nt): validation data
#   trans(bool): are the mcmc predictions transposed
#
# Returns:
#   score(nprobs): a single quantile score per quantile per site
################################################################
Quant_Score_Site <- function(preds, probs, validate, trans = FALSE) {
    nt <- ncol(validate) # number of prediction days
    np <- nrow(validate) # number of prediction sites
    nprobs <- length(probs) # number of quantiles to find quantile score

    # we get the predicted quantile for each site nprobs x np
    # for each site, we estimate the quantiles over all times
    if (trans) {
        pred.quants <- apply(preds, 3, quantile, probs = probs, na.rm = T)
    } else {
        pred.quants <- apply(preds, 2, quantile, probs = probs, na.rm = T)
    }

    pred.quants <- t(pred.quants) # need np x nprobs for proper matrix subtraction
    scores.sites <- array(NA, dim = c(nprobs, np, nt)) # QS per site per time

    # we need to figure out how many times the site did or didn't exceed the prediction
    for (q in 1:nprobs) {
        diff <- pred.quants[, q] - validate # np x nt matrix of differences
        i <- diff >= 0 # diff >= 0 means qhat is larger
        scores.sites[q, , ] <- 2 * (i - probs[q]) * diff
    }

    # for each site, average over times
    scores <- apply(scores.sites, c(1, 2), mean, na.rm = T) # nprobs x np
    scores <- t(scores) # np x nprobs
    return(scores) # np x nprobs
}

################################################################
# Arguments:
#   preds(iters, num_missing_points): mcmc predictions at
#                         validation locations
#   probs(nprobs): sample quantiles for scoring
#   validate(numeric): validation data
#
# Returns:
#   score(nprobs): a single quantile score per quantile
################################################################
# Quant_Score_general <- function(preds, probs, validate) {
#     num_missing_points <- length(validate) # number of prediction points
#     nprobs <- length(probs) # number of quantiles to find quantile score

#     # we get the predicted quantile over all sites and times
#     pred.quants <- quantile(preds, probs = probs, na.rm = T)

#     scores <- array(NA, dim = c(nprobs, num_missing_points))

#     for (q in 1:nprobs) {
#         diff <- pred.quants[q] - validate # num_missing_points vector of differences
#         i <- diff >= 0 # diff >= 0 means qhat is larger
#         scores[q, ] <- 2 * (i - probs[q]) * diff
#     }

#     scores <- apply(scores, 1, mean, na.rm = T) # average over times
#     return(scores)
# }


################################################################
# Update date: 2026-01-13
# Author: Yi-Xuan Xie
# Description: Generalized quantile score function that works
#              for predictions in any shape (e.g., iters x np x nt,
#              iters x length(preds), etc.)
# Arguments:
#   preds(iters, np, nt) or (iters, length(preds)): mcmc predictions at
#                         validation locations
#   probs(nprobs): sample quantiles for scoring
#   validate(np, nt) or (length(validate)): validation data
#
# Returns:
#   score(nprobs): a single quantile score per quantile
################################################################
quant_score_general <- function(preds, probs, validate) {
    dim_val <- if (is.null(dim(validate))) length(validate) else dim(validate)
    nprobs <- length(probs) # number of quantiles to find quantile score

    # we get the predicted quantile over all sites and times
    pred.quants <- quantile(preds, probs = probs, na.rm = T)

    scores <- array(NA, dim = c(nprobs, dim_val))

    scores <- sapply(1:nprobs, function(q) {
        diff <- pred.quants[q] - validate
        i <- diff >= 0
        2 * (i - probs[q]) * diff
    }, simplify = "array")

    scores <- aperm(scores, c(length(dim(scores)), seq_len(length(dim(scores)) - 1)))

    scores <- apply(scores, 1, mean, na.rm = T) # average over times
    return(scores)
}


################################################################
# Arguments:
#   preds(iters, yp, nt): mcmc predictions at validation
#                         locations
#   thresholds(nthreshs): sample thresholds for scoring
#   validate(np, nt): validation data
#
# Returns:
#   score(nthreshs): a single Brier score per threshold
################################################################
# Brier_Score_general <- function(preds, thresholds, validate) {
#     nthreshs <- length(thresholds)
#     scores <- rep(NA, nthreshs)

#     for (b in 1:nthreshs) {
#         pat <- apply((preds > thresholds[b]), 2, mean) # np x nt

#         i <- validate > thresholds[b] # np x nt
#         scores[b] <- mean((i - pat)^2, na.rm = T)
#     }

#     return(scores)
# }

################################################################
# Update date: 2026-01-13
# Author: Yi-Xuan Xie
# Description: Generalized Brier score function that works
#              for predictions in any shape (e.g., iters x np x nt,
#              iters x length(preds), etc.)
# Arguments:
#   preds(iters, np, nt) or (iters, length(preds)): mcmc predictions at
#                         validation locations
#   thresholds(nthreshs): sample thresholds for scoring
#   validate(np, nt) or (length(validate)): validation data
#
# Returns:
#   score(nthreshs): a single Brier score per threshold
################################################################
brier_score_general <- function(preds, thresholds, validate) {
    dims <- dim(preds)
    # Legacy reference:
    # temp <- seq_len(length(dim(preds)))[-1]
    # scores <- sapply(1:nthreshs, function(b) {
    #     pat <- apply((preds > thresholds[b]), temp, mean)
    #     i <- validate > thresholds[b]
    #     mean((i - pat)^2, na.rm = TRUE)
    # })
    pred_mat <- matrix(preds, nrow = dims[1], ncol = prod(dims[-1]))
    validate_vec <- as.vector(validate)

    scores <- vapply(seq_along(thresholds), function(b) {
        pat_vec <- colMeans(pred_mat > thresholds[b])
        i <- validate_vec > thresholds[b]
        mean((i - pat_vec)^2, na.rm = TRUE)
    }, numeric(1))

    return(scores)
}
