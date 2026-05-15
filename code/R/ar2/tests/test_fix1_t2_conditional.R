# ============================================================================
# Unit test for Fix 1: the t == 2 branch of get_ar2_conditional_params must
# return the *stationary* AR(2) conditional density
#
#     Y_2 | Y_1 = y_1  ~  N( gamma_1 * y_1,  1 - gamma_1^2 ),     gamma_1 = phi_1 / (1 - phi_2)
#
# which is what bivariate-Gaussian conditioning gives under the stationary
# joint (Y_1, Y_2) ~ N_2(0, [[1, gamma_1], [gamma_1, 1]]).
#
# The original (pre-Fix 1) implementation returned
#
#     N( phi_1 * y_1,  sigma_eps^2 ),
#
# which is the AR(2) recursion form extrapolated backward with a fictitious
# Y_0 = 0. It coincides with the stationary form only when phi_2 = 0.
#
# This test compares the function's returned (mean, sd) at t = 2 to the
# closed-form stationary values for several (phi_1, phi_2) configurations
# spanning the interior and the near-boundary regions of the stationarity
# triangle. The y_1 value is fixed and known so we can construct the
# expected target deterministically.
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

# ---- Configurations: each row is (phi1, phi2, y1)
configs <- rbind(
    c( 0.7, -0.4,  1.2),
    c( 0.4,  0.3,  0.8),
    c(-0.3, -0.2, -0.5),
    c( 0.8, -0.6,  0.3),
    c( 0.0,  0.5,  1.0),    # phi1 = 0 -> stationary mean is 0
    c( 0.3,  0.0,  2.0)     # phi2 = 0 -> stationary and recursion coincide
)
colnames(configs) <- c("phi1", "phi2", "y1")

tol <- 1e-10   # closed-form comparison; tight tolerance is appropriate

cat(sprintf("%-6s | %-6s | %-6s | %-10s %-10s | %-10s %-10s | %s\n",
            "phi1", "phi2", "y1",
            "expect mu", "got mu",
            "expect sd", "got sd",
            "match"))
cat(strrep("-", 95), "\n", sep = "")

all_ok <- TRUE

for (i in seq_len(nrow(configs))) {
    phi1 <- configs[i, "phi1"]
    phi2 <- configs[i, "phi2"]
    y1   <- configs[i, "y1"]

    moments <- get_ar2_standard_moments(phi1, phi2)
    gamma1  <- moments$gamma1                  # phi1 / (1 - phi2)

    # Stationary target for (Y_2 | Y_1 = y1):
    expect_mu <- gamma1 * y1
    expect_sd <- sqrt(1 - gamma1^2)

    got <- get_ar2_conditional_params(t = 2, lag1 = y1,
                                       phi1 = phi1, phi2 = phi2,
                                       moments = moments)

    ok_mu <- abs(got$mean - expect_mu) < tol
    ok_sd <- abs(got$sd   - expect_sd) < tol
    ok    <- ok_mu && ok_sd

    cat(sprintf("%+6.2f | %+6.2f | %+6.2f | %+10.6f %+10.6f | %10.6f %10.6f | %s\n",
                phi1, phi2, y1,
                expect_mu, got$mean,
                expect_sd, got$sd,
                if (ok) "yes" else "NO"))
    if (!ok) all_ok <- FALSE
}

if (all_ok) {
    cat("\nPASS: get_ar2_conditional_params returns the stationary form at t = 2.\n")
    quit(status = 0)
} else {
    cat("\nFAIL: the t == 2 branch does not match the stationary conditional.\n")
    cat("      Expected mean  : gamma_1 * y_1            = phi_1 / (1 - phi_2) * y_1\n")
    cat("      Expected sd    : sqrt(1 - gamma_1^2)\n")
    cat("      Fix             : update the t == 2 branch in\n")
    cat("                       code/R/ar2/auxfunctions_ar2.R\n")
    quit(status = 1)
}
