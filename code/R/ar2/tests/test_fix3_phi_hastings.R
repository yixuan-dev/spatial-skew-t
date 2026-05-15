# ============================================================================
# Regression test for the updatePhiAR2TS sampling scheme.
#
# History
# -------
# An earlier change to updatePhiAR2TS added a Hastings-style adjustment
#     log P(S_j; cur) - log P(S_j; can)
# to the acceptance ratio, on the (mistaken) ground that silent rejection of
# out-of-region candidates makes the proposal asymmetric in the Hastings
# sense. A detailed-balance proof in tex/where_we_are/where_we_are.tex
# (Prop. 3.3) shows that for a target pi supported on S and a symmetric
# untruncated proposal q_0, the Metropolis kernel
#
#     accept with prob  min(1, pi(can) / pi(cur))     if can in S,
#     stay otherwise                                   (silent rejection)
#
# already satisfies detailed balance for pi. NO Hastings term is needed.
# Adding one targets pi / P(S; .) instead, shifting mass toward the boundary
# of S (Prop. 3.4).
#
# The patch has been reverted. This test guards against re-introduction in
# two layers:
#
#   (1) STRUCTURAL: parse update_params_ar2.R and assert that the body of
#       updatePhiAR2TS contains no calls to pnorm(.) (the cdf calls that
#       would compute log P(S; .)). This is the strong regression assertion.
#
#   (2) SANITY: run a phi-only Markov chain with an interior truth and
#       verify the posterior median sits near the truth. This catches a
#       different family of regressions (e.g. likelihood or prior bugs in
#       updatePhiAR2TS) but is not the primary protection against the
#       Hastings re-introduction.
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

# ---------------------------------------------------------------------------
# (1) Structural assertion
# ---------------------------------------------------------------------------

extract_fn_lines <- function(src.lines, fn.name) {
    # Find the line "<fn.name> <- function(" and return the range of lines
    # spanning the function body via brace counting. We count "{" and "}"
    # robustly using nchar(gsub(...)); the more obvious gregexpr-length
    # idiom returns 1 for no-match (because gregexpr emits -1 as a length-1
    # vector to signal absence), which would collapse depth to zero at
    # line 1 and silently report a 1-line function body.
    start <- grep(paste0("^", fn.name, "\\s*<-\\s*function"), src.lines)
    if (length(start) == 0) {
        stop(sprintf("could not locate function %s", fn.name))
    }
    start <- start[1]
    depth <- 0
    end   <- NA_integer_
    seen.open <- FALSE
    count_char <- function(ch, s) nchar(gsub(paste0("[^", ch, "]"), "", s))

    for (i in seq(start, length(src.lines))) {
        ln <- sub("#.*$", "", src.lines[i])  # strip line comments
        opens  <- count_char("\\{", ln)
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

src <- readLines(src_path)
phi.body <- extract_fn_lines(src, "updatePhiAR2TS")

# Strip comments before testing for pnorm calls (we only want the live code)
phi.body.code <- sub("#.*$", "", phi.body)
pnorm.hits <- grep("\\bpnorm\\s*\\(", phi.body.code)

cat("=== Structural check (updatePhiAR2TS) ===\n")
cat(sprintf("Function body: %d lines\n", length(phi.body)))
cat(sprintf("pnorm(.) calls in live code: %d  (expected 0)\n",
            length(pnorm.hits)))
if (length(pnorm.hits) > 0) {
    cat("Offending lines:\n")
    for (i in pnorm.hits) {
        cat(sprintf("  %s\n", trimws(phi.body[i])))
    }
}

struct.ok <- (length(pnorm.hits) == 0)

# ---------------------------------------------------------------------------
# (2) Behavioural sanity check at interior truth
# ---------------------------------------------------------------------------

run_phi_chain <- function(updater, ystar, init, mh, n.iter, n.burn, prior.sd) {
    phi1 <- init[1]; phi2 <- init[2]
    acc1 <- att1 <- acc2 <- att2 <- 0
    k1 <- k2 <- numeric(n.iter)
    for (iter in seq_len(n.iter)) {
        r <- updater(
            data = ystar, phi1 = phi1, phi2 = phi2, day.mar = 2,
            att.phi1 = att1, acc.phi1 = acc1, mh.phi1 = mh,
            att.phi2 = att2, acc.phi2 = acc2, mh.phi2 = mh,
            prior.mean = 0, prior.sd = prior.sd
        )
        phi1 <- r$phi1; phi2 <- r$phi2
        k1[iter] <- phi1; k2[iter] <- phi2
        acc1 <- r$acc.phi1; att1 <- r$att.phi1
        acc2 <- r$acc.phi2; att2 <- r$att.phi2
    }
    keep <- (n.burn + 1):n.iter
    list(phi1 = k1[keep], phi2 = k2[keep])
}

phi1.t <- 0.30   # interior of the stationarity triangle
phi2.t <- 0.20

n.rep    <- 8
nt       <- 80
n        <- 4
mh       <- 0.08
n.iter   <- 5000
n.burn   <- 1000
prior.sd <- 0.5

cat("\n=== Sanity: posterior recovery at interior truth ===\n")
cat(sprintf("Truth (phi1, phi2) = (%.2f, %.2f); %d replicate datasets\n",
            phi1.t, phi2.t, n.rep))

bias.acc <- matrix(NA_real_, n.rep, 2)
for (r in seq_len(n.rep)) {
    set.seed(1000 + r)
    ystar <- simulate_ar2_standard(nt = nt, n = n,
                                    phi1 = phi1.t, phi2 = phi2.t)
    set.seed(2000 + r)
    out <- run_phi_chain(updatePhiAR2TS, ystar,
                         init = c(phi1.t, phi2.t),
                         mh = mh, n.iter = n.iter, n.burn = n.burn,
                         prior.sd = prior.sd)
    bias.acc[r, ] <- c(median(out$phi1) - phi1.t,
                       median(out$phi2) - phi2.t)
}

mean.bias <- colMeans(bias.acc)
cat(sprintf("Averaged median bias: phi1 = %+6.4f, phi2 = %+6.4f\n",
            mean.bias[1], mean.bias[2]))

tol.bias <- 0.05
sanity.ok <- max(abs(mean.bias)) < tol.bias
cat(sprintf("Sanity: |averaged bias|_inf < %.3f ? %s\n",
            tol.bias, if (sanity.ok) "yes" else "NO"))

# ---------------------------------------------------------------------------
# Outcome
# ---------------------------------------------------------------------------

cat("\n=== Outcome ===\n")
cat(sprintf("Structural (no pnorm in updatePhiAR2TS) : %s\n",
            if (struct.ok) "PASS" else "FAIL"))
cat(sprintf("Sanity     (interior recovery)          : %s\n",
            if (sanity.ok) "PASS" else "FAIL"))

if (struct.ok && sanity.ok) {
    cat("\nPASS: updatePhiAR2TS uses silent rejection and recovers interior truth.\n")
    quit(status = 0)
} else if (!struct.ok) {
    cat("\nFAIL (structural): pnorm() found inside updatePhiAR2TS. This is the\n")
    cat("       signature of the wrongly-shipped Hastings adjustment. Remove\n")
    cat("       the log P(S; .) terms; silent rejection alone is correct under\n")
    cat("       a symmetric untruncated proposal.\n")
    quit(status = 1)
} else {
    cat("\nFAIL (sanity): chain does not recover an interior truth. This is\n")
    cat("       not the Hastings-revert bug; look elsewhere in updatePhiAR2TS\n")
    cat("       (likelihood / prior / proposal mechanics).\n")
    quit(status = 1)
}
