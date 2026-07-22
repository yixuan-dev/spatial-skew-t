# Posterior prediction on the Cincinnati grid for setting 204, the analogue of
# predict-all-71.R. Reads everything from the fit file written by
# us-all-full-204.R (fit, S.o, X.o, S.p, X.p) -- nothing is rebuilt here, so
# the MRTS basis cannot drift between fit and prediction.
#
# Differences from predict-all-71.R, all forced by setting 204:
#   * nknots = 1: fit$tau / fit$z are (draws x nt) matrices, fit$knots is
#     NULL, and the partition is trivial -- no mem() lookups.
#   * skew = TRUE with the z.init fix live: z and lambda come from the chain
#     instead of being hard-coded to zero. mu = X beta + lambda * z_g.
#
# Usage:  Rscript predict-cincy-204.R [dev]
options(warn = 2)
library(fields)
library(SpatialTools)
library(mvtnorm)

source('./package_load.R', chdir = TRUE)

args <- commandArgs(trailingOnly = TRUE)
suffix <- if (length(args) > 0 && tolower(args[1]) == "dev") "-dev" else ""
infile     <- paste0("results/us-all-full-204", suffix, ".RData")
outputfile <- paste0("us-all-pred-cincy-204", suffix, ".RData")
load(infile)

nreps <- dim(fit$tau)[1]
ns    <- nrow(S.o)
np    <- nrow(S.p)
nt    <- ncol(fit$tau)
stopifnot(dim(X.o)[3] == ncol(fit$beta), dim(X.p)[3] == ncol(fit$beta))
cat("draws:", nreps, " sites:", ns, " grid:", np, " days:", nt, "\n")

d11       <- rdist(S.p, S.p)
d22       <- rdist(S.o, S.o)
diag(d11) <- 0
diag(d22) <- 0
d12       <- rdist(S.p, S.o)
cov.model <- "matern"
x.beta    <- matrix(0, ns, nt)
taug      <- matrix(0, ns, nt)
zg        <- matrix(0, ns, nt)

# storage
y.pred <- array(0, c(nreps, np, nt))

set.seed(1)
for (i in 1:nreps) {
  # get values for current iteration from mcmc (nknots = 1 shapes)
  rho    <- fit$rho[i]
  nu     <- fit$nu[i]
  gamma  <- fit$gamma[i]
  beta   <- fit$beta[i, ]
  tau    <- matrix(fit$tau[i, ], nrow = 1)   # 1 x nt
  z      <- matrix(fit$z[i, ], nrow = 1)     # 1 x nt
  lambda <- fit$lambda[i]
  y      <- fit$y[i, , ]  # want to use imputed y not true y

  # precision matrix
  C <- gamma * simple.cov.sp(D=d22, sp.type=cov.model, sp.par=c(1, rho),
                             error.var=0, smoothness=nu, finescale.var=0)
  diag(C) <- 1
  prec <- chol2inv(chol(C))

  # update x.beta and taug; single partition, so knot 1 everywhere
  for (t in 1:nt) {
    x.beta[, t] <- X.o[, t, ] %*% beta
    taug[, t]   <- tau[1, t]
    zg[, t]     <- z[1, t]
  }

  mu  <- x.beta + lambda * zg
  res <- y - mu

  y.pred[i, , ] <- predictY_cont_lambda(d11 = d11, d12 = d12,
                                        cov.model = cov.model, rho = rho,
                                        nu = nu, gamma = gamma, res = res,
                                        beta = beta, tau = tau, taug = taug,
                                        z = z, prec = prec, lambda = lambda,
                                        s.pred = S.p, x.pred = X.p,
                                        knots = NULL)

  if (i %% 500 == 0) {
    print(paste("Iter", i))
  }
}
save(y.pred, S.p, cincy, file = outputfile)
cat("saved:", outputfile, "\n")
