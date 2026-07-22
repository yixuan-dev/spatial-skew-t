# Full-data fit for setting 204 (skew-t, nknots = 1, thresh = 0, CMAQ = yes,
# TS = yes, MRTS k = 10), the map-pipeline analogue of us-all-full-71.R.
#
# Differences from us-all-full-71.R:
#   * skew = TRUE with the temporal-z block -- requires the z.init fix in
#     mcmc_cont_lambda.R (see tex/z_init_bug); before that fix every legacy
#     skew-t + TS fit ran with z frozen at 0.
#   * MRTS k = 10 covariates appended to the design matrix (mrts_basis.R,
#     shared with us-all-run.R).
#   * The prediction target is the Cincinnati window, not the southeast: the
#     grid and its design matrix (CMAQ + MRTS pred basis) are built here and
#     saved with the fit so predict-cincy-204.R reuses them verbatim.
#
# US_ALL_RUN_MODE=dev gives a short chain and writes to a -dev suffixed file;
# the production output is results/us-all-full-204.RData.
rm(list = ls())
source("./package_load.R", chdir = TRUE)
source("./mrts_basis.R")
load("../ozone_data.RData")

RUN_MODE <- tolower(Sys.getenv("US_ALL_RUN_MODE", unset = "prod"))
if (!RUN_MODE %in% c("dev", "prod")) {
  stop("US_ALL_RUN_MODE must be one of: dev, prod", call. = FALSE)
}
mcmc_ctrl <- switch(RUN_MODE,
  dev  = list(iters = 2000,  burn = 1000,  update = 200),
  prod = list(iters = 30000, burn = 25000, update = 500)
)

setting      <- 204
method       <- "t"
skew         <- TRUE
nknots       <- 1
mrts.k       <- 10
keep.knots   <- FALSE   # nknots = 1: partition is trivial, knots never stored
threshold    <- 0
thresh.quant <- FALSE
temporalw    <- TRUE
temporalz    <- TRUE
temporaltau  <- TRUE
beta.init    <- 0
tau.init     <- 0.05    # runner's pick_tau_init for non-gaussian methods
outputfile   <- paste0("results/us-all-full-", setting,
                       if (RUN_MODE == "dev") "-dev" else "", ".RData")

# Cincinnati, OH (39.1031 N, 84.5120 W) in the CMAQ Lambert conformal system
# (+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=+b=6370000), rescaled
# by 1000 like the site coordinates. Window is +/- 200 km.
cincy.x  <- 1.0681
cincy.y  <- -0.0257
cincy.hw <- 0.20

# rescale x and y coordinates to make easier to work with
x <- x / 1000
y <- y / 1000

S       <- cbind(x[s[, 1]], y[s[, 2]])  # expands the grid of x, y
excl    <- which(rowMeans(is.na(Y)) > 0.50)  # remove where we're missing 50%
index   <- index[-excl]
Y       <- Y[-excl, ]
S       <- S[-excl, ]
CMAQ.cs <- (CMAQ - mean(CMAQ)) / sd(CMAQ)  # center and scale CMAQ data
cmaq    <- CMAQ.cs[index, ]  # extract cmaq for sites

#### Some site locations are duplicated which prevents us from having a positive
#### definite covariance matrix. So we find the duplicated sites and slightly
#### to make covariance matrix positive definite
set.seed(548837)  # jitter
d <- rdist(S)
same <- which(d == 0, arr.ind = TRUE)
same <- same[same[, 1] != same[, 2], ]
while (nrow(same) > 0) {
  S[same[1, 1], ] <- S[same[1, 1], ] + rnorm(1, 0, 0.00001)
  S[same[1, 2], ] <- S[same[1, 2], ] + rnorm(1, 0, 0.00001)
  d <- rdist(S)
  same <- which(d == 0, arr.ind = TRUE)
  same <- same[same[, 1] != same[, 2], ]
  print(nrow(same))
}

# prediction grid around Cincinnati
S.grid     <- expand.grid(x, y)
keep.these <- (S.grid[, 1] > cincy.x - cincy.hw & S.grid[, 1] < cincy.x + cincy.hw) &
              (S.grid[, 2] > cincy.y - cincy.hw & S.grid[, 2] < cincy.y + cincy.hw)
S.p    <- as.matrix(S.grid[keep.these, ])
cmaq.p <- CMAQ.cs[keep.these, ]
cat("Cincinnati grid:", nrow(S.p), "points\n")

# make design matrices: intercept + CMAQ + MRTS basis
nt <- ncol(Y)
X  <- array(1, dim = c(nrow(cmaq), nt, 2))
X.p <- array(1, dim = c(nrow(cmaq.p), nt, 2))
for (t in 1:nt) {
  X[, t, 2]   <- cmaq[, t]
  X.p[, t, 2] <- cmaq.p[, t]
}

mrts_cov <- build_mrts_covariates(S_train = S, S_pred = S.p, nt = nt, k = mrts.k)
cat("MRTS source:", mrts_cov$source,
    "| cols(original/kept/dropped)=", mrts_cov$original_cols, "/",
    mrts_cov$kept_cols, "/", mrts_cov$dropped_cols, "\n")

p_base <- dim(X)[3]
p_mrts <- dim(mrts_cov$train)[3]
X.ext   <- array(NA_real_, dim = c(dim(X)[1], nt, p_base + p_mrts))
X.p.ext <- array(NA_real_, dim = c(dim(X.p)[1], nt, p_base + p_mrts))
X.ext[, , seq_len(p_base)]   <- X
X.p.ext[, , seq_len(p_base)] <- X.p
for (j in seq_len(p_mrts)) {
  X.ext[, , p_base + j]   <- mrts_cov$train[, , j]
  X.p.ext[, , p_base + j] <- mrts_cov$pred[, , j]
}
X   <- X.ext
X.p <- X.p.ext
cat("total covariates:", dim(X)[3], "\n")

run_started_at <- Sys.time()
start <- proc.time()

y.o <- Y
X.o <- X
S.o <- S

set.seed(setting * 100)
tic.set <- proc.time()
fit <- mcmc(y=y.o, s=S.o, x=X.o,
            method=method, skew=skew, keep.knots=keep.knots,
            min.s=c(-2.25, -1.60), max.s=c(2.35, 1.30),
            thresh.all=threshold, thresh.quant=thresh.quant, nknots=nknots,
            iters=mcmc_ctrl$iters, burn=mcmc_ctrl$burn,
            update=mcmc_ctrl$update, iterplot=F,
            beta.init=beta.init, tau.init=tau.init,
            gamma.init=0.5, rho.init=1, rho.upper=5, nu.init=0.5, nu.upper=10,
            temporaltau=temporaltau, temporalw=temporalw, temporalz=temporalz)
toc.set <- proc.time()
time.set <- (toc.set - tic.set)[3]

elap.time.val <- (proc.time() - start)[3]
runtime_finished_at <- run_started_at + as.numeric(elap.time.val)
runtime_info <- list(
  schema_version = 1L,
  runner = "us-all-full-204.R",
  output_file = outputfile,
  setting = as.integer(setting),
  method = method,
  skew = isTRUE(skew),
  nknots = as.integer(nknots),
  keep_knots = isTRUE(keep.knots),
  threshold = as.numeric(threshold),
  thresh_quant = isTRUE(thresh.quant),
  temporalw = isTRUE(temporalw),
  temporalz = isTRUE(temporalz),
  temporaltau = isTRUE(temporaltau),
  mrts_k = as.integer(mrts.k),
  mrts_source = mrts_cov$source,
  mrts_kept_cols = as.integer(mrts_cov$kept_cols),
  z_init_fix = "mcmc_cont_lambda.R z.init NULL -> 0.6745/sqrt(tau.init); see tex/z_init_bug",
  started_at_utc = format(as.POSIXct(run_started_at, tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  finished_at_utc = format(as.POSIXct(runtime_finished_at, tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  elapsed_sec = unname(as.numeric(elap.time.val)),
  fit_elapsed_sec = unname(as.numeric(time.set)),
  control = list(iters = as.integer(mcmc_ctrl$iters), burn = as.integer(mcmc_ctrl$burn),
                 update = as.integer(mcmc_ctrl$update), thin = 1L)
)

# Everything predict-cincy-204.R needs: the jittered training sites, both
# finished design matrices (CMAQ + MRTS basis, observed and grid sides), and
# the grid itself. No rebuilding, so no basis drift between fit and predict.
mrts_meta <- list(k = mrts.k, source = mrts_cov$source,
                  keep_idx = mrts_cov$keep_idx,
                  kept_cols = mrts_cov$kept_cols)
cincy <- list(x = cincy.x, y = cincy.y, hw = cincy.hw, keep.these = which(keep.these))
save(fit, runtime_info, S.o, X.o, S.p, X.p, mrts_meta, cincy, file = outputfile)
cat("saved:", outputfile, "\n")
