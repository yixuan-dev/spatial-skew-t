rm(list = ls())
source("./package_load.R", chdir = TRUE)
load("../ozone_data.RData")

setting <- 71
method <- "t"
nknots <- 10
keep.knots <- TRUE
threshold <- 75
tau.init <- 0.05
thresh.quant <- FALSE
skew <- FALSE
temporalw <- TRUE
temporalz <- TRUE
temporaltau <- TRUE
beta.init <- 0
tau.init <- 1
outputfile <- paste("results/us-all-full-", setting, ".RData", sep="")

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

# make design matrix
nt <- ncol(Y)
X  <- array(1, dim=c(nrow(cmaq), nt, 2))
for (t in 1:nt) {
  X[, t, 2] <- cmaq[, t]
}

run_started_at <- Sys.time()
start <- proc.time()

y.o <- Y
X.o <- X
S.o <- S

set.seed(setting * 100)
tic.set <- proc.time()
fit <- mcmc(y=y.o, s=S.o, x=X.o, # x.pred=X.p, s.pred=S.p,
            method=method, skew=skew, keep.knots=keep.knots,
            min.s=c(-2.25, -1.60), max.s=c(2.35, 1.30),
            thresh.all=threshold, thresh.quant=thresh.quant, nknots=nknots,
            iters=30000, burn=25000, update=500, iterplot=F,
            beta.init=beta.init, tau.init=tau.init,
            gamma.init=0.5, rho.init=1, rho.upper=5, nu.init=0.5, nu.upper=10,
            temporaltau=temporaltau, temporalw=temporalw, temporalz=temporalz)
toc.set <- proc.time()
time.set <- (toc.set - tic.set)[3]

elap.time.val <- (proc.time() - start)[3]
avg.time.val <- elap.time.val
runtime_finished_at <- run_started_at + as.numeric(elap.time.val)
runtime_info <- list(
  schema_version = 1L,
  runner = "us-all-full-71.R",
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
  started_at_utc = format(as.POSIXct(run_started_at, tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  finished_at_utc = format(as.POSIXct(runtime_finished_at, tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  elapsed_sec = unname(as.numeric(elap.time.val)),
  fit_elapsed_sec = unname(as.numeric(time.set)),
  control = list(iters = 30000L, burn = 25000L, update = 500L, thin = 1L)
)
save(fit, runtime_info, file=outputfile)
