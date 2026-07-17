#########################################################################
# Sanity checks for the non-t settings 13-21 in simdata_nonsta.RData.
#
# HISTORY: the extension was first generated as simdata_nonsta_ext.RData
# and verified against the then-frozen simdata_nonsta.RData on 2026-07-17
# -- settings 1-12 were IDENTICAL bit-for-bit (y, truemean.field, tau.t,
# z.t, knots.t, settings.nonsta, s, x, W.varying, W.nsdep all PASS) --
# and then merged back over simdata_nonsta.RData.  The bit-level
# comparison is therefore retired; what remains below are the invariants
# of the 13-21 extension itself, re-checkable at any time.
#
# READ-ONLY.
#########################################################################

if (basename(getwd()) != "simstudy") setwd("d:/Github/spatial-skew-t/code/analysis/simstudy")

new <- new.env()
load("simdata_nonsta.RData", envir = new)

ok <- function(cond, label) {
  cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", label))
  isTRUE(cond)
}
all_ok <- TRUE
chk <- function(cond, label) all_ok <<- ok(cond, label) && all_ok

cat("=== settings 13-21 sanity ===\n")
ns <- new$ns; nt <- new$nt; nsets <- new$nsets; s <- new$s
chk(dim(new$y)[4] == 21, "y carries 21 settings")

# every fixed surface: sd 11, orthogonal to [1, s1, s2]
X0 <- cbind(1, s)
for (nm in names(new$g.nont)) {
  g <- new$g.nont[[nm]]
  chk(abs(sd(g) - 11) < 1e-8, sprintf("g.nont$%s sd = 11", nm))
  chk(max(abs(crossprod(X0, g))) < 1e-6, sprintf("g.nont$%s orthogonal to [1,s1,s2]", nm))
}
chk(abs(sd(new$G.gpmean[, 1]) - 11) < 1e-8 &&
    max(abs(crossprod(X0, new$G.gpmean))) < 1e-6, "G.gpmean columns sd 11 + orthogonal")

# truemean.field contract
tm13 <- new$truemean.field[[13]]
chk(all(tm13 == 10 + new$g.nont$franke),               "truemean.field[[13]] = 10 + g.franke (all sets)")
chk(!identical(new$truemean.field[[17]][, 1],
               new$truemean.field[[17]][, 2]),          "truemean.field[[17]] varies per dataset")
chk(all(new$truemean.field[[17]] == 10 + new$G.gpmean), "truemean.field[[17]] = 10 + G.gpmean")
chk(identical(new$truemean.field[[13]], new$truemean.field[[18]]),
    "settings 13/18 share the same true mean")

# NA contract for tau/z/knots
chk(all(vapply(13:21, function(i) all(is.na(new$tau.t[[i]])), logical(1))),
    "tau.t[[13..21]] all NA")

# error scale and distribution.  Residual r = y - truemean has
# Var = sigma.g^2 (unit-diagonal C.stat); pooled over sites x time x 3 sets.
r13 <- sapply(1:3, function(k) new$y[, , k, 13] - tm13[, k])
chk(abs(sd(r13) - new$sigma.g) < 0.05,
    sprintf("setting 13 residual sd = %.3f ~ sigma.g = %g", sd(r13), new$sigma.g))
sk <- function(v) mean((v - mean(v))^3) / sd(v)^3
r18 <- sapply(1:3, function(k) new$y[, , k, 18] - tm13[, k])
chk(abs(sk(r13)) < 0.05, sprintf("setting 13 residual skewness ~ 0 (%.3f)", sk(r13)))
chk(sk(r18) > 0.5,       sprintf("setting 18 residual right-skewed (%.3f)", sk(r18)))

# the 13/18 A/B pairing: same underlying Z, only the transform differs.
# Invert the SAS transform on 18's standardised residual -> must equal 13's Z.
z13 <- r13 / new$sigma.g
z18 <- sinh((asinh(r18 / new$sigma.g * new$sas.mom$sd + new$sas.mom$mean) -
             new$sas.par$eps) / new$sas.par$delta)
chk(max(abs(z13 - z18)) < 1e-10, "13/18 share identical Gaussian innovations Z")

cat(sprintf("\n%s\n", if (all_ok) "ALL CHECKS PASSED" else "*** SOME CHECKS FAILED ***"))
if (!all_ok) quit(status = 1)
