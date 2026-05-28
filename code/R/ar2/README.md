# ar2 — AR(2) temporal MCMC backend

## Bug fixes

### `z.init` default changed from 0 to 1 (`mcmc_ar2.R`)

**Symptom:** When fitting methods with AR(2) temporal structure on `z`
(`ar2_z = TRUE`, e.g. methods 7 and 8 in the simulation study), the
posterior mean and SD of `phi1.z` and `phi2.z` were exactly 0 across
all MCMC iterations — the chain never moved.

**Root cause:** `z` was initialised to 0 by default.  The AR(2) z-updater
(`updateZTS_AR2`) works in the copula-transformed space
`z.star = hn.cop(z, sig) = qnorm(phn(z, sig))`.  At `z = 0`,
`phn(0, sig) = 0` so `z.star = qnorm(0) = -Inf`.  Every MH proposal
then drew `can.z.star = rnorm(1, -Inf, mh)` → `NaN`, making the
accept/reject ratio `R = NaN`, which is silently rejected by
`!is.na(R)`.  As a result:

1. `z` was permanently frozen at 0 (zg = 0 every iteration).
2. `updatePhiAR2TS` received `z.star = -Inf`, computed
   `cur.ll = -Inf` and `can.ll = NaN`, giving `R = NaN` →
   `phi1.z` and `phi2.z` stayed at their initial value of 0 forever.

Note: `phi1.tau` / `phi2.tau` were unaffected because `tau` is
initialised to 1 (positive), so `tau.star` is always finite.
`phi1.w` / `phi2.w` were also unaffected for `nknots > 1` because
`knots.star` is the probit transform of spatial coordinates (finite).

**Fix:** Changed the default in the `mcmc()` function signature from
`z.init = 0` to `z.init = NULL`.  When `NULL`, the initialisation block
computes the **median of the marginal half-normal**:

```r
z.init <- 0.6745 / sqrt(tau.init)   # 0.6745 = qnorm(0.75)
```

This is principled because `hn.cop(median(HN(sigma)), sigma) = qnorm(0.5) = 0`
exactly, so `z.star` starts at 0 — the centre of the copula space —
regardless of `tau.init`.  Callers who pass `z.init` explicitly are
unaffected.
