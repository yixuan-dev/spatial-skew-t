# US-all-auto: Design Notes

---

## Primary Request and Intent

Build a model **auto-selection pipeline** for spatial ozone exceedance models.
The goal is to automatically choose the best candidate model configuration
(method, knots, threshold, covariates) without manual inspection of every
cross-validation table from the parent `US-all/` analysis.

The pipeline fits each candidate setting **once** on a fixed training split
and scores it on a held-out validation split using the Brier score at a target
exceedance probability. The setting with the lowest validation Brier is
selected and carried forward for final evaluation on the test split.

This replaces ad-hoc visual / table inspection and makes the selection
criterion explicit and reproducible.

---

## Key Technical Concepts

### 1. Single-fold fit — no cross-validation

Each candidate setting is fit **once** on `split.lst$train` (300 sites) and
evaluated on `split.lst$val` (100 sites). This is intentionally simpler than
the k-fold CV in `US-all/`: selection speed matters more than variance
reduction at this stage. The held-out `split.lst$test` (400 sites) is reserved
for final model evaluation after the best setting is chosen.

### 2. Site-partition structure

`split.lst` is the single source of truth for the data split. It contains
**integer row-indices** into the shared `Y`, `X`, `S` arrays:

```r
split.lst <- list(
  train = int[1:300],   # row indices into Y / X / S
  val   = int[1:100],
  test  = int[1:400]
)
```

Saved in `us-all-setup-auto.RData` together with `Y`, `X`, `S`,
`beta.init`, and `tau.init`. Reproducible via seed 2024 (default).

### 3. Posterior predictive layout — `fit$yp`

When `x.pred` / `s.pred` are passed to `mcmc()`, the returned object carries:

```
fit$yp   :  draws × n_val × ntime
```

`autoselect.R` reads this array directly; no re-fitting or post-processing is
needed to compute predictive exceedance probabilities.

### 4. Brier score at exceedance probability

The selection criterion is the **Brier score** at a target quantile `p`:

```
BS(p) = mean over (site j, time t) of
        ( P̂(Y_jt > u_p) − 1{Y_jt > u_p} )²
```

where `u_p = quantile(Y_all, p)` is computed from **all 800+ sites** (not
just the val split) to avoid split-induced threshold bias.
`P̂(Y_jt > u_p)` is the fraction of posterior draws exceeding `u_p`.

### 5. MCMC backends — both export `mcmc()`

Two backends are available, controlled by `US_ALL_MCMC_BACKEND` (default
`"ar2"`):

| Backend   | Files loaded by `package_load.R` | Main function |
| --------- | --------------------------------- | ------------- |
| `ar2`     | `ar2/mcmc_ar2.R` + helpers        | `mcmc()`      |
| `legacy`  | `mcmc_cont_lambda.R` + helpers    | `mcmc()`      |

Both backends expose **the same function name** `mcmc()`. `package_load.R`
loads only one set to prevent name conflicts. `run_settings_val.R` resolves the
function after loading:

```r
if (!exists("mcmc", envir = .GlobalEnv, inherits = TRUE))
  stop(...)
run_mcmc <- get("mcmc", envir = .GlobalEnv, inherits = TRUE)
```

Use `ar2` for all settings that include temporal structure (TS / AR2 lanes).

### 6. MRTS spatial covariates

Settings with `mrts > 0` append Markov-random-field thin-plate spline (MRTS)
basis columns to `X`. The basis is built via `autoFRK::mrts()` (loaded
unconditionally in `package_load.R` so a missing install fails at startup, not
mid-run). Columns with near-zero variance across training sites are dropped
before appending.

Each MRTS column is replicated across all time points, so the appended slice
of `X` has shape `[n_sites, ntime, k_kept]`.

### 7. Environment independence between settings

`run_one_setting()` is an ordinary R function. All data subsets (`y_tr`,
`x_tr`, `x_vl`, MRTS arrays) are local to each call and garbage-collected on
return. The RNG seed is set inside the function:

```r
set.seed(setting_id * 100 + 1L)
```

so each setting is deterministic and independent of run order.

### 8. Output contract per setting

Each `fits/val-<N>.RData` contains exactly two objects:

| Object         | Content                                              |
| -------------- | ---------------------------------------------------- |
| `fit`          | MCMC result; `fit$yp` = draws × n_val × ntime        |
| `runtime_info` | Named list: setting id, method, backend, timing, ... |

`autoselect.R` expects only `fit$yp` and the split from
`us-all-setup-auto.RData`; it does not depend on `runtime_info`.
