# time_block_forecast / simstudy

Simulation study for the **time-block holdout** redesign of
`tex/time_block_strategy/ar2_rethink.tex`. It tests whether AR(2)
temporal pooling improves prediction of extremes once the evaluation is
no longer structurally blind to temporal dependence.

## Why this study exists

The Morris baseline study (`code/analysis/simstudy/`) scored AR(2) and
i.i.d. skew-t models with the Brier score under a **spatial-only**
holdout and found a null result. `ar2_rethink.tex` proves this is an
artifact of the experiment, not of temporal pooling:

- **Proposition 1** — the Brier score depends on the predictive only
  through marginal exceedance probabilities.
- **Proposition 2** — under a spatial-only holdout the AR(2) and i.i.d.
  priors induce the *same* time-`t` marginal, so a marginal score
  cannot separate them; the lower-dimensional model wins the
  bias–variance contest.
- **Section 4** — the cure is to hold out **time blocks** and forecast
  forward from the posterior of the seam state, scoring by lead time.

This study implements that redesign.

## Pipeline (mirrors `code/analysis/simstudy/`)

1. **Data generation** — `setup.R` → `simdata.RData`
2. **Fit + forecast** — `run-settings.R` → `results/<setting>-<method>-<dataset>.RData`
3. **Post-fit pipeline**
   - `scores.R` (Stage 1) → `output/results/scores<setting>.RData`
   - `tables.R` (Stage 2) → CSV tables + `simresults<setting>.RData`
   - `plots.R`  (Stage 3) → PDF lead-time curves

`time_block_helpers.R` carries the shared CLI / filename / seed / catalog
helpers and the forecasting function. The AR(2) model code is imported
through the Morris study's loader, `../../simstudy/ar2_load.R` — this
study keeps no private copy of the backend.

## Scripts

| script                 | role                                                                   |
| ---------------------- | ---------------------------------------------------------------------- |
| `setup.R`              | generate `simdata.RData` (**implemented + run**)                       |
| `time_block_helpers.R` | CLI / filename / seed / catalogs / block geometry / `forecast_block()` |
| `run-settings.R`       | per-block prefix fit + Algorithm 1 forecast driver                     |
| `scores.R`             | Stage 1 — lead-time CRPS/Brier + energy + variogram scores             |
| `tables.R`             | Stage 2 — lead-time curve tables, relative curve, crossing lead        |
| `plots.R`              | Stage 3 — lead-time curve PDFs (SE band, h*, crossing-lead marker)     |

## Data settings (`--setting`, the data-generating axis)

Every setting is the same family — skew-t, `dist = "t"`, `K = 1` knot,
`lambda = 3` — and differs only in the AR(2) coefficient pair `phi`
shared by the latent processes `tau*, z*, w*` (eq. 2 of the note):

| setting | label    | `phi = (phi1, phi2)` | role                                               |
| ------- | -------- | -------------------- | -------------------------------------------------- |
| 1       | iid      | (0.00, 0.00)         | null control — AR(2) must **not** beat i.i.d. here |
| 2       | weak     | (0.12, -0.05)        | small spectral radius, short memory horizon        |
| 3       | moderate | (0.60, -0.30)        | the note's reference case (Remark 5)               |
| 4       | strong   | (0.80, -0.35)        | long memory horizon                                |

With `K = 1` the Voronoi membership is constant, so `w` carries no
signal and the temporal effect lives entirely in `tau` and `z` — this
isolates the AR(2) signal exactly as Remark 8 recommends.

## Method catalog (`--methods`, the analysis axis)

| method | label                                                                                  |
| ------ | -------------------------------------------------------------------------------------- |
| 1      | Skew-t, K=1, i.i.d. in time (baseline)                                                 |
| 2      | Skew-t, K=1, AR(2) temporal (tau, z, w)                                                |
| 3      | Skew-t, K=1, AR(2) temporal, fixed membership (Remark 8 ablation; ≡ method 2 when K=1) |

## Time-block holdout geometry

Stored in `simdata.RData`:

- `nt = 200` — long record so several near-independent seams fit.
- `block_H = 15` — forecast horizon (lead times `1..H`), within the
  practical range `H ∈ [5,15]` of Definition 7.
- `block_seams = c(50, 80, 110, 140, 170)` — five expanding-window
  blocks. Block `b` fits on the contiguous prefix `y[, 1:T_o]` and
  forecasts the window `(T_o, T_o + H]`. Seams are spread across the
  record so the five seam states are near-independent.

## Forecasting (`forecast_block()`, Algorithm 1 of the note)

For each MCMC draw `m`:

1. extract the seam state `(X_{T_o-1}, X_{T_o})` from draw `m`'s imputed
   latent trajectory — **never** from the stationary distribution
   (Corollary 1: that erases the AR(2) signal);
2. recurse the AR(2) on the latent Gaussian scale with Yule-Walker
   innovation variance (eq. 5), so the marginal variance stays 1;
3. apply the copula transforms (eqs. 8–9) **after** the recursion;
4. draw the spatial field.

Averaging the predictive over draws propagates parameter *and*
seam-state uncertainty (Proposition 3). The i.i.d. baseline instead
draws each latent slice from the `N(0,1)` marginal, reproducing the
stationary predictive of Corollary 1.

## Scoring (Section 4.3)

- **Lead-time curve** `S_bar(h)` — univariate proper scores (CRPS, and
  Brier at a grid of high thresholds) evaluated separately at each lead
  `h`, with a standard error from the dispersion of the `B` per-block
  means (Definition 9). Expected signature of a genuine effect: the
  AR(2) curve lies below the i.i.d. baseline at short leads and
  converges to it near the memory horizon `h*` — the **crossing lead**
  is the headline quantity.
- **Joint summary** — energy score and variogram score over the
  per-site length-`H` time vector, one number per method.

## Usage

```powershell
$R = "C:\Program Files\R\R-4.5.1\bin\Rscript.exe"

# Stage 0: generate data (already done)
& $R .\setup.R

# Stage A: fit + forecast (smoke test: setting 3, dataset 1, both methods)
& $R .\run-settings.R --setting=3 "1" 1 "1:2"

# Stage 1 + 2 + 3: scores + tables + plots
& $R .\scores.R --setting=3 --datasets="1" --methods="1:2"
& $R .\tables.R --setting=3
& $R .\plots.R  --setting=3
```

## Status

- `setup.R` — implemented and executed; `simdata.RData` generated
  (144 sites × 200 times × 50 datasets × 4 settings).
- `run-settings.R`, `scores.R`, `tables.R`, `plots.R`,
  `time_block_helpers.R` — implemented to the design above; not yet
  validated by a run. Run the
  smoke test before scaling out, and check the three diagnostics of
  Remark 9 (`H=0` parity, `phi=0` collapse, long-horizon `N(0,1)`).
