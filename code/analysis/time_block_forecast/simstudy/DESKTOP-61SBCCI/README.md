# Experiment A (5 blocks), settings 5 and 7, on DESKTOP-61SBCCI

Tooling for the time-block forecast study's multi-block arm. Everything in
this folder is *tooling*; no data lives here except the driver transcript.

## What this run is

Settings **5** (near-unit-root AR(2), phi = (0.15, 0.80)) and **7**
(ARFIMA(0, 0.45, 0)), datasets **1:10**, methods **1** (i.i.d. in time),
**2** (AR(2) temporal) and **4** (AR(1) temporal), across all **5**
expanding-window blocks (seams 50/80/110/140/170, H = 15).

60 cells, 300 MCMC fits, `iters = 20000`, `burn = 10000`, ~1.8 h per cell,
~110 core-hours. Each cell is an 810 MB file, so the fits are **deleted
after scoring**; peak disk is one fit batch (~12 GB at `-ChunkSize 5`).

Before this run, Experiment A existed only for setting 4 at n = 5 datasets.

## Run it

```powershell
cd code\analysis\time_block_forecast\simstudy
.\DESKTOP-61SBCCI\expA_hn_driver.ps1 -DryRun      # plan + calls, runs nothing
.\DESKTOP-61SBCCI\expA_hn_driver.ps1              # workers auto-capped
```

`-Workers 0` (the default) caps at `min(CPU - 1, floor((RAM_GB - 6) / 2.0))`
and prints both caps; RAM is usually the binding one. Other knobs:
`-Settings`, `-Datasets`, `-Methods`, `-ChunkSize`, `-EsDraws`,
`-DiskFloorGB`, `-SkipEnvCheck`, `-NoCommit`.

Interrupted? Re-run the identical command. A dataset whose cache exists and
passes the gate is skipped entirely; valid fits are skipped by the
inventory; a fit truncated by a kill is detected by size and structure and
re-run.

## The lambda prior, and why there is no guard

Every fit runs with `--hn`, i.e. `lambda ~ HN(0, 20)` via the backend's
`lambda.positive` flag (commit `4ec3628`). This removes the discrete sign
reflection of the `(beta0, lambda, z)` ridge at the model level. The
block-1 counterfactual (`block1_positive_control/hn_prior_experiment/`)
showed 13/40 attempt-0 failures fall to 4/40 with scores statistically
indistinguishable from the guarded study.

**Assertion C is not used in this run.** There is no fit-level guard, no
reseed loop, no delete gate on C, and no cell exclusion rule: the primary
analysis uses all 50 (dataset, block) cells per setting. Since 2026-08-03
the A/A'/B/C flags are not recorded at all; `fit.diag` carries the 14
numeric summaries (`lambda`, `sdz_ratio`, `sd_lead_max`, the phis, ...)
per block, and `chunk_sanity_tbf.R` reports `lambda < 0` counts as a
descriptive ledger only.

Consequence: these caches are **not cell-comparable** to settings 1-4,
whose fits used the N(0, 20) prior and therefore a different RNG stream.
The `hn_prior` scalar inside every cache marks this, and the gate refuses
any cache where it is not `TRUE`.

## Pipeline per dataset

```
inventory_tbf.R  -> which (method, dataset) fits are missing, as run-settings.R calls
run-settings.R   -> results/<S>-<m>-<d>.RData  (810 MB, 5 blocks)
scores.R         -> output/results/scores<S>_d<d>.RData
chunk_sanity_tbf.R -> the gate; non-zero exit means the fits are NOT deleted
                    -> then results/<S>-<m>-<d>.RData is removed by exact name
merge_score_caches.R (reused from ../../simstudy/DESKTOP-61SBCCI/)
                 -> output/results/scores<S>.RData, verifying every input exactly
expA_threeway.R, tables.R, plots.R -> the tracked artifacts
```

## What survives the fit deletion

The fits are unrecoverable once deleted. Everything summarised while
`fit` was still in memory (the per-block `fit_diag` and `post_summary`)
rides inside the fit file and is carried into the score cache by
scores.R as the `fit.diag` and `post.summary` arrays -- the cache is the
durable copy. (Until 2026-08-03 run-settings.R also mirrored these, plus
a per-(site, lead) predictive summary, to per-cell files under
`output/diag/`; the mirrors were a second copy of the same rows and the
predictive summary had no consumer, so the mechanism was removed. Use
`-KeepFits` when a future re-scoring is anticipated.)

## What the gate checks before deleting anything

`chunk_sanity_tbf.R` refuses, among other things, a cache whose `--datasets`
or `--methods` were omitted (`scores.R` defaults to 1:50 and 1:2, and either
would delete unscored fits), a cache with any NA in the score arrays, a
cache where `hn_prior` is not `TRUE`, a cell that "finished" in under
30 minutes (`options(warn = 2)` makes any warning fatal inside `mcmc()`),
and NA anywhere in the `fit.diag` / `post.summary` arrays -- the durable
copies of the fit-time summaries.

## What returns in git

The 20 per-dataset caches, the 2 merged caches,
`output/tables/expA_*.csv`, `joint_summary{5,7}.csv`, the lead-curve
PDFs, this folder and the driver transcript. Under 5 MB. Fits never
enter git; `git add -f` is needed because the root `.gitignore` covers
`*.RData`.

Then, back on the machine that owns the report:

```powershell
git pull
Rscript ..\..\simstudy\DESKTOP-61SBCCI\merge_score_caches.R <tmp> output\results\scores5_d*.RData
# and diff <tmp> against the returned scores5.RData -- a free cross-machine check
```
