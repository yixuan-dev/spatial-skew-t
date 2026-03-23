# US-all ozone experiments: full-done reproduction checklist

## Scope
This checklist reproduces the models used in the supplemental "Best performing models" table,
with the primary success criterion set to **ranking and conclusion consistency** (not bitwise-equal decimals).

## Phase A (completed)
- [x] `package_load.R` portability fix (optional private source + optional BLAS threading function)
- [x] `us-all-setup.R` Windows-safe `quartz()` fallback
- [x] `us-all-setup.RData` generated/updated

## Runtime assumptions
- Windows Rscript path used in this repo: `C:\Program Files\R\R-4.5.1\bin\Rscript.exe`
- Working directory: `d:\Github\spatial-skew-t\code\analysis\ozone\US-all`
- Input data file: `..\ozone_data.RData`

## Reproduction target (full done settings)

Run the same settings used by `us-all-results.R`:

```r
done <- c(1:5, 7:9, 11:13, 15:17, 33:36, 38:41, 43:46, 51:74)
```

Total settings: **52**.

## Batch partition (full run, recommended)

### Batch F1 (non-TS base block)
- settings: `1:5, 7:9, 11:13, 15:17`

### Batch F2 (non-TS extended knots)
- settings: `33:36, 38:41, 43:46`

### Batch F3 (TS block 1)
- settings: `51:56`

### Batch F4 (TS block 2)
- settings: `57:62`

### Batch F5 (TS block 3)
- settings: `63:68`

### Batch F6 (TS block 4)
- settings: `69:74`

## Post-batch checks
After each batch, verify files exist:
- `results/us-all-<setting>.RData`

### Expected checks by batch
- F1: `results/us-all-1.RData` ... `results/us-all-17.RData` (only listed F1 settings)
- F2: `results/us-all-33.RData` ... `results/us-all-46.RData` (only listed F2 settings)
- F3: `results/us-all-51.RData` ... `results/us-all-56.RData`
- F4: `results/us-all-57.RData` ... `results/us-all-62.RData`
- F5: `results/us-all-63.RData` ... `results/us-all-68.RData`
- F6: `results/us-all-69.RData` ... `results/us-all-74.RData`

## Final aggregation (full run)
After all F1-F6 settings finish, run:
- `us-all-results.R`

Outputs:
- `us-all-results.RData` (from `savelist`)
- Console output for per-threshold best settings (from `which(min(...))` blocks)
- Derived objects: `brier.score.mean`, `bs.mean.ref.gau`, and top-2 selection in `score.compare`

## Acceptance criteria (primary)

The reproduction is accepted when:
- Top-1 / Top-2 ordering for
	`q(0.90), q(0.95), q(0.98), q(0.99), q(0.995)`
	matches the paper table in **ranking and conclusion**.
- Mapped `(TS, K, T)` from `settings.csv` is consistent with the supplemental table rows.

Numerical values may differ slightly due to platform/R/BLAS/RNG differences.

## Notes
- `us-all-results-quick.R` can still be used as a preflight sanity check, but it is not the final criterion.
- Prefer checkpointing after each batch and rerun only missing settings.

---

## AR2 smoke phase (comparison workflow gate)

This phase validates AR2 wiring and output contract before any large AR2 run.

### Scope
- AR2 smoke settings: `101, 102, 103`
- Loader: `ar2_load.R` (keeps `package_load.R` unchanged)
- Expected outputs:
	- `results/us-all-101.RData`
	- `results/us-all-102.RData`
	- `results/us-all-103.RData`

### Run options
- By batch alias:
	- `Rscript repro_batch/batch_auto_check.R --batch=AR2_SMOKE --label=AR2-SMOKE --stop-on-batch-failure`
- Equivalent explicit settings:
	- `Rscript repro_batch/batch_auto_check.R --settings=101,102,103 --label=AR2-SMOKE --stop-on-batch-failure`

### Acceptance checks (smoke gate)
- `AR2-SMOKE-run-status-latest.csv` reports all settings with `run_status = ok`.
- Each `results/us-all-<setting>.RData` exists for `101:103`.
- Each fit object has prediction samples available for scoring:
	- `fit[[1]][['yp']]` and `fit[[2]][['yp']]` are present and non-empty.

### Next step after passing smoke
- Proceed to build/extend AR2 settings for full comparison and aggregate with a dedicated results script (without altering baseline `us-all-results.R`).
