# ozone_prop — ozone analysis with the prop backend

Mirrors `code/analysis/ozone/US-all` but swaps the legacy `mcmc()` backend for
the prop backend in `code/R/prop/`. Same data, same 2-fold CV split.

## Files

| File                       | Role                                                                  |
| -------------------------- | --------------------------------------------------------------------- |
| `ozone_data.RData`         | Raw ozone inputs (CMAQ, Y, grids, borders).                           |
| `ozone-prop-setup.RData`   | Preprocessed `Y`, `X`, `S`, `cv.lst`, `beta.init`, `tau.init`.        |
| `prop_load.R`              | Loads setup data + sources `code/R/prop/*.R`.                         |
| `ozone-prop-run.R`         | Runner. Reads `settings.csv`, runs prop `mcmc()` per fold.            |
| `settings.csv`             | **(awaiting user)** experiment grid.                                  |
| `results/`                 | Per-setting outputs (created on first run).                           |

## Important: `mrts` column → `prop_k`

The reference `us-all-run.R` uses the `mrts` column to build extra
covariate columns (autoFRK basis) appended to `X`. The **prop backend
does not take MRTS covariates**; instead the prop model has its own
basis controlled by the `prop_k` argument of `mcmc()`.

In this directory the `mrts` column of `settings.csv` is reinterpreted:

- value (positive integer) → passed as `prop_k = <value>` to `mcmc()`
- empty / non-positive → **error** (every prop setting must specify a rank)

## Settings columns honored

| Column    | Meaning                                                                |
| --------- | ---------------------------------------------------------------------- |
| `setting` | ID; numeric prefix used as seed base                                   |
| `method`  | `gaussian` / `t` / `skew-t` (no `max-stable` — prop unsupported)       |
| `knots`   | `nknots` (gaussian must be 1)                                          |
| `thresh`  | `thresh.all` (numeric, `thresh.quant=FALSE`)                           |
| `CMAQ`    | `yes` use full `X`; `no` keep intercept only                           |
| `TS`      | must be `no` — prop backend has no temporal blocks                     |
| `ar2`     | must be empty/no — prop backend has no AR2 blocks                      |
| `mrts`    | **prop_k** (positive integer, required)                                |

## Environment variables

| Variable                       | Default     | Meaning                                  |
| ------------------------------ | ----------- | ---------------------------------------- |
| `OZONE_PROP_RUN_MODE`          | `prod`      | `dev` → 2000/1000/200; `prod` → 30000/25000/500 |
| `OZONE_PROP_RESULTS_DIR`       | `results`   | output directory                         |
| `OZONE_PROP_SETTINGS`          | (unset)     | settings spec, e.g. `1:5` or `1,3,7`     |
| `OZONE_PROP_SETTINGS_FILE`     | `settings.csv` | path to settings file                 |
| `OZONE_PROP_COV_UPDATE_EVERY`  | `1`         | `prop_cov_update_every` for `mcmc()`     |

## Quick start (PowerShell)

```powershell
Set-Location "D:\Github\spatial-skew-t\code\analysis\ozone_prop"
$cand = Get-ChildItem 'C:\Program Files\R' -Directory | Sort-Object Name -Descending | Select-Object -First 1
$rs = Join-Path $cand.FullName 'bin\Rscript.exe'

# smoke test
$env:OZONE_PROP_RUN_MODE = "dev"
& $rs ozone-prop-run.R 1

# full prod run on a range
Remove-Item Env:\OZONE_PROP_RUN_MODE -ErrorAction SilentlyContinue
& $rs ozone-prop-run.R 1:10
```

Output: `results/ozone-prop-<setting>.RData` containing `fit` (length 2,
one per CV fold) and `runtime_info`.
