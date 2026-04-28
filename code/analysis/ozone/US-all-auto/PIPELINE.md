# US-all-auto Pipeline

Auto-selection experiment for spatial ozone exceedance models. Fits each
candidate model **once** on a training split and selects the best by Brier
score on a held-out validation split. No cross-validation.

---

## Pipeline Overview

```
us-all-setup-auto.R   →  us-all-setup-auto-<tr>-<va>-<te>.RData  (data + split)
settings-auto.R       →  settings-auto.csv                        (candidate grid)
run_settings_val.R    →  fits/val-<N>.RData                       (one fit/setting)
autoselect.R          →  output/…/tables + .RData                 (best by val Brier)
```

### Step 1 — Build the data split

```powershell
Rscript us-all-setup-auto.R
```

Loads `us-all-setup.RData` (copied from `../US-all/`), partitions the 800
CV-covered sites into train / validation / test, and saves
`us-all-setup-auto-<N_TRAIN>-<N_VAL>-<N_TEST>.RData` containing:

| Object                  | Description                           |
| ----------------------- | ------------------------------------- |
| `Y`, `X`, `S`           | Full dataset (all 800+ sites)         |
| `beta.init`, `tau.init` | MCMC initialisation values            |
| `split.lst$train`       | 300 site row-indices (into Y / X / S) |
| `split.lst$val`         | 100 site row-indices                  |
| `split.lst$test`        | 400 site row-indices                  |

Split is controlled by environment variables (see table below). Default seed
is 2024, so the split is reproducible.

| Env var                     | Default | Effect                        |
| --------------------------- | ------- | ----------------------------- |
| `US_ALL_AUTOSELECT_N_TRAIN` | `300`   | Number of training sites      |
| `US_ALL_AUTOSELECT_N_VAL`   | `100`   | Number of validation sites    |
| `US_ALL_AUTOSELECT_N_TEST`  | `400`   | Number of test sites          |
| `US_ALL_AUTOSELECT_SEED`    | `2024`  | Random seed for the partition |

### Step 2 — Generate the settings grid

```powershell
Rscript settings-auto.R
```

Reads `../US-all/settings.csv` and `../US-all/output/us-all/tables/comparison_full_table.csv`,
scores each setting by composite tail Brier, and writes `settings-auto.csv`
with tier / lane annotations.

| Tier | N   | Description                                         |
| ---- | --- | --------------------------------------------------- |
| 0    | 1   | Gaussian reference (setting 1)                      |
| 1    | 41  | Core candidates — top evidence, recommended default |
| 2    | 38  | Extended coverage — TS / AR2 / MRTS variants        |
| 3    | 36  | Non-TS baseline settings, not recommended           |

### Step 3 — Fit candidate models on training split

```powershell
Rscript run_settings_val.R
```

For each requested setting: fits MCMC on `split.lst$train`, generates
posterior predictive draws at `split.lst$val`, saves `fits/val-<N>.RData`
with:

- `fit$yp` — draws × n_val × ntime posterior predictive matrix
- `runtime_info` — metadata (timing, model spec, run mode)

### Step 4 — Select best model

```powershell
Rscript autoselect.R
```

Loads each `fits/val-<N>.RData`, computes Brier score at target exceedance
probability, picks the setting with the lowest validation Brier. Writes:

- `output/us-all-auto/tables/autoselect_validation_scores.csv`
- `output/us-all-auto/results/us-all-auto-results.RData`

---

## run_settings_val.R — Environment Variables

### Required

| Env var                      | Format                              | Example                |
| ---------------------------- | ----------------------------------- | ---------------------- |
| `US_ALL_AUTOSELECT_SETTINGS` | Comma-separated IDs or `a:b` ranges | `'8,12,15:16,204:206'` |

Integers and range tokens can be mixed freely. The script processes them in
ascending order and skips any setting not found in `settings-auto.csv`.

### Optional

| Env var                  | Default                               | Values          | Effect                                                                                                                                             |
| ------------------------ | ------------------------------------- | --------------- | -------------------------------------------------------------------------------------------------------------------------------------------------- |
| `US_ALL_MCMC_BACKEND`    | `ar2`                                 | `legacy`, `ar2` | Selects the MCMC implementation. `ar2` adds AR(2) temporal priors; required for TS / AR2 lane settings.                                            |
| `US_ALL_VAL_RUN_MODE`    | `prod`                                | `prod`, `dev`   | `prod`: 30 000 iters / 25 000 burn. `dev`: 2 000 / 1 000 for quick sanity checks.                                                                  |
| `US_ALL_VAL_RESULTS_DIR` | `fits/`                               | any path        | Output directory for `val-<N>.RData` files. Created automatically if it does not exist.                                                            |
| `US_ALL_VAL_SETUP_FILE`  | `us-all-setup-auto-200-200-400.RData` | any `.RData`    | Setup file from `us-all-setup-auto.R`. Default filename encodes the train/val/test split (200/200/400). Accepts absolute or script-relative paths. |

---

## Common Run Recipes

All commands assume:
```powershell
cd D:\Github\spatial-skew-t\code\analysis\ozone\US-all-auto
$env:US_ALL_MCMC_BACKEND = 'ar2'
```

### Validate pipeline (dev mode, 2 settings)

```powershell
$env:US_ALL_VAL_RUN_MODE        = 'dev'
$env:US_ALL_AUTOSELECT_SETTINGS = '111,206'
Rscript run_settings_val.R
```

### Run Tier 1 — production (41 settings)

```powershell
$env:US_ALL_VAL_RUN_MODE        = 'prod'
$env:US_ALL_AUTOSELECT_SETTINGS = '3,8,12:13,15:16,29,34:35,38:39,41:42,51,54:55,58,60:62,67:68,70,74,101,105,108,111:112,114,116:118,120,124,204:206,214:216'
Rscript run_settings_val.R
```

### Run Tier 2 (38 settings)

```powershell
$env:US_ALL_AUTOSELECT_SETTINGS = '52:53,56:57,59,63:66,69,71:73,102:104,106:107,109:110,113,115,119,121:123,201:203,207:209,211:213,217:219'
Rscript run_settings_val.R
```

### Run Tier 1 + Tier 2 together (79 settings)

```powershell
$env:US_ALL_AUTOSELECT_SETTINGS = '3,8,12:13,15:16,29,34:35,38:39,41:42,51:74,101:124,201:209,211:219'
Rscript run_settings_val.R
```

### Run specific settings

```powershell
$env:US_ALL_AUTOSELECT_SETTINGS = '111,204,205,206'
Rscript run_settings_val.R
```

### Re-run a single failed setting

```powershell
$env:US_ALL_AUTOSELECT_SETTINGS = '58'
Rscript run_settings_val.R
```

### Use a custom split and fits directory

Set both scripts to the same setup file and fits directory.

```powershell
$env:US_ALL_VAL_SETUP_FILE          = 'us-all-setup-auto-300-100-400.RData'
$env:US_ALL_VAL_RESULTS_DIR         = 'fits-300-100-400'
$env:US_ALL_AUTOSELECT_SETTINGS     = '111,204:206'
Rscript run_settings_val.R

$env:US_ALL_AUTOSELECT_SETUP_FILE   = 'us-all-setup-auto-300-100-400.RData'
$env:US_ALL_VAL_RESULTS_DIR         = 'fits-300-100-400'
$env:US_ALL_AUTOSELECT_SETTINGS     = '111,204:206'
Rscript autoselect.R
```

---

## Tier 1 Setting IDs (quick reference)

```
3, 8, 12, 13, 15, 16, 29, 34, 35, 38, 39, 41, 42,
51, 54, 55, 58, 60, 61, 62, 67, 68, 70, 74,
101, 105, 108, 111, 112, 114, 116, 117, 118, 120, 124,
204, 205, 206, 214, 215, 216
```

Range-token form: `3,8,12:13,15:16,29,34:35,38:39,41:42,51,54:55,58,60:62,67:68,70,74,101,105,108,111:112,114,116:118,120,124,204:206,214:216`

> Note: 214–216 (MRTS k=10 + AR2) have `has_result = FALSE` in settings-auto.csv — results pending.

---

## autoselect.R — Environment Variables

| Env var                          | Default                     | Effect                                                                                              |
| -------------------------------- | --------------------------- | --------------------------------------------------------------------------------------------------- |
| `US_ALL_AUTOSELECT_SETTINGS`     | — (required)                | Same format as run_settings_val.R; must match the settings that have been fitted.                   |
| `US_ALL_AUTOSELECT_TARGET_PROBS` | `0.90,0.95,0.98,0.99,0.995` | Comma-separated exceedance probabilities.                                                           |
| `US_ALL_AUTOSELECT_OBJECTIVE`    | `each`                      | `each` = best setting per target prob; `mean` = best by mean Brier across all probs.               |
| `US_ALL_SUMMARY_DRAWS`           | `5000`                      | Posterior draw cap for Brier computation (randomly subsampled if fit has more).                    |
| `US_ALL_VAL_RESULTS_DIR`         | `fits/`                     | Fits directory. **Must match** `US_ALL_VAL_RESULTS_DIR` used in run_settings_val.R.               |
| `US_ALL_AUTOSELECT_SETUP_FILE`   | `us-all-setup-auto.RData`   | Setup file path. **Must match** the setup used in run_settings_val.R (`US_ALL_VAL_SETUP_FILE`). Accepts absolute or script-relative paths. |

### Example autoselect run

> **Setup file warning:** `autoselect.R` defaults to `us-all-setup-auto.RData`
> but `run_settings_val.R` defaults to `us-all-setup-auto-200-200-400.RData`.
> Always set `US_ALL_AUTOSELECT_SETUP_FILE` to match `US_ALL_VAL_SETUP_FILE`.

```powershell
$env:US_ALL_AUTOSELECT_SETTINGS     = '3,8,12:13,15:16,29,34:35,38:39,41:42,51,54:55,58,60:62,67:68,70,74,101,105,108,111:112,114,116:118,120,124,204:206,214:216'
$env:US_ALL_AUTOSELECT_TARGET_PROBS = '0.90,0.95,0.98,0.99,0.995'
$env:US_ALL_AUTOSELECT_OBJECTIVE    = 'each'
$env:US_ALL_AUTOSELECT_SETUP_FILE   = 'us-all-setup-auto-200-200-400.RData'
Rscript autoselect.R
```

---

## File Layout

```
US-all-auto/
├── us-all-setup.RData                       # copied from ../US-all/ (source data)
├── us-all-setup-auto-<tr>-<va>-<te>.RData  # generated by us-all-setup-auto.R
├── settings-auto.csv                        # generated by settings-auto.R
├── fits/                                    # default fits dir (US_ALL_VAL_RESULTS_DIR)
│   ├── val-8.RData
│   ├── val-12.RData
│   └── ...                                  # one file per fitted setting
└── output/
    └── us-all-auto/
        ├── tables/
        │   └── autoselect_validation_scores_<suffix>.csv
        └── results/
            └── us-all-auto-results.RData
```
