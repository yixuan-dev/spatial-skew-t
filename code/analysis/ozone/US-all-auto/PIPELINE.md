# US-all-auto Pipeline

Auto-selection experiment for spatial ozone exceedance models. Fits each
candidate model **once** on a training split and selects the best by Brier
score on a held-out validation split. No cross-validation.

---

## Pipeline Overview

```
us-all-setup-auto.R   →  us-all-setup-auto.RData   (data + 300/100/400 split)
settings-auto.R       →  settings-auto.csv          (candidate settings grid)
run_settings_val.R    →  fits/val-<N>.RData         (one MCMC fit per setting)
autoselect.R          →  output/…/tables + .RData   (best setting by val Brier)
```

### Step 1 — Build the data split

```powershell
Rscript us-all-setup-auto.R
```

Loads `us-all-setup.RData` (copied from `../US-all/`), partitions the 800
CV-covered sites into train / validation / test, and saves
`us-all-setup-auto.RData` containing:

| Object | Description |
|--------|-------------|
| `Y`, `X`, `S` | Full dataset (all 800+ sites) |
| `beta.init`, `tau.init` | MCMC initialisation values |
| `split.lst$train` | 300 site row-indices (into Y / X / S) |
| `split.lst$val` | 100 site row-indices |
| `split.lst$test` | 400 site row-indices |

Split is controlled by environment variables (see table below). Default seed
is 2024, so the split is reproducible.

| Env var | Default | Effect |
|---------|---------|--------|
| `US_ALL_AUTOSELECT_N_TRAIN` | `300` | Number of training sites |
| `US_ALL_AUTOSELECT_N_VAL` | `100` | Number of validation sites |
| `US_ALL_AUTOSELECT_N_TEST` | `400` | Number of test sites |
| `US_ALL_AUTOSELECT_SEED` | `2024` | Random seed for the partition |

### Step 2 — Generate the settings grid

```powershell
Rscript settings-auto.R
```

Reads `../US-all/settings.csv` and `../US-all/output/us-all/tables/comparison_full_table.csv`,
scores each setting by composite tail Brier, and writes `settings-auto.csv`
with tier / lane annotations.

| Tier | N | Description |
|------|---|-------------|
| 0 | 1 | Gaussian reference (setting 1) |
| 1 | 29 | Core candidates — top evidence, recommended default |
| 2 | 37 | Extended coverage — TS / AR2 / MRTS variants |
| 3 | 41 | Non-TS baseline settings, not recommended |

### Step 3 — Fit candidate models on training split

```powershell
Rscript run_settings_val.R
```

For each requested setting: fits MCMC on `split.lst$train` (300 sites),
generates posterior predictive draws at `split.lst$val` (100 sites), saves
`fits/val-<N>.RData` with:

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

| Env var | Format | Example |
|---------|--------|---------|
| `US_ALL_AUTOSELECT_SETTINGS` | Comma-separated IDs or `a:b` ranges | `'8,12,15:16,204:206'` |

Integers and range tokens can be mixed freely. The script processes them in
ascending order and skips any setting not found in `settings-auto.csv`.

### Optional

| Env var | Default | Options | Effect |
|---------|---------|---------|--------|
| `US_ALL_MCMC_BACKEND` | `legacy` | `legacy`, `ar2` | MCMC function to use. Use `ar2` for all settings including TS and AR2 lanes. |
| `US_ALL_VAL_RUN_MODE` | `prod` | `prod`, `dev` | MCMC length. `prod` = 30 000 iters / 25 000 burn. `dev` = 2 000 / 1 000 (quick check). |
| `US_ALL_VAL_RESULTS_DIR` | `fits/` | any path | Override the output directory for `val-<N>.RData` files. |

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

### Run Tier 1 — production (29 settings)

```powershell
$env:US_ALL_VAL_RUN_MODE        = 'prod'
$env:US_ALL_AUTOSELECT_SETTINGS = '8,12,15:16,38:39,41:43,55:56,58,61:62,65,67:68,70,74,105,111:112,117:118,120,124,204:206'
Rscript run_settings_val.R
```

### Run Tier 2 (37 settings)

```powershell
$env:US_ALL_AUTOSELECT_SETTINGS = '51:54,57,59:60,63:64,66,69,71:73,101:104,106:110,113:116,119,121:123,201:203,207:209'
Rscript run_settings_val.R
```

### Run Tier 1 + Tier 2 together (66 settings)

```powershell
$env:US_ALL_AUTOSELECT_SETTINGS = '8,12,15:16,38:39,41:43,51:74,101:124,201:209'
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

---

## Tier 1 Setting IDs (quick reference)

```
8, 12, 15, 16, 38, 39, 41, 42, 43, 55, 56, 58, 61, 62, 65,
67, 68, 70, 74, 105, 111, 112, 117, 118, 120, 124, 204, 205, 206
```

Range-token form: `8,12,15:16,38:39,41:43,55:56,58,61:62,65,67:68,70,74,105,111:112,117:118,120,124,204:206`

---

## autoselect.R — Environment Variables

| Env var | Default | Effect |
|---------|---------|--------|
| `US_ALL_AUTOSELECT_SETTINGS` | — (required) | Same format as run_settings_val.R; must match the settings that have been fitted |
| `US_ALL_AUTOSELECT_TARGET_PROBS` | `0.97` | Comma-separated exceedance probabilities, e.g. `'0.95,0.97,0.99'` |
| `US_ALL_AUTOSELECT_OBJECTIVE` | `single` | `single` = minimise Brier at first prob only; `mean` = minimise mean across all probs |
| `US_ALL_SUMMARY_DRAWS` | `400` | Number of posterior draws to use for Brier computation (subsampled if fit has more) |
| `US_ALL_VAL_RESULTS_DIR` | `fits/` | Must match the directory used in run_settings_val.R |

### Example autoselect run

```powershell
$env:US_ALL_AUTOSELECT_SETTINGS     = '8,12,15:16,38:39,41:43,55:56,58,61:62,65,67:68,70,74,105,111:112,117:118,120,124,204:206'
$env:US_ALL_AUTOSELECT_TARGET_PROBS = '0.97'
$env:US_ALL_AUTOSELECT_OBJECTIVE    = 'single'
Rscript autoselect.R
```

---

## File Layout

```
US-all-auto/
├── us-all-setup.RData          # copied from ../US-all/ (source data)
├── us-all-setup-auto.RData     # generated by us-all-setup-auto.R
├── settings-auto.csv           # generated by settings-auto.R
├── fits/
│   ├── val-8.RData
│   ├── val-12.RData
│   └── ...                     # one file per fitted setting
└── output/
    └── us-all-auto/
        ├── tables/
        │   └── autoselect_validation_scores.csv
        └── results/
            └── us-all-auto-results.RData
```
