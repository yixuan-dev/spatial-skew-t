# US-all settings recap

This note summarizes the **main differences** across all settings listed in `settings.csv`.

## Column meanings (what changes across settings)

- `method`: model family (`gaussian`, `max-stable`, `skew-t`, `t`)
- `knots`: number of spatial knots
- `thresh`: threshold level (typically 0, 50, 75, 85)
- `CMAQ`: whether CMAQ covariates are used (`yes` in all listed settings except the `a` variants)
- `TS`: temporal model switch (`yes` or `no`)
- `ar2`: AR(2) temporal extension flag (`yes` for settings 101–124)
- `rerun` / `running`: execution bookkeeping (not model-spec differences)

---

## High-level setting families

## 1) Baseline/early comparison settings (1–6a)

- `1`: Gaussian baseline (`knots=1`, `thresh=0`, `TS=no`)
- `2`: Max-stable baseline (`knots=1`, `thresh=75`, `TS=no`)
- `3,4,5,6`: CMAQ=`yes`, TS=`no`, `knots=1`, comparing:
  - `skew-t` at `thresh=0,50,85`
  - `t` at `thresh=75`
- `5a,6a`: same model forms as 5/6 but **CMAQ=`no`** (ablation checks)

## 2) Non-TS spatial grid expansion (7–50)

- TS remains `no`
- Expands knot counts from `{2,3,4,5,6,7,8,9,10,11,15}`
- Main comparison pattern:
  - `skew-t` at lower thresholds (`0`, `50`, some `75`/`85`)
  - `t` at `75` (and some `85` variants for `t` in 19–26)
- `9a,10a,13a,14a`: paired **CMAQ=`no`** counterparts of selected settings

## 3) TS-enabled production grid (51–74)

- Same core model comparison idea as non-TS grid, but with **`TS=yes`**
- Typical triplet per knot count:
  - `skew-t, thresh=0`
  - `skew-t, thresh=50`
  - `t, thresh=75`
- Knot counts covered: `{1,5,6,7,8,9,10,15}`

## 4) AR(2) extension grid (101–124)

- Mirrors TS grid structure of 51–74, with **`ar2=yes`**
- `CMAQ=yes`, `TS=yes`, `rerun=yes` throughout
- Same triplet per knot count `{1,5,6,7,8,9,10,15}`:
  - `skew-t` at `thresh=0`
  - `skew-t` at `thresh=50`
  - `t` at `thresh=75`

---

## Canonical triplet pattern (used repeatedly)

For many knot counts `k`, the settings follow:

1. `skew-t`, `knots=k`, `thresh=0`
2. `skew-t`, `knots=k`, `thresh=50`
3. `t`,      `knots=k`, `thresh=75`

This appears in:
- TS AR(1)-style block: `51–74`
- TS AR(2) block: `101–124`

---

## Practical interpretation

If you compare results by blocks:

- **Method effect**: `skew-t` vs `t` vs early baselines (`gaussian`, `max-stable`)
- **Threshold effect**: mostly `0/50` (skew-t) vs `75` (t), plus some `85` stress settings
- **Knot effect**: increasing spatial complexity via `knots`
- **Temporal effect**: `TS=no` (7–50) vs `TS=yes` (51–74)
- **AR order effect**: TS block (51–74) vs AR(2) block (101–124)
- **CMAQ covariate effect**: compare specific `a` variants (`CMAQ=no`) against their paired settings

---

## Quick reference by setting ranges

- `1–2`: special baselines (`gaussian`, `max-stable`)
- `3–50`: TS=`no` exploration and ablations
- `51–74`: TS=`yes` main grid
- `101–124`: TS=`yes` + `ar2=yes` extension
