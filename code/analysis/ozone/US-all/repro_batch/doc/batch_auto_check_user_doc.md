# `batch_auto_check.R` User Guide

This document explains how to use:

- `repro_batch/batch_auto_check.R`

for batch execution and health-check of `US-all` experiment scripts (`us-all-<setting>.R`).

Current implementation executes settings through the shared runner:

- `us-all-run.R <setting>`

---

## What this script is for

`batch_auto_check.R` helps you:

1. Run a set of settings (or predefined batch) automatically.
2. Skip already completed results unless forced.
3. Do dry checks without execution.
4. Save run status and failure summaries to `repro_batch/output/`.
5. Run full baseline batches `F1 -> F13` (settings 1-74) in one command.

---

## Expected project structure

The script assumes it is run from or near:

- `code/analysis/ozone/US-all/repro_batch/batch_auto_check.R`

and it writes/reads from:

- results directory: `code/analysis/ozone/US-all/results/`
- output directory: `code/analysis/ozone/US-all/repro_batch/output/`

---

## Predefined batch map

Built-in batches in the script:

- `F1`: `1:2` (2 settings)
- `F2`: `3:8` (6 settings)
- `F3`: `9:14` (6 settings)
- `F4`: `15:20` (6 settings)
- `F5`: `21:26` (6 settings)
- `F6`: `27:32` (6 settings)
- `F7`: `33:38` (6 settings)
- `F8`: `39:44` (6 settings)
- `F9`: `45:50` (6 settings)
- `F10`: `51:56` (6 settings)
- `F11`: `57:62` (6 settings)
- `F12`: `63:68` (6 settings)
- `F13`: `69:74` (6 settings)
- `F14`: `101:108` (8 AR2 settings)
- `F15`: `109:116` (8 AR2 settings)
- `F16`: `117:124` (8 AR2 settings)
- `AR2_SMOKE`: `101,102,103` (smoke test)
- `A1`: alias of `AR2_SMOKE`

Default run (no args) uses `F1` settings.

---

## Command-line options

### Main selection options

- `--settings=...`  
  Run specific integer settings, e.g. `--settings=33,34,35`

- `--batch=<name>`  
  Run a predefined batch from the map above.

- `--full-done`  
  Run baseline sequence `F1 -> F13` (all non-AR2 settings 1-74).

> Rules:
>
> - `--settings` and `--batch` cannot be used together.
> - `--full-done` cannot be combined with `--settings` or `--batch`.

### Behavior options

- `--dry-run`  
  Check file presence only; do not run scripts.

- `--force`  
  Re-run even if `results/us-all-<setting>.RData` already exists.

- `--label=<text>`  
  Override run label used in output filenames.

### Full-done failure policy options

- `--stop-on-batch-failure`  
  In `--full-done` mode, stop after first failed batch.

- `--continue-on-error`  
  Explicitly continue through all batches even if one fails.

> You cannot use `--stop-on-batch-failure` and `--continue-on-error` together.

---

## Recommended usage patterns

## 1) Quick dry check for a batch

Use this to verify missing outputs before launching long runs.

`Rscript repro_batch/batch_auto_check.R --batch=F4 --dry-run`

## 2) Run one predefined batch

`Rscript repro_batch/batch_auto_check.R --batch=F4 --label=F4`

## 3) Run custom settings

`Rscript repro_batch/batch_auto_check.R --settings=57,58,59 --label=custom-F4`

## 4) Full baseline pipeline (F1–F13)

`Rscript repro_batch/batch_auto_check.R --full-done`

## 5) Full pipeline but stop on first failed batch

`Rscript repro_batch/batch_auto_check.R --full-done --stop-on-batch-failure`

## 6) AR2 smoke check

`Rscript repro_batch/batch_auto_check.R --batch=AR2_SMOKE --label=AR2-SMOKE --stop-on-batch-failure`

---

## Windows example (explicit Rscript path)

If `Rscript` is not on `PATH`, use explicit executable path:

`& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" ".\repro_batch\batch_auto_check.R" --batch=F6 --label=F6 --stop-on-batch-failure`

---

## Output files and where to read them

For each run label `<LABEL>`, script writes:

- timestamped status CSV:  
  `repro_batch/output/<LABEL>-run-status-YYYYmmdd-HHMMSS.csv`
- latest status CSV (overwrite):  
  `repro_batch/output/<LABEL>-run-status-latest.csv`
- timestamped failures text:  
  `repro_batch/output/<LABEL>-failures-YYYYmmdd-HHMMSS.txt`
- latest failures text (overwrite):  
  `repro_batch/output/<LABEL>-failures-latest.txt`

In `--full-done`, each checkpoint batch (`F1`…`F13`) writes its own latest files, and a final aggregated `full-done` report is also written.

---

## `run_status` values (status CSV)

- `ok`: runner ran and result file exists
- `already_present`: result existed; script not rerun (unless `--force`)
- `dry_present`: dry-run saw result file present
- `dry_missing`: dry-run saw result file missing
- `script_missing`: `us-all-run.R` not found
- `missing_result`: runner returned but expected `.RData` not found
- `error`: runner raised error and no result file exists
- `error_result_exists`: runner error occurred but result file exists

Failure statuses tracked in summary are:

- `script_missing`
- `missing_result`
- `error`
- `error_result_exists`
- `dry_missing`

---

## Typical workflow

1. Start with dry run (`--dry-run`) for your target batch.
2. Execute batch without `--dry-run`.
3. Inspect `*-run-status-latest.csv` and `*-failures-latest.txt`.
4. Re-run only failed settings with `--settings=...`.
5. For complete baseline reproduction, run `--full-done`.

---

## Common pitfalls

- Mixing mutually exclusive options (`--settings` + `--batch`, or `--full-done` + `--batch`).
- Forgetting `--force` when you actually want to overwrite existing results.
- Running from an unexpected working directory without a valid script path context.
- Non-integer values in `--settings` (only integers are accepted).

---

## Minimal troubleshooting checklist

- Check `repro_batch/output/<LABEL>-run-status-latest.csv` first.
- If many `script_missing`: verify `us-all-run.R` exists in project root.
- If many `missing_result`: inspect `us-all-run.R` and `settings.csv` mapping for save path/result object logic.
- If many `error`: open `error_message` column in status CSV and rerun those settings individually.

---

## Notes on reproducibility

- The script is orchestration + health-check only; scientific comparability is determined by the underlying model scripts and scoring aggregation logic.
- For same-basis model comparison, pair this runner with your comparison aggregator workflow (e.g., baseline fixed + proposed lane).
