# expA_hn_driver.ps1 -- one-command driver for the time-block Experiment A
# run (5 expanding-window blocks) on settings 5 and 7, datasets 1:10,
# methods 1 (i.i.d.) / 2 (AR(2)) / 4 (AR(1)), under lambda ~ HN(0, 20).
# Fits stay on this machine and are deleted after scoring; the score
# caches and the diagnostics ledger are committed to git at the end.
#
# Usage (from code/analysis/time_block_forecast/simstudy):
#   .\DESKTOP-61SBCCI\expA_hn_driver.ps1 -DryRun      # print the plan, run nothing
#   .\DESKTOP-61SBCCI\expA_hn_driver.ps1              # full run, workers auto
#   .\DESKTOP-61SBCCI\expA_hn_driver.ps1 -Workers 10 -Settings "5"
#   .\DESKTOP-61SBCCI\expA_hn_driver.ps1 -KeepFits    # gate but never delete fits
#
# Interrupted? Re-run the same command. Datasets whose cache exists and
# passes the gate are skipped entirely; valid fits are skipped by the
# inventory; a fit truncated by a kill is detected and re-run.
#
# Per fit batch of $ChunkSize datasets: fit missing cells -> then per
# dataset: score -> gate -> delete that dataset's 3 fits (~2.4 GB). Disk
# high-water mark is one batch, ~12 GB at the default ChunkSize.
# Per setting: merge the 10 per-dataset caches into scores<S>.RData; the
# merge script verifies every input is reproduced exactly.
#
# Assertion C is NOT a gate here: the HN prior is the protection against
# the lambda reflected ridge, and C_spread / B_truth are recorded as a
# descriptive ledger only. Nothing is excluded on their basis.

param(
    [string]$Settings   = "5,7",
    [string]$Datasets   = "1:10",
    [string]$Methods    = "c(1,2,4)",
    [int]$Workers       = 0,          # 0 = auto (CPU and RAM aware)
    [double]$GbPerWorker = 2.0,
    [int]$ChunkSize     = 5,
    [int]$EsDraws       = 1000,
    [int]$DiskFloorGB   = 40,
    [switch]$DryRun,
    [switch]$SkipEnvCheck,
    [switch]$NoCommit,
    # Keep every fit on disk after gating instead of deleting it. Full
    # campaign is ~49 GB (60 cells x 0.81 GB) -- affordable insurance
    # against a future re-scoring, which the 2026-07-31 threshold-rule
    # change showed can otherwise force a refit (the pred_* mirrors cannot
    # resolve exceedance probabilities at a new threshold).
    [switch]$KeepFits
)

$ErrorActionPreference = 'Stop'
$toolDir  = $PSScriptRoot
$study    = Split-Path $toolDir -Parent
Set-Location $study

$stamp = Get-Date -Format 'yyyyMMdd-HHmmss'
Start-Transcript -Path (Join-Path $toolDir "expA_hn_driver_$stamp.log") -Append | Out-Null

# the backend commit that introduced lambda.positive; without it --hn is a
# silently ignored argument
$HN_COMMIT = '4ec3628'
$FIT_GB    = 0.81      # one 5-block cell on disk

function Invoke-Rscript([string]$CmdLine) {
    Write-Host ">> $CmdLine"
    Invoke-Expression $CmdLine
    if ($LASTEXITCODE -ne 0) { throw "command failed (exit $LASTEXITCODE): $CmdLine" }
}

function Get-FreeGB {
    $drive = (Get-Location).Drive
    if ($null -eq $drive) { return $null }   # UNC path, no drive letter
    [math]::Round((Get-PSDrive $drive.Name).Free / 1GB, 1)
}

function Assert-DiskFloor {
    $free = Get-FreeGB
    if ($null -eq $free) { Write-Warning "cannot determine free disk space (UNC path?)"; return }
    if ($free -lt $DiskFloorGB) { throw "free disk $free GB is below the $DiskFloorGB GB floor -- stopping" }
}

# ---- parse settings / datasets / batches --------------------------------
$settingList = $Settings -split ',' | ForEach-Object { [int]$_.Trim() }
if ($Datasets -notmatch '^(\d+):(\d+)$') { throw "-Datasets must look like 1:10" }
$dFirst = [int]$Matches[1]; $dLast = [int]$Matches[2]
if ($dLast -lt $dFirst) { throw "-Datasets range is empty" }
$allDatasets = $dFirst..$dLast
$batches = @()
for ($a = $dFirst; $a -le $dLast; $a += $ChunkSize) {
    $b = [math]::Min($a + $ChunkSize - 1, $dLast)
    $batches += ,@($a, $b)
}
$methodIds = (Invoke-Expression ($Methods -replace '^c', '@')) | ForEach-Object { [int]$_ }
$nMethods = $methodIds.Count

# ---- machine survey + worker cap ---------------------------------------
$cpu = [Environment]::ProcessorCount
$ramGB = [math]::Round((Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory / 1GB, 1)
$capCpu = [math]::Max(1, $cpu - 1)
$capRam = [math]::Max(1, [math]::Floor(($ramGB - 6) / $GbPerWorker))
if ($Workers -le 0) { $Workers = [math]::Min($capCpu, $capRam) }

$cells = $settingList.Count * $allDatasets.Count * $nMethods
Write-Host "=== Experiment A (5 blocks), lambda ~ HN(0,20) ==="
Write-Host "host       : $env:COMPUTERNAME"
Write-Host "settings   : $($settingList -join ', ')"
Write-Host "datasets   : ${dFirst}:${dLast} in fit batches of $ChunkSize -> $(($batches | ForEach-Object { "$($_[0]):$($_[1])" }) -join ', ')"
Write-Host "methods    : $Methods"
Write-Host "cpu / ram  : $cpu logical cores / $ramGB GB"
Write-Host "workers    : $Workers  (cap: cpu $capCpu, ram $capRam at $GbPerWorker GB/worker)"
if ($capRam -lt $capCpu) {
    Write-Warning "RAM is the binding constraint: $capRam workers, not $capCpu. Fitting will not use all cores."
}
Write-Host "disk free  : $(Get-FreeGB) GB (floor $DiskFloorGB GB; a $ChunkSize-dataset batch peaks at ~$([math]::Round($ChunkSize * $nMethods * $FIT_GB, 1)) GB of fits)"
Write-Host "est. total : $cells cells x ~1.8 h = ~$([math]::Round($cells * 1.8, 0)) core-hours"
Write-Host "es_draws   : $EsDraws (energy/variogram only; CRPS and Brier keep all draws)"
if ($DryRun) { Write-Host "*** DRY RUN: printing the plan and the calls, nothing is executed ***" }

# ---- preflight ----------------------------------------------------------
foreach ($f in @('simdata.RData', 'run-settings.R', 'scores.R', 'tables.R',
                 'time_block_helpers.R', 'fit_diag_utils.R',
                 'collect_diag.R', 'expA_threeway.R',
                 '../../simstudy/DESKTOP-61SBCCI/merge_score_caches.R')) {
    if (-not (Test-Path $f)) { throw "$f not found -- git pull incomplete?" }
}
if (-not (Get-Command Rscript -ErrorAction SilentlyContinue)) { throw "Rscript not found on PATH -- install R first" }

# --hn is inert unless the backend carries lambda.positive
git merge-base --is-ancestor $HN_COMMIT HEAD 2>$null
if ($LASTEXITCODE -ne 0) {
    throw "backend commit $HN_COMMIT (lambda.positive) is not in HEAD -- git pull origin master first"
}
Write-Host "git sha    : $(git rev-parse --short HEAD)  (includes $HN_COMMIT)"

if (-not $DryRun -and -not $SkipEnvCheck) {
    Write-Host "== env check (packages + Rtools + C++ smoke test) =="
    Invoke-Rscript 'Rscript ../../simstudy/env_check.R'
}
Assert-DiskFloor

# ---- main loop ----------------------------------------------------------
foreach ($S in $settingList) {
    Write-Host "`n#### setting $S ####"
    $caches = @()

    foreach ($batch in $batches) {
        $a = $batch[0]; $b = $batch[1]
        $dsSpec = "${a}:${b}"

        # which datasets of this batch still need work?
        $todo = @()
        foreach ($d in $a..$b) {
            $cache = "output/results/scores${S}_d${d}.RData"
            if (Test-Path $cache) {
                Rscript "$toolDir/chunk_sanity_tbf.R" $cache $S "$d" "$Methods" $EsDraws
                if ($LASTEXITCODE -eq 0) { Write-Host "== s$S d$d already scored and gated, skipping"; continue }
                Write-Warning "existing $cache failed the gate -- redoing dataset $d"
                Remove-Item $cache -Force
            }
            $todo += $d
        }
        if ($todo.Count -eq 0) { continue }

        # -- fit: inventory -> missing calls -> execute (up to 3 rounds) --
        $callsFile = Join-Path $toolDir "calls_s${S}_d${a}-${b}.txt"
        $missing = -1
        for ($round = 1; $round -le 3; $round++) {
            $inv = & Rscript "$toolDir/inventory_tbf.R" $S $dsSpec "$Methods" $Workers $callsFile 2>&1
            if ($LASTEXITCODE -ne 0) { $inv | Write-Host; throw "inventory_tbf.R failed for setting $S $dsSpec" }
            $inv | Write-Host
            $missing = [int]($inv | Select-String '^MISSING (\d+)$').Matches[0].Groups[1].Value
            if ($missing -eq 0) { break }
            if ($DryRun) {
                Write-Host "-- DryRun: would execute these calls --"
                Get-Content $callsFile | Write-Host
                break
            }
            foreach ($call in Get-Content $callsFile) {
                Assert-DiskFloor
                Invoke-Rscript $call
            }
        }
        if ($DryRun) { continue }
        if ($missing -ne 0) { throw "setting $S batch $dsSpec still has $missing missing fits after 3 rounds" }

        # -- score / gate / delete, one dataset at a time -----------------
        # Scoring costs ~1 min per dataset at EsDraws = 1000, so it stays
        # sequential; the finer unit means each committed cache is an
        # independently verified atom and disk is freed sooner.
        foreach ($d in $a..$b) {
            $cache = "output/results/scores${S}_d${d}.RData"
            $caches += $cache
            if (Test-Path $cache) { continue }   # gated above

            Invoke-Rscript ("Rscript scores.R --setting=$S --methods=`"$Methods`" " +
                            "--datasets=`"$d`" --es_draws=$EsDraws --out=`"$cache`"")

            Invoke-Rscript "Rscript `"$toolDir/chunk_sanity_tbf.R`" $cache $S $d `"$Methods`" $EsDraws"

            if ($KeepFits) {
                Write-Host "== s$S d$d gated, fits kept (-KeepFits; disk free $(Get-FreeGB) GB)"
            } else {
                $removed = 0
                foreach ($m in $methodIds) {
                    $f = "results/$S-$m-$d.RData"     # exact name, never a wildcard
                    if (Test-Path $f) { Remove-Item $f -Force; $removed++ }
                }
                Write-Host "== s$S d$d gated, $removed fits deleted (disk free $(Get-FreeGB) GB)"
            }
        }
    }

    if ($DryRun) { continue }

    # -- merge the per-dataset caches -> final cache for this setting -----
    $caches = $allDatasets | ForEach-Object { "output/results/scores${S}_d$_.RData" }
    $finalCache = "output/results/scores${S}.RData"
    Invoke-Rscript ("Rscript `"../../simstudy/DESKTOP-61SBCCI/merge_score_caches.R`" " +
                    "$finalCache " + ($caches -join ' '))
    # brier_threshold_basis + the anyNA checks close merge_score_caches.R's
    # silent NA-fill path: a dataset-array absent from one input is NA-filled,
    # not rejected, so the merged cache must be re-checked for completeness.
    Invoke-Rscript ("Rscript -e `"load('$finalCache'); " +
                    "stopifnot(identical(as.integer(datasets), ${dFirst}L:${dLast}L), isTRUE(hn_prior), " +
                    "identical(brier_threshold_basis, 'full_series'), " +
                    "!anyNA(crps.lead), !anyNA(brier.lead), !anyNA(brier.lead.blockq), " +
                    "!anyNA(pexceed.mean), !anyNA(brier.thresholds), !anyNA(exceed.rate.lead)); " +
                    "cat('final cache OK: n =', length(datasets), 'datasets,', length(methods), 'methods\n')`"")

    # -- downstream artifacts (cheap, and they travel back in git) --------
    Invoke-Rscript "Rscript collect_diag.R --setting=$S"
    Invoke-Rscript "Rscript expA_threeway.R --setting=$S"
    Invoke-Rscript "Rscript tables.R --setting=$S"
    Invoke-Rscript "Rscript plots.R --setting=$S"
}

if ($DryRun) { Write-Host "`n*** DRY RUN complete ***"; Stop-Transcript | Out-Null; exit 0 }
if ($NoCommit) { Write-Host "`n-NoCommit: skipping the git step"; Stop-Transcript | Out-Null; exit 0 }

# ---- commit the caches and the ledger (fits never enter git) ------------
Write-Host "`n== committing score caches and diagnostics =="
$toAdd = @()
foreach ($S in $settingList) {
    $toAdd += "output/results/scores${S}.RData"
    $toAdd += (Get-ChildItem "output/results/scores${S}_d*.RData" | ForEach-Object { "output/results/$($_.Name)" })
    $toAdd += (Get-ChildItem "output/diag" -Filter "diag_${S}-*.csv" | ForEach-Object { "output/diag/$($_.Name)" })
    $toAdd += (Get-ChildItem "output/diag" -Filter "post_${S}-*.csv" | ForEach-Object { "output/diag/$($_.Name)" })
    $toAdd += (Get-ChildItem "output/diag" -Filter "*_${S}.csv" | ForEach-Object { "output/diag/$($_.Name)" })
    $toAdd += (Get-ChildItem "output/tables" -Filter "*${S}.csv" | ForEach-Object { "output/tables/$($_.Name)" })
    $toAdd += (Get-ChildItem "output/plots" -Filter "*set${S}.pdf" | ForEach-Object { "output/plots/$($_.Name)" })
}
# pred_*.RData is ~47 MB of predictive summaries: kept on this machine as
# the insurance against a future re-scoring, deliberately NOT committed.
git add -f $toAdd
git add $toolDir
git commit -m "Experiment A (5 blocks): settings $Settings x methods $Methods x datasets $Datasets, lambda~HN(0,20), run on $env:COMPUTERNAME"
Write-Host "`nDONE. Now push the results back with:`n  git push origin master"
Stop-Transcript | Out-Null
