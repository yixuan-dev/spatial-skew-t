# n30_driver.ps1 -- one-command driver for the nonsta n=30 repeats run
# (settings 20 transform_mix + 16 sigmoid_front, datasets 11-30) on
# DESKTOP-61SBCCI. Fits stay on this machine only; score caches are
# committed to git at the end.
#
# Usage (any of):
#   .\DESKTOP-61SBCCI\n30_driver.ps1                 # full run, workers = all logical CPUs
#   .\DESKTOP-61SBCCI\n30_driver.ps1 -DryRun         # print the plan + calls, run nothing
#   .\DESKTOP-61SBCCI\n30_driver.ps1 -Workers 12 -Settings "20"
#
# Interrupted? Just re-run the same command: existing valid fits are
# skipped (inventory), finished chunk caches are skipped (sanity), and
# the newest possibly-corrupt fits are detected and re-run.
#
# Per chunk of $ChunkSize datasets: fit missing cells -> score ->
# rename cache to scores<S>_nonsta_d<a>-<b>.RData -> sanity -> delete the
# chunk's fits (~18 GB per 5-dataset chunk stays the disk high-water mark).
# Per setting: merge n10.bak + chunk caches -> scores<S>_nonsta.RData
# (n=30); the merge script itself verifies every input is reproduced
# exactly (this is the 1:10 == n10.bak regression gate).

param(
    [string]$Settings   = "20,16",
    [string]$Datasets   = "11:30",
    [int]$Workers       = [Environment]::ProcessorCount,
    [int]$ChunkSize     = 5,
    [int]$DiskFloorGB   = 25,
    [switch]$DryRun,
    [switch]$SkipEnvCheck
)

$ErrorActionPreference = 'Stop'
$toolDir  = $PSScriptRoot
$simstudy = Split-Path $toolDir -Parent
Set-Location $simstudy

$stamp = Get-Date -Format 'yyyyMMdd-HHmmss'
Start-Transcript -Path (Join-Path $toolDir "n30_driver_$stamp.log") -Append | Out-Null

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

# ---- parse settings / datasets / chunks --------------------------------
$settingList = $Settings -split ',' | ForEach-Object { [int]$_.Trim() }
if ($Datasets -notmatch '^(\d+):(\d+)$') { throw "-Datasets must look like 11:30" }
$dFirst = [int]$Matches[1]; $dLast = [int]$Matches[2]
if ($dLast -lt $dFirst) { throw "-Datasets range is empty" }
$chunks = @()
for ($a = $dFirst; $a -le $dLast; $a += $ChunkSize) {
    $b = [math]::Min($a + $ChunkSize - 1, $dLast)
    $chunks += ,@($a, $b)
}
$KGrid = 'c(0,3,5,8,10,12,15,20,25,30,40,50)'
$cellsPerDataset = 5 * 12

Write-Host "=== n30 driver ==="
Write-Host "settings   : $($settingList -join ', ')"
Write-Host "datasets   : ${dFirst}:${dLast} in chunks of $ChunkSize -> $(($chunks | ForEach-Object { "$($_[0]):$($_[1])" }) -join ', ')"
Write-Host "workers    : $Workers (logical CPUs: $([Environment]::ProcessorCount))"
Write-Host "disk free  : $(Get-FreeGB) GB (floor $DiskFloorGB GB; a $ChunkSize-dataset chunk peaks at ~$($ChunkSize * 18) GB of fits)"
Write-Host "est. total : $($settingList.Count * ($dLast - $dFirst + 1) * $cellsPerDataset) fits x ~10.7 min avg / $Workers workers"
if ($DryRun) { Write-Host "*** DRY RUN: printing calls only, nothing is executed ***" }

# ---- preflight ----------------------------------------------------------
if (-not (Test-Path 'simdata_nonsta.RData')) { throw "simdata_nonsta.RData not found -- git pull incomplete?" }
# datasets starting past 10 build on the verified n=10 cache; a run from
# dataset 1 rebuilds everything and merges chunks only
$useBak = $dFirst -gt 10
if ($useBak) {
    foreach ($S in $settingList) {
        if (-not (Test-Path "output/results/scores${S}_nonsta.n10.bak.RData")) {
            throw "output/results/scores${S}_nonsta.n10.bak.RData not found -- git pull incomplete?"
        }
    }
}
if (-not (Get-Command Rscript -ErrorAction SilentlyContinue)) { throw "Rscript not found on PATH -- install R first" }
if (-not $DryRun -and -not $SkipEnvCheck) {
    Write-Host "== env check (packages + Rtools + C++ smoke test + pipeline load) =="
    Invoke-Rscript 'Rscript env_check.R'
}
Assert-DiskFloor

# ---- main loop ----------------------------------------------------------
foreach ($S in $settingList) {
    Write-Host "`n#### setting $S ####"
    $chunkCaches = @()

    foreach ($chunk in $chunks) {
        $a = $chunk[0]; $b = $chunk[1]
        $dsSpec = "${a}:${b}"
        $chunkCache = "output/results/scores${S}_nonsta_d${a}-${b}.RData"
        $chunkCaches += $chunkCache

        # resume: a chunk whose cache already exists and passes sanity is done
        if (Test-Path $chunkCache) {
            Rscript "$toolDir/chunk_sanity.R" $chunkCache $S $dsSpec
            if ($LASTEXITCODE -eq 0) { Write-Host "== chunk $dsSpec already scored, skipping"; continue }
            Write-Warning "existing $chunkCache failed sanity -- redoing this chunk"
            Remove-Item $chunkCache -Force
        }

        # -- fit: inventory -> missing calls -> execute (up to 3 rounds) --
        $callsFile = Join-Path $toolDir "calls_s${S}_d${a}-${b}.txt"
        $missing = -1
        for ($round = 1; $round -le 3; $round++) {
            $inv = & Rscript "$toolDir/inventory_gen.R" $S $dsSpec $Workers $callsFile 2>&1
            if ($LASTEXITCODE -ne 0) { $inv | Write-Host; throw "inventory_gen.R failed for setting $S $dsSpec" }
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
        if ($missing -ne 0) { throw "setting $S chunk $dsSpec still has $missing missing fits after 3 rounds" }

        # -- score the chunk, then move the fixed-path output aside --
        Invoke-Rscript ("Rscript scores.R --setting=$S --data=simdata_nonsta.RData " +
                        "--methods=`"1:5`" --datasets=`"$dsSpec`" --mrts_k=`"$KGrid`"")
        Move-Item "output/results/scores${S}_nonsta.RData" $chunkCache -Force

        # -- sanity gate, then free the chunk's disk --
        Invoke-Rscript "Rscript `"$toolDir/chunk_sanity.R`" $chunkCache $S $dsSpec"
        $removed = 0
        foreach ($d in $a..$b) {
            foreach ($pat in @("${S}-*-${d}.RData", "${S}-*-${d}-K*.RData")) {
                $files = Get-ChildItem "results_nonsta" -Filter $pat -ErrorAction SilentlyContinue
                $removed += $files.Count
                $files | Remove-Item -Force
            }
        }
        Write-Host "== chunk $dsSpec done: cache $chunkCache, deleted $removed fit files (disk free $(Get-FreeGB) GB)"
    }

    if ($DryRun) { continue }

    # -- merge n10.bak + chunks -> final n=30 cache (built-in regression) --
    $finalCache = "output/results/scores${S}_nonsta.RData"
    $inputs = if ($useBak) {
        @("output/results/scores${S}_nonsta.n10.bak.RData") + $chunkCaches
    } else { $chunkCaches }
    Invoke-Rscript ("Rscript `"$toolDir/merge_score_caches.R`" $finalCache " + ($inputs -join ' '))
    $expFirst = if ($useBak) { 1 } else { $dFirst }
    Invoke-Rscript ("Rscript -e `"load('$finalCache'); stopifnot(identical(datasets, ${expFirst}L:${dLast}L)); " +
                    "cat('final cache OK: n =', length(datasets), 'datasets\n')`"")
}

if ($DryRun) { Write-Host "`n*** DRY RUN complete ***"; Stop-Transcript | Out-Null; exit 0 }

# ---- commit score files (fits never enter git) ---------------------------
Write-Host "`n== committing score caches =="
$toAdd = @()
foreach ($S in $settingList) {
    $toAdd += "output/results/scores${S}_nonsta.RData"
    $toAdd += (Get-ChildItem "output/results/scores${S}_nonsta_d*-*.RData" | ForEach-Object {
        "output/results/$($_.Name)" })
}
git add -f $toAdd
git add $toolDir
git commit -m "Add n=30 nonsta score caches for settings $Settings (datasets $Datasets run on DESKTOP-61SBCCI)"
Write-Host "`nDONE. Now push the results back with:`n  git push origin master"
Stop-Transcript | Out-Null
