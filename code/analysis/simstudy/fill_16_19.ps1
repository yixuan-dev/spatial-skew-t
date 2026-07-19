# =====================================================================
# fill_16_19.ps1  --  bring non-t settings 16, 19 to full coverage
#                     (n=10, K in {0,3,5,8,10,12,15,20,25,30,40,50})
#                     then score and DELETE the fits (disk-limited).
#
# PORTABLE: run from any machine that has this repo folder and Rscript
# on PATH.  Paths are resolved relative to this script's own location,
# so just copy the repo (or link/mount it) and run:
#
#     pwsh -File code/analysis/simstudy/fill_16_19.ps1
#
# Optional parameters:
#     -Settings 16,19     which settings to fill (default 16,19)
#     -Workers 0          PSOCK workers for run-settings.R.
#                         0 (default) = AUTO: logical cores - 1, additionally
#                         capped at floor(RAM_GB / 3) so each worker has
#                         ~3 GB headroom. Pass an explicit number to override.
#     -MsThreads 2        max-stable threads (unused here, default 2)
#     -MinFreeGB 25       abort a setting if free space below this
#     -KeepFits           do NOT delete fits after scoring
#     -SkipEnvCheck       skip env_check.R (packages/Rtools/compile test)
#     -SkipTest           skip the single-fit serial test run
#
# Startup sequence (fail fast, nothing destructive before the real run):
#   1. env_check.R  -- packages present (auto-installs), Rtools/make works,
#                      C++ kernel compiles, ar2 pipeline loads.
#   2. test run     -- ONE fit, serial (1 worker): first setting, dataset 4,
#                      method 1, K=0.  Catches data/path/compile issues in
#                      ~6 min instead of failing 8 workers deep.
#   3. full run     -- the parallel increment per setting, as planned.
#
# Per setting it (a) fits only the INCREMENT over the existing pilot
# (datasets 1:3 x K{0,30,40,50} already on disk are reused, not refit --
# the seed does not depend on K), (b) backs up the old n=3 score cache,
# (c) re-scores to n=10 / full-K (OVERWRITING scores<ID>_nonsta.RData),
# (d) verifies the overlap cells are bit-identical to the pilot, and only
# then (e) deletes the fits.  A verify failure keeps the fits + backup.
# =====================================================================

param(
    [int[]] $Settings  = @(16, 19),
    [int]   $Workers   = 0,        # 0 = auto-detect from this machine's CPU/RAM
    [int]   $MsThreads = 2,
    [int]   $MinFreeGB = 25,
    [switch] $KeepFits,
    [switch] $SkipEnvCheck,
    [switch] $SkipTest
)

$ErrorActionPreference = 'Stop'
Set-Location $PSScriptRoot                       # = code/analysis/simstudy

if ($Workers -lt 1) {
    # Auto-size to THIS machine: one PSOCK worker per logical core, minus one
    # core for the OS/master process, capped so each worker keeps ~3 GB RAM.
    $cores  = [Environment]::ProcessorCount
    $ramGB  = 64   # assume plenty if RAM cannot be queried
    try {
        $ramGB = [math]::Floor((Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory / 1GB)
    } catch {}
    $Workers = [math]::Max(1, [math]::Min($cores - 1, [math]::Floor($ramGB / 3)))
    Write-Host "auto workers: $Workers  (logical cores = $cores, RAM = ${ramGB} GB)" -ForegroundColor Yellow
}

$data       = 'simdata_nonsta.RData'
$resultsDir = Join-Path $PSScriptRoot 'results_nonsta'
$cacheDir   = Join-Path $PSScriptRoot 'output/results'
$verifyR    = Join-Path $PSScriptRoot 'non_stationary/verify_n10_overlap.R'
$fullK      = 'c(3,5,8,10,12,15,20,25,30,40,50)'   # datasets 4:10 (all new)
$pilotNewK  = 'c(3,5,8,10,12,15,20,25)'            # datasets 1:3 (K0/30/40/50 exist)
$scoreK     = 'c(0,3,5,8,10,12,15,20,25,30,40,50)'

function Get-FreeGB {
    # results_nonsta is a subdir of the script dir, i.e. same volume.
    # Works for drive letters AND UNC paths (\\server\share\...): PSDrive has
    # no Free for UNC, so fall back to the Win32 GetDiskFreeSpaceEx via .NET.
    try {
        $d = (Get-Item $PSScriptRoot).PSDrive
        if ($null -ne $d -and $null -ne $d.Free) { return [math]::Round($d.Free / 1GB, 0) }
    } catch {}
    try {
        $free = 0L; $total = 0L; $totalFree = 0L
        $sig = @'
[System.Runtime.InteropServices.DllImport("kernel32.dll", SetLastError = true, CharSet = System.Runtime.InteropServices.CharSet.Auto)]
public static extern bool GetDiskFreeSpaceEx(string lpDirectoryName, out long lpFreeBytesAvailable, out long lpTotalNumberOfBytes, out long lpTotalNumberOfFreeBytes);
'@
        $k32 = Add-Type -MemberDefinition $sig -Name 'DiskFree' -Namespace 'Win32' -PassThru
        if ($k32::GetDiskFreeSpaceEx($PSScriptRoot, [ref]$free, [ref]$total, [ref]$totalFree)) {
            return [math]::Round($free / 1GB, 0)
        }
    } catch {}
    Write-Warning "cannot determine free space for $PSScriptRoot -- skipping the disk guard."
    return [int]::MaxValue
}
function Rrun([string]$label, [string[]]$rargs) {
    Write-Host "[$(Get-Date -Format HH:mm:ss)] $label" -ForegroundColor Cyan
    & Rscript @rargs
    if ($LASTEXITCODE -ne 0) { throw "FAILED: $label (exit $LASTEXITCODE)" }
}

if (-not (Get-Command Rscript -ErrorAction SilentlyContinue)) {
    throw "Rscript not found on PATH."
}

# ---- 1. environment check (packages, Rtools, compile, pipeline load) ----
if (-not $SkipEnvCheck) {
    Rrun "environment check (env_check.R)" @('env_check.R')
}

# ---- 2. single-fit serial test run --------------------------------------
# One real fit (first setting, dataset 4, method 1, K=0) with 1 worker.
# It is part of the needed set anyway; the full run refits it with the same
# seed, so the only cost is ~6 minutes -- cheap insurance before going wide.
if (-not $SkipTest) {
    $s0 = $Settings[0]
    Rrun "TEST RUN: setting $s0, dataset 4, method 1, K=0, serial" `
        @('run-settings.R', "--data=$data", "--setting=$s0", '4', '1', "$MsThreads", '1')
    $probe = Join-Path $resultsDir "$s0-1-4.RData"
    if (-not (Test-Path $probe)) { throw "test run produced no output ($probe missing)." }
    Write-Host "TEST RUN OK -> $probe  (proceeding to the parallel run)" -ForegroundColor Green
}

foreach ($s in $Settings) {
    $free = Get-FreeGB
    Write-Host "==== setting $s : free ${free} GB ====" -ForegroundColor Yellow
    if ($free -lt ($MinFreeGB + 151)) {
        throw "setting ${s}: only ${free} GB free; need ~$($MinFreeGB + 151) GB for the increment. Free space and retry."
    }

    # (a) fit the increment ------------------------------------------------
    Rrun "setting $s : datasets 4:10 baseline (K=0)" `
        @('run-settings.R', "--data=$data", "--setting=$s", '4:10', "$Workers", "$MsThreads", '1:5')
    Rrun "setting $s : datasets 4:10 full-K sweep" `
        @('run-settings.R', "--data=$data", "--setting=$s", '4:10', "$Workers", "$MsThreads", '1:5', "$fullK")
    Rrun "setting $s : datasets 1:3 new-K only" `
        @('run-settings.R', "--data=$data", "--setting=$s", '1:3', "$Workers", "$MsThreads", '1:5', "$pilotNewK")

    # (b) back up the old n=3 cache ----------------------------------------
    $cache  = Join-Path $cacheDir "scores${s}_nonsta.RData"
    $backup = Join-Path $cacheDir "scores${s}_nonsta.pre_n10.bak.RData"
    if (-not (Test-Path $cache)) { throw "setting ${s}: expected pilot cache $cache not found." }
    Copy-Item $cache $backup -Force

    # (c) re-score to n=10 / full-K (overwrites $cache) --------------------
    Rrun "setting $s : score n=10 full-K" `
        @('scores.R', "--setting=$s", "--data=$data", '--methods=1:5', '--datasets=1:10', "--mrts_k=$scoreK")

    # (d) verify overlap vs pilot ------------------------------------------
    Write-Host "[$(Get-Date -Format HH:mm:ss)] setting $s : verifying overlap" -ForegroundColor Cyan
    & Rscript $verifyR "$s" $backup $cache
    if ($LASTEXITCODE -ne 0) {
        throw "setting ${s}: overlap verification FAILED. Kept fits + backup ($backup). Investigate before deleting."
    }

    # (e) delete fits (unless -KeepFits); drop the now-redundant backup ----
    if ($KeepFits) {
        Write-Host "setting $s : -KeepFits set, fits retained." -ForegroundColor Green
    } else {
        $n = (Get-ChildItem $resultsDir -Filter "$s-*.RData").Count
        Get-ChildItem $resultsDir -Filter "$s-*.RData" | Remove-Item -Force -Confirm:$false
        Write-Host "setting $s : deleted $n fit files. Free now $(Get-FreeGB) GB." -ForegroundColor Green
    }
    Remove-Item $backup -Force
    Write-Host "==== setting $s DONE (scores${s}_nonsta.RData is now n=10 / full-K) ====" -ForegroundColor Green
}

Write-Host "ALL SETTINGS DONE." -ForegroundColor Green
Write-Host "Next (on the repo machine): regenerate tables + recompile --" -ForegroundColor Gray
Write-Host "  Rscript non_stationary/bestK_tables.R" -ForegroundColor Gray
Write-Host "  Rscript non_stationary/brier_byK_table.R" -ForegroundColor Gray
