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
#     -Workers 8          PSOCK workers for run-settings.R (default 8)
#     -MsThreads 2        max-stable threads (unused here, default 2)
#     -MinFreeGB 25       abort a setting if free space below this
#     -KeepFits           do NOT delete fits after scoring
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
    [int]   $Workers   = 8,
    [int]   $MsThreads = 2,
    [int]   $MinFreeGB = 25,
    [switch] $KeepFits
)

$ErrorActionPreference = 'Stop'
Set-Location $PSScriptRoot                       # = code/analysis/simstudy

$data       = 'simdata_nonsta.RData'
$resultsDir = Join-Path $PSScriptRoot 'results_nonsta'
$cacheDir   = Join-Path $PSScriptRoot 'output/results'
$verifyR    = Join-Path $PSScriptRoot 'non_stationary/verify_n10_overlap.R'
$fullK      = 'c(3,5,8,10,12,15,20,25,30,40,50)'   # datasets 4:10 (all new)
$pilotNewK  = 'c(3,5,8,10,12,15,20,25)'            # datasets 1:3 (K0/30/40/50 exist)
$scoreK     = 'c(0,3,5,8,10,12,15,20,25,30,40,50)'

function Get-FreeGB {
    # results_nonsta is a subdir of the script dir, i.e. same volume.
    [math]::Round((Get-PSDrive (Get-Item $PSScriptRoot).PSDrive.Name).Free / 1GB, 0)
}
function Rrun([string]$label, [string[]]$rargs) {
    Write-Host "[$(Get-Date -Format HH:mm:ss)] $label" -ForegroundColor Cyan
    & Rscript @rargs
    if ($LASTEXITCODE -ne 0) { throw "FAILED: $label (exit $LASTEXITCODE)" }
}

if (-not (Get-Command Rscript -ErrorAction SilentlyContinue)) {
    throw "Rscript not found on PATH."
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
