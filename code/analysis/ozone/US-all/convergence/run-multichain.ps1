# convergence/run-multichain.ps1
# Launch chains 2-4 of 04_multichain_run.R as three parallel detached
# Rscript processes for one setting. Logs go to
# output/us-all/logs/convergence/multichain-<setting>-chain<c>.log.
#
# Usage (from anywhere):
#   .\convergence\run-multichain.ps1 -Setting 55
#   .\convergence\run-multichain.ps1 -Setting 111 -RunMode dev
param(
    [Parameter(Mandatory = $true)][int]$Setting,
    [ValidateSet('dev', 'prod')][string]$RunMode = 'prod'
)

$rs = (Get-Command Rscript -ErrorAction SilentlyContinue).Source
if (-not $rs) { $rs = 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' }
if (-not (Test-Path $rs)) { throw "Rscript not found: $rs" }

$usall = Split-Path -Parent $PSScriptRoot
$logdir = Join-Path $usall 'output\us-all\logs\convergence'
New-Item -ItemType Directory -Force $logdir | Out-Null

foreach ($c in 2..4) {
    $log = Join-Path $logdir "multichain-$Setting-chain$c.log"
    $err = Join-Path $logdir "multichain-$Setting-chain$c.err.log"
    Start-Process -FilePath $rs `
        -ArgumentList @('convergence\04_multichain_run.R', "$Setting", "$c", $RunMode) `
        -WorkingDirectory $usall `
        -RedirectStandardOutput $log `
        -RedirectStandardError $err `
        -WindowStyle Hidden
    Write-Host "launched setting $Setting chain $c ($RunMode) -> $log"
}
Write-Host "monitor with: Get-Content -Tail 5 $logdir\multichain-$Setting-chain*.log"
