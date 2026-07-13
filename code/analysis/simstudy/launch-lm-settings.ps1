<#
  launch-lm-settings.ps1 -- full-iters run of the long-memory settings 17-19.

  Runs run-settings.R once per setting (each call fans 50 datasets x methods
  {1,7,9} across -Workers PSOCK workers at the study's default iters=20000,
  burn=10000). Each fit writes results/<setting>-<method>-<dataset>.RData, so
  the run is resumable: rerun with a reduced -Datasets range to skip completed
  ones.

  Examples:
    .\launch-lm-settings.ps1                       # settings 17-19, datasets 1:50, 6 workers
    .\launch-lm-settings.ps1 -Datasets "1:10"      # quick partial
    .\launch-lm-settings.ps1 -Settings 19 -Workers 8
#>
param(
  [int[]] $Settings = @(17, 18, 19),
  [string] $Datasets = "1:50",
  [int]    $Workers  = 6,
  [string] $Methods  = "(1,7,9)"
)

$R = "C:\Program Files\R\R-4.5.1\bin\Rscript.exe"
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $here
if (-not (Test-Path logs)) { New-Item -ItemType Directory logs | Out-Null }

$stamp = Get-Date -Format "yyyyMMdd-HHmmss"
Write-Output "[$(Get-Date -Format o)] START long-memory run: settings=$($Settings -join ','), datasets=$Datasets, methods=$Methods, workers=$Workers"

foreach ($s in $Settings) {
  $log = "logs/run-lm-setting-$s-$stamp.log"
  Write-Output "[$(Get-Date -Format o)] === setting $s -> $log ==="
  & $R .\run-settings.R --setting=$s $Datasets $Workers 1 $Methods *>&1 | Tee-Object -FilePath $log
  if ($LASTEXITCODE -ne 0) {
    Write-Output "[$(Get-Date -Format o)] setting $s FAILED (exit $LASTEXITCODE); continuing to next."
  } else {
    Write-Output "[$(Get-Date -Format o)] setting $s DONE."
  }
}

Write-Output "[$(Get-Date -Format o)] ALL LONG-MEMORY SETTINGS COMPLETE."
