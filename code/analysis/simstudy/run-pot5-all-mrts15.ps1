$ErrorActionPreference = "Stop"

$simDir = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $simDir

$logDir = Join-Path $simDir "run-logs"
New-Item -ItemType Directory -Force -Path $logDir | Out-Null

$rscript = "C:\Program Files\R\R-4.5.1\bin\Rscript.exe"
$mrtsMethods = "1:5"
$mrtsK = "15"
$scripts = @(
    "pot5-a.R",
    "pot5-b.R",
    "pot5-c.R",
    "pot5-d-1.R",
    "pot5-d-2.R",
    "pot5-d-3.R",
    "pot5-d-4.R",
    "pot5-d-5.R",
    "pot5-e.R"
)

foreach ($script in $scripts) {
    Write-Output ("[{0}] START {1}" -f (Get-Date -Format o), $script)
    if ($script -eq "pot5-e.R") {
        $env:SIMSTUDY_MRTS_METHODS = $mrtsMethods
        $env:SIMSTUDY_MRTS_K = $mrtsK
    } else {
        Remove-Item Env:SIMSTUDY_MRTS_METHODS -ErrorAction SilentlyContinue
        Remove-Item Env:SIMSTUDY_MRTS_K -ErrorAction SilentlyContinue
    }

    & $rscript (Join-Path $simDir $script)

    Write-Output ("[{0}] END {1}" -f (Get-Date -Format o), $script)
}
