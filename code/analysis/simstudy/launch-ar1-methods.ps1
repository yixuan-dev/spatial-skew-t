# Launch AR(1) temporal methods 9 and 10 for data settings 9..15.
#
# Each setting runs one call to run-settings.R with datasets 1:10 and methods (9,10).
# Fits are written to results/<setting>-<method>-<dataset>.RData (default simdata.RData).
#
# Usage (from this directory):
#   .\launch-ar1-methods.ps1
#   .\launch-ar1-methods.ps1 -Workers 4 -Settings 9..15 -Datasets "1:10"
#   .\launch-ar1-methods.ps1 -Workers 8 -SkipExisting
#   .\launch-ar1-methods.ps1 -DryRun
#
# On exit (success, failure, or Ctrl+C), simstudy R / PSOCK worker processes are
# cleaned up automatically. Cursor R helpServer processes are left running.
#
param(
    [string[]] $Settings = 9..15,
    [string] $Datasets = "1:10",
    [string] $Methods = "(9,10)",
    [int] $Workers = 8,
    [int] $MsThreads = 2,
    [string] $DataPath = "./simdata.RData",
    [string] $Rscript = "C:\Program Files\R\R-4.5.1\bin\Rscript.exe",
    [switch] $SkipExisting,
    [switch] $DryRun,
    [switch] $NoAutoCleanup
)

$ErrorActionPreference = "Stop"
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $ScriptDir

function Clear-SimstudyRProcesses {
    param([string] $Reason = "exit")
    $targets = Get-CimInstance Win32_Process -ErrorAction SilentlyContinue | Where-Object {
        ($_.CommandLine -like '*run-settings.R*') -or
        ($_.CommandLine -like '*parallel:::.workRSOCK*') -or
        ($_.CommandLine -like '*parallel:::.slaveRSOCK*') -or
        ($_.Name -eq 'powershell.exe' -and $_.CommandLine -like '*launch-ar1-methods*' -and $_.ProcessId -ne $PID)
    }
    if (-not $targets) { return 0 }

    $count = 0
    foreach ($proc in $targets) {
        Stop-Process -Id $proc.ProcessId -Force -ErrorAction SilentlyContinue
        $count++
    }
    if ($script:WriteLog) {
        & $script:WriteLog "CLEANUP ($Reason): stopped $count simstudy process(es)"
    } else {
        Write-Host "[cleanup] stopped $count simstudy process(es) ($Reason)"
    }
    return $count
}

if (-not $NoAutoCleanup) {
    Register-EngineEvent -SourceIdentifier PowerShell.Exiting -Action {
        Clear-SimstudyRProcesses -Reason "PowerShell.Exiting" | Out-Null
    } | Out-Null
}

if (-not (Test-Path $Rscript)) {
    $Rscript = (Get-Command Rscript -ErrorAction SilentlyContinue).Source
    if (-not $Rscript) {
        throw "Rscript not found. Pass -Rscript or add R to PATH."
    }
}

if (-not (Test-Path $DataPath)) {
    throw "Data file not found: $DataPath"
}

function Get-ResultsDir {
    param([string] $DataPath)
    $base = [System.IO.Path]::GetFileNameWithoutExtension($DataPath)
    if ($base -eq "simdata") { return "results" }
    $suffix = $base -replace "^simdata", ""
    if (-not $suffix) { return "results_$base" }
    return "results$suffix"
}

function Test-SettingComplete {
    param(
        [int] $SettingId,
        [string] $ResultsDir,
        [int[]] $MethodIds,
        [int[]] $DatasetIds
    )
    foreach ($m in $MethodIds) {
        foreach ($d in $DatasetIds) {
            $path = Join-Path $ResultsDir "$SettingId-$m-$d.RData"
            if (-not (Test-Path $path)) { return $false }
        }
    }
    return $true
}

function Resolve-DatasetIds {
    param([string] $Spec)
    $expr = $Spec.Trim()
    if ($expr -match "^(\d+):(\d+)$") {
        return [int]$Matches[1]..[int]$Matches[2]
    }
    if ($expr -match "^\((.+)\)$") {
        return ($Matches[1] -split "," | ForEach-Object { [int]$_.Trim() })
    }
    if ($expr -match "^c\((.+)\)$") {
        return ($Matches[1] -split "," | ForEach-Object { [int]$_.Trim() })
    }
    throw "Unsupported datasets spec: $Spec (use e.g. 1:10 or (1,2,6))"
}

$ResultsDir = Get-ResultsDir $DataPath
$MethodIds = @(9, 10)
$DatasetIds = Resolve-DatasetIds $Datasets

$LogDir = Join-Path $ScriptDir "logs"
if (-not (Test-Path $LogDir)) {
    New-Item -ItemType Directory -Path $LogDir | Out-Null
}

$RunTag = Get-Date -Format "yyyyMMdd-HHmmss"
$SummaryLog = Join-Path $LogDir "launch-ar1-methods-$RunTag.log"

function Write-Log {
    param([string] $Message)
    $line = "[$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')] $Message"
    Write-Host $line
    Add-Content -Path $SummaryLog -Value $line
}
$script:WriteLog = ${function:Write-Log}

if (-not $NoAutoCleanup) {
    $orphans = Clear-SimstudyRProcesses -Reason "startup"
    if ($orphans -gt 0) {
        Write-Log "Removed $orphans orphan simstudy process(es) from a prior run"
    }
}

Write-Log "launch-ar1-methods.ps1"
Write-Log "Rscript=$Rscript"
Write-Log "DataPath=$DataPath ResultsDir=$ResultsDir"
Write-Log "Settings=$($Settings -join ',') Datasets=$Datasets Methods=$Methods Workers=$Workers"

$planned = @()
$skipped = @()
$exitCode = 0

try {
foreach ($setting in $Settings) {
    $settingId = [int]$setting
    if ($SkipExisting -and (Test-SettingComplete -SettingId $settingId -ResultsDir $ResultsDir -MethodIds $MethodIds -DatasetIds $DatasetIds)) {
        Write-Log "SKIP setting=$settingId (all expected outputs exist)"
        $skipped += $settingId
        continue
    }

    $args = @(
        "./run-settings.R",
        "--data=$DataPath",
        "--setting=$settingId",
        $Datasets,
        "$Workers",
        "$MsThreads",
        $Methods
    )

    $settingLog = Join-Path $LogDir "launch-ar1-setting-$settingId-$RunTag.log"
    Write-Log "START setting=$settingId -> $settingLog"
    $planned += $settingId

    if ($DryRun) {
        Write-Log "DRY-RUN: & `"$Rscript`" $($args -join ' ')"
        continue
    }

    $started = Get-Date
    try {
        # R package load messages go to stderr; do not treat them as PowerShell errors.
        $prevEap = $ErrorActionPreference
        $ErrorActionPreference = "Continue"
        & $Rscript @args 2>&1 | Tee-Object -FilePath $settingLog
        $rExit = $LASTEXITCODE
        $ErrorActionPreference = $prevEap
        if ($rExit -ne 0) {
            throw "run-settings.R exited with code $rExit"
        }
        $elapsed = (Get-Date) - $started
        Write-Log ("DONE setting=$settingId elapsed={0:hh\:mm\:ss}" -f $elapsed)
    } catch {
        Write-Log "FAIL setting=$settingId $_"
        $exitCode = 1
        break
    }
}

Write-Log "Finished. planned=$($planned -join ',') skipped=$($skipped -join ',')"
if ($DryRun) {
    Write-Log "Dry run only; no fits were produced."
}
} finally {
    if (-not $NoAutoCleanup) {
        Clear-SimstudyRProcesses -Reason "finally" | Out-Null
    }
}

exit $exitCode
