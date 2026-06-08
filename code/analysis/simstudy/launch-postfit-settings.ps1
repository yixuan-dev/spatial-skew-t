# Run the simstudy post-fit pipeline for each data setting in sequence:
#   scores.R -> tables.R -> plots.R -> posterior.R
#
# Usage (from this directory):
#   .\launch-postfit-settings.ps1
#   .\launch-postfit-settings.ps1 -Settings 9..15 -Methods "c(1:5,7:10)" -Datasets "1:10"
#   .\launch-postfit-settings.ps1 -MrtsK "0"
#   .\launch-postfit-settings.ps1 -SkipExisting
#   .\launch-postfit-settings.ps1 -DryRun
#
param(
    [string[]] $Settings = 9..15,
    [string] $Methods = "c(1:5,7:10)",
    [string] $Datasets = "1:10",
    [string] $MrtsK = "0",
    [string] $DataPath = "./simdata.RData",
    [string] $Rscript = "C:\Program Files\R\R-4.5.1\bin\Rscript.exe",
    [switch] $SkipExisting,
    [switch] $DryRun
)

$ErrorActionPreference = "Stop"
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $ScriptDir

if (-not (Test-Path $Rscript)) {
    $Rscript = (Get-Command Rscript -ErrorAction SilentlyContinue).Source
    if (-not $Rscript) {
        throw "Rscript not found. Pass -Rscript or add R to PATH."
    }
}

if (-not (Test-Path $DataPath)) {
    throw "Data file not found: $DataPath"
}

$dataSuffix = ""
$base = [System.IO.Path]::GetFileNameWithoutExtension($DataPath)
if ($base -ne "simdata") {
    $dataSuffix = if ($base -match "^simdata(.+)$") { $Matches[1] } else { "_$base" }
}

$LogDir = Join-Path $ScriptDir "logs"
if (-not (Test-Path $LogDir)) {
    New-Item -ItemType Directory -Path $LogDir | Out-Null
}

$RunTag = Get-Date -Format "yyyyMMdd-HHmmss"
$SummaryLog = Join-Path $LogDir "launch-postfit-$RunTag.log"

function Write-Log {
    param([string] $Message)
    $line = "[$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')] $Message"
    Write-Host $line
    Add-Content -Path $SummaryLog -Value $line
}

function Test-SettingPipelineComplete {
    param([int] $SettingId)
    $simresults = Join-Path $ScriptDir "output/results/simresults$SettingId$dataSuffix.RData"
    $posterior = Join-Path $ScriptDir "output/results/posterior$SettingId$dataSuffix.RData"
    $posteriorCsv = Join-Path $ScriptDir "output/tables/posterior_summary$SettingId$dataSuffix.csv"
    (Test-Path $simresults) -and (Test-Path $posterior) -and (Test-Path $posteriorCsv)
}

function Invoke-RStage {
    param(
        [int] $SettingId,
        [string] $StageName,
        [string] $RScriptFile,
        [string[]] $StageArgs,
        [string] $LogFile
    )

    $cmd = @($RScriptFile) + $StageArgs
    Write-Log "  [$StageName] setting=$SettingId"
    if ($DryRun) {
        Write-Log "  DRY-RUN: & `"$Rscript`" $($cmd -join ' ')"
        return
    }

    $started = Get-Date
    $prevEap = $ErrorActionPreference
    $ErrorActionPreference = "Continue"
    & $Rscript @cmd 2>&1 | Tee-Object -FilePath $LogFile
    $rExit = $LASTEXITCODE
    $ErrorActionPreference = $prevEap
    if ($rExit -ne 0) {
        throw "${StageName} exited with code $rExit (setting=$SettingId)"
    }
    $elapsed = (Get-Date) - $started
    Write-Log ("  [$StageName] done setting=$SettingId elapsed={0:hh\:mm\:ss}" -f $elapsed)
}

Write-Log "launch-postfit-settings.ps1"
Write-Log "Rscript=$Rscript DataPath=$DataPath"
Write-Log "Settings=$($Settings -join ',') Methods=$Methods Datasets=$Datasets MrtsK=$MrtsK"

$planned = @()
$skipped = @()
$exitCode = 0

try {
foreach ($setting in $Settings) {
    $settingId = [int]$setting
    if ($SkipExisting -and (Test-SettingPipelineComplete -SettingId $settingId)) {
        Write-Log "SKIP setting=$settingId (pipeline outputs already exist)"
        $skipped += $settingId
        continue
    }

    Write-Log "START setting=$settingId"
    $planned += $settingId
    $settingLogBase = Join-Path $LogDir "launch-postfit-setting-$settingId-$RunTag"

    $commonFlags = @("--data=$DataPath", "--setting=$settingId")
    $scorePosteriorFlags = @(
        "--methods=$Methods",
        "--datasets=$Datasets"
    )
    if ($MrtsK -ne "") {
        $scorePosteriorFlags += "--mrts_k=$MrtsK"
    }

    try {
        Invoke-RStage -SettingId $settingId -StageName "scores" `
            -RScriptFile "./scores.R" `
            -StageArgs ($commonFlags + $scorePosteriorFlags) `
            -LogFile "$settingLogBase-scores.log"

        Invoke-RStage -SettingId $settingId -StageName "tables" `
            -RScriptFile "./tables.R" `
            -StageArgs $commonFlags `
            -LogFile "$settingLogBase-tables.log"

        Invoke-RStage -SettingId $settingId -StageName "plots" `
            -RScriptFile "./plots.R" `
            -StageArgs $commonFlags `
            -LogFile "$settingLogBase-plots.log"

        Invoke-RStage -SettingId $settingId -StageName "posterior" `
            -RScriptFile "./posterior.R" `
            -StageArgs ($commonFlags + $scorePosteriorFlags) `
            -LogFile "$settingLogBase-posterior.log"

        Write-Log "DONE setting=$settingId"
    } catch {
        Write-Log "FAIL setting=$settingId $_"
        $exitCode = 1
        break
    }
}

Write-Log "Finished. planned=$($planned -join ',') skipped=$($skipped -join ',')"
if ($DryRun) {
    Write-Log "Dry run only; no pipeline outputs were produced."
}
} finally {
    # no-op placeholder for symmetry with launch-ar1-methods.ps1
}

exit $exitCode
