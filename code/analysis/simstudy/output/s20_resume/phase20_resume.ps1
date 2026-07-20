$ErrorActionPreference = 'Stop'
Set-Location d:\Github\spatial-skew-t\code\analysis\simstudy
$scratch = 'd:\Github\spatial-skew-t\code\analysis\simstudy\output\s20_resume'

$calls = Get-Content 'output/s20_resume_calls.txt'
$i = 0
foreach ($call in $calls) {
    $i++
    $free = [math]::Round((Get-PSDrive D).Free / 1GB, 0)
    if ($free -lt 25) { throw "disk below 25 GB before call $i -- stopping" }
    Write-Output ("[resume {0}/{1}] (disk {2} GB) {3}" -f $i, $calls.Count, $free, $call)
    Invoke-Expression $call
    if ($LASTEXITCODE -ne 0) { throw "resume call $i failed: $call" }
}

$n = (Get-ChildItem results_nonsta -Filter '20-*.RData').Count
Write-Output "[fits] setting 20 files on disk: $n (expect 600)"
if ($n -lt 600) { throw "expected 600 fits, found $n" }

# ---- score n=10 full grid with cache protection ----
$orig = "output/results/scores20_nonsta.RData"
$bak  = "output/results/scores20_nonsta.n3.bak.RData"
if (-not (Test-Path $bak)) { Copy-Item $orig $bak }
Rscript scores.R --setting=20 --data=simdata_nonsta.RData --methods="1:5" --datasets="1:10" --mrts_k="c(0,3,5,8,10,12,15,20,25,30,40,50)"
if ($LASTEXITCODE -ne 0) { throw "scores.R failed (old n=3 cache safe in $bak)" }

# ---- regression check: new cache must reproduce the old n=3 cells ----
Rscript "$scratch\s20_regression.R"
if ($LASTEXITCODE -ne 0) {
    Copy-Item $bak $orig -Force
    throw "regression check FAILED -- old n=3 cache restored, fits kept for inspection"
}

Write-Output "PHASE 20 RESUME DONE (regression passed; fits kept until user deletes)"
