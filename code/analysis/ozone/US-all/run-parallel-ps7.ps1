Set-Location "D:\Github\spatial-skew-t\code\analysis\ozone\US-all"

$settings = (1..74) | Where-Object { $_ -ne 2 }
$throttle = 4               # 同時最多幾個進程（建議先 2~4）
$backend = "ar2"        # legacy / ar2
$runMode = "prod"          # dev / prod
$results = "results_new"   # 預設輸出資料夾
$wd = "D:\Github\spatial-skew-t\code\analysis\ozone\US-all"
$rscript = (Get-Command Rscript -ErrorAction Stop).Source

$settings | ForEach-Object -Parallel {
    $s = $_
    Set-Location $using:wd
    $env:US_ALL_MCMC_BACKEND = $using:backend
    $env:US_ALL_RUN_MODE = $using:runMode
    $env:US_ALL_RESULTS_DIR = $using:results

    & $using:rscript "us-all-run.R" $s
} -ThrottleLimit $throttle
