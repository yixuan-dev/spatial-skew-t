Set-Location "D:\Github\spatial-skew-t\code\analysis\ozone\US-all"

$settings = 1..74      # 或自己改成 1..124
$throttle = 4               # 同時最多幾個進程（建議先 2~4）
$backend = "legacy"        # legacy / ar2
$runMode = "prod"          # dev / prod
$results = "results_new"   # 預設輸出資料夾

$settings | ForEach-Object -Parallel {
    $s = $_
    Set-Location "D:\Github\spatial-skew-t\code\analysis\ozone\US-all"
    $env:US_ALL_MCMC_BACKEND = $using:backend
    $env:US_ALL_RUN_MODE = $using:runMode
    $env:US_ALL_RESULTS_DIR = $using:results
    Rscript us-all-run.R $s
} -ThrottleLimit $throttle
