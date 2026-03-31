#!/usr/bin/env pwsh
# ============================================================================
# verify-equivalence.ps1
# 一鍵執行等價性驗證，顯示詳細報告
# ============================================================================

Set-Location "D:\Github\spatial-skew-t\code\analysis\ozone\US-all"

Write-Host ""
Write-Host "啟動等價性驗證..." -ForegroundColor Cyan

# 尋找 Rscript
$r_path = (Get-ChildItem "C:\Program Files\R" -Recurse -Filter Rscript.exe -ErrorAction SilentlyContinue | Select-Object -First 1).FullName

if ($null -eq $r_path) {
    Write-Host "✗ 找不到 Rscript.exe，請確認 R 已安裝" -ForegroundColor Red
    exit 1
}

Write-Host "✓ 找到 Rscript：$r_path" -ForegroundColor Green
Write-Host ""

# 執行驗證腳本
& $r_path --vanilla verify-equivalence.R

Write-Host ""
Write-Host "驗證完成。" -ForegroundColor Green
