# US-all MCMC 收斂診斷

論文 ozone 比較(`sec:ozonedata`,30,000 iterations / 25,000 burn-in / thin=1)
的收斂診斷:trace、R-hat、ESS、Geweke。範圍是論文表 `tab:ozone-top2-brier`
點名的 9 個 headline settings:**8, 51, 55, 58, 68, 111, 120, 124, 204**。

## 兩層設計

1. **Pass 1(單鏈,現有存檔即可)** — 每個 (setting, fold) 的存檔只有
   post-burn 5000 筆 draws 的單條 chain(`fit[[1]]`/`fit[[2]]` 是兩個 CV
   fold,不是平行鏈)。把 5000 筆切 4 段當虛擬鏈,算 rank-normalized split
   R-hat(Vehtari et al. 2021)+ coda ESS + coda Geweke,並輸出 trace 圖。
2. **Pass 2(多鏈,重跑)** — 代表 setting **55**(legacy AR1 線:skew-t,
   K=5, T=50)與 **111**(ar2 線:skew-t, K=7, T=50,論文 AR2 頭牌)。生產
   fold-1 chain 當 chain 1(seed = setting×100+1,fold 切分不耗 RNG),另跑
   3 條不同 seed + overdispersed inits 的 chain,算真正 4 鏈 R-hat。

## 腳本與執行順序

| 腳本 | 作用 |
|---|---|
| `00_conv_lib.R` | 共用函式庫(登錄表、R-hat、cache、chain 規格),source 無副作用 |
| `01_extract_chains.R` | 從 `results/us-all-<N>.RData`(每檔 0.9–1.8 GB)抽 scalar chains 存小 cache(`US_ALL_CONV_SETTINGS` 可縮範圍;`US_ALL_CONV_FORCE=1` 重抽) |
| `02_diagnostics_table.R` | 產 `convergence_diagnostics.csv`(long)+ `convergence_summary.csv` |
| `03_trace_plots.R` | 每 setting×fold 一張多面板 trace PNG(標題嵌 R-hat/ESS) |
| `04_multichain_run.R` | 重跑一條 chain:`Rscript convergence/04_multichain_run.R <setting> <chain 2-4> [dev|prod]` |
| `run-multichain.ps1` | 平行送出某 setting 的 chains 2–4(`.\convergence\run-multichain.ps1 -Setting 55`) |
| `05_multichain_diagnostics.R` | 4 鏈 R-hat/ESS 表 + 疊圖 trace + 密度圖 |

順序:01 → 02 → 03(幾分鐘);04×3(每鏈約 3–5.5 小時,三鏈平行過夜)→ 05。
所有腳本都用 PowerShell `Rscript` 執行(此機器上 Bash 跑 R 會 segfault),
工作目錄自動設為 US-all。

## 輸出

- cache:`output/us-all/results/convergence/chains-us-all-<N>.rds`、
  `multichain-<N>-fold1-chain<c>.rds`
- 表:`output/us-all/tables/convergence_diagnostics.csv`、
  `convergence_summary.csv`、`convergence_multichain.csv`
- 圖:`output/us-all/plots/convergence/trace-us-all-<N>-fold<d>.png`、
  `trace-multichain-<N>-fold1.png`、`density-multichain-<N>-fold1.png`
- log:`output/us-all/logs/convergence/multichain-<N>-chain<c>.log`

## 門檻(Vehtari et al. 2021)

- **R-hat**:< 1.01 pass;1.01–1.05 marginal;> 1.05 fail。summary 同時報
  1.01 / 1.05 / 1.10 三種門檻的 flag 數,讀者可自選(1.1 是舊 Gelman–Rubin
  慣例,現已公認過鬆)。報 `max(bulk, folded)`,bulk/tail 兩欄都保留。
- **ESS**:≥ 400 pass;< 100 severe。
- **Geweke**:|z| > 3 flag;|z| > 2 只計數——pass 1 約 250 個檢定,H₀ 下
  本來就會出現十幾個 |z|>2,要跟 R-hat/ESS 一起解讀,不能單看。

## 診斷了哪些鏈

每 fit 的 scalar 參數:`beta[*]`(204 的 MRTS 版有 11 欄)、`tau.alpha`、
`tau.beta`、`rho`、`nu`、`gamma`、`lambda`(skew 才有)、`phi.z/w/tau`
(有時間結構才有;ar2 為 φ₁、φ₂ 兩欄),外加兩條對 knot 重標記不變的場摘要
`mean.log.tau`、`mean.z`。**per-knot 的 tau/z 鏈刻意不診斷**:knots 可交換,
那些鏈會 label-switch,R-hat 沒有意義。

## 重大發現:TS 生產跑的 z 場整場凍結(2026-07-22 稽核)

Pass 1 的 degenerate 偵測揭露:**所有帶時間結構的生產配適(51, 55, 58,
111, 120, 204;68/124 無 z)之 z ≡ 0、φ_z ≡ 0,整整 30,000 次迭代一步未動。**

機制:`temporalz` 的 z 更新走 half-normal copula(`update_params.R:151`
`z.star = hn.cop(z, sig)`)。生產跑未傳 `z.init`,legacy 預設 `z.init = 0`
→ `z.star = qnorm(0) = -Inf` → 每步 Metropolis ratio 皆 NaN →
`if (!is.na(R))` 永遠跳過,z 與 φ_z(以 z.star 為資料更新)全數凍結。
ar2 backend 的 `z.init` 預設遲至 2026-05-28(commit f7770f6)才改為
half-normal 中位數,ozone 生產跑(2026 年 3–4 月)在修正之前,同樣凍結;
legacy sampler 的 `z.init = 0` 預設至今未修。

後果:
- **λ 不進 likelihood**(mu = Xβ + λ·z,z ≡ 0),整條 λ 鏈在抽先驗
  N(0, 20²)——cache 實測 λ 的 sd ≈ 20、mean ≈ 0,與先驗一致;對照
  setting 8(無 TS,z 正常更新)λ = −0.54 ± 0.07,由資料收緊。
- **預測不受 λ 汙染**:`predictY_cont_lambda`(update_params.R:1190)的
  skew 項是 `lambda * z[gp,t]`,乘上的正是凍結的 0,故 yp 與 Brier 分數
  等同 symmetric-t 模型的輸出,沒有額外雜訊。
- 換句話說,論文 ozone 比較中所有「skew-t + TS」候選其實是 symmetric-t
  配適;skew-t 與 t 的標籤差異只在無 TS 的 setting(如 8)真實存在。
- Pass 2 的重跑鏈(chains 2–4)顯式傳 `z.init > 0`,z.star 有限、更新
  正常,目標是真正的 skew-t posterior;與凍結的生產 chain 1 比,λ/mean.z
  的 4 鏈 R-hat 必然爆表——這是正確的診斷結果,不是腳本錯誤。
  `convergence_multichain.csv` 另有 `rhat_reruns` 欄(只比 chains 2–4),
  回答「解凍後的鏈彼此是否收斂」。

另:`phi.w ≡ 0` 只出現在 K=1 的 51/204,屬結構性平凡結果——partition
更新整塊被 `nknots > 1` 守衛(mcmc_cont_lambda.R:352),K=1 沒有 w 過程。

## 已知限制

- **無 acceptance rate**:兩個 sampler 的 `acc/att` 是區域變數、不在回傳
  值裡,存檔無從回推;不補救。
- **burn-in 段看不到**:sampler 回傳前就丟掉 burn-in draws,trace 只能畫
  post-burn 段(x 軸 25001–30000)。
- **legacy φ 起點寫死 0**(mcmc_cont_lambda.R:257/261/273):setting 55 的
  φ 只能透過 tau/z 場的 init 間接 overdisperse;setting 111(ar2)的 φ 可
  直接 overdisperse(`phi.*.init = c(φ1, φ2)`,已驗證在平穩三角內)。
- **ar2 backend 的 `tau.alpha`/`tau.beta` 存的是 ×2 值**(對齊稿件的 a/2,
  b/2),表格帶 `backend` 欄,勿跨 backend 直接比數值。
- **MRTS setting 204 的多鏈重跑未實作**(要另複製 `build_mrts_covariates`
  的 X 擴充邏輯);要加時擴充 `04_multichain_run.R` 即可,cache/診斷端
  已通用。
- 04 複製了 `us-all-run.R:378-491`(commit 51fa383)的 call 建構,並用
  `stopifnot` 釘住維度/iters/burn/threshold/nknots,兩邊漂移會直接炸掉
  而不是默默跑錯 posterior。
- 04 預設不抽 held-out prediction(`US_ALL_CONV_PRED=1` 可開):兩個
  sampler 都只在 `iter > burn` 抽 prediction 且不回饋參數更新,參數鏈的
  目標分布不變,省 ~10–20% 時間與記憶體。
