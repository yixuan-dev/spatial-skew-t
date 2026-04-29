# US-all 指標計算說明

這份文件說明 `refactor_results` 目前輸出的各種指標，在程式碼中的實際計算方式。

核心原則只有一條：

- 主分數一律先在每個 fold 的 validation data 上各自計算，再在 summary 階段跨 folds 取平均

目前實作主要分成三層：

- 單一 fold 的分數計算：`01_score_engine.R`
- 跨 folds 的平均、標準誤、relative score：`01_score_engine.R`
- 匯出成 CSV：`02_comparison_tables.R`、`12_run_proposed_results.R`

## 記號

設：

- \(d = 1, \dots, D\)：cross-validation fold
- \(j\)：fold 內的 validation observation index
- \(Y_j^{(d)}\)：第 \(d\) 個 fold 中，第 \(j\) 個觀測值
- \(X_{k,j}^{(d)}\)：posterior predictive 第 \(k\) 次抽樣對同一 observation 的預測
- \(m\)：posterior predictive draws 數
- \(\tau\)：quantile level 或 threshold quantile
- \(q_\tau\)：以整體資料 `Y` 算出的 threshold value

若 observation 或預測不是有限值，該 observation 不會被計入該分數。

## 1. Quantile score

對每個 quantile level \(\tau\)，先從 posterior predictive draws 取預測分位數：

\[
\hat q_{\tau,j}^{(d)} = \text{Quantile}\left(X_{1:m,j}^{(d)}, \tau\right)
\]

接著計算 pinball loss 形式的 quantile score：

\[
QS_{\tau,j}^{(d)}
=
2 \left(\mathbf 1\{\hat q_{\tau,j}^{(d)} \ge Y_j^{(d)}\} - \tau \right)
\left(\hat q_{\tau,j}^{(d)} - Y_j^{(d)}\right)
\]

fold-level 分數是該 fold 內 observation 平均：

\[
QS_{\tau}^{(d)}
=
\frac{1}{n_d}
\sum_j QS_{\tau,j}^{(d)}
\]

最後匯出的：

\[
\text{score\_mean}_{\tau}
=
\frac{1}{D}\sum_{d=1}^D QS_{\tau}^{(d)}
\]

\[
\text{score\_se}_{\tau}
=
\frac{\operatorname{sd}(QS_{\tau}^{(1)}, \dots, QS_{\tau}^{(D)})}{\sqrt{D}}
\]

對應程式：

- `code/R/auxfunctions.R` 的 `QuantScore()`
- `refactor_results/01_score_engine.R` 的 `compute_us_all_scores()` 與 `summarize_us_all_scores()`

## 2. Brier score

對每個 threshold value \(q_\tau\)，先用 posterior draws 估計超標機率：

\[
\hat p_{\tau,j}^{(d)}
=
\frac{1}{m}\sum_{k=1}^m \mathbf 1\{X_{k,j}^{(d)} > q_\tau\}
\]

事件定義採用嚴格大於：

\[
o_{\tau,j}^{(d)} = \mathbf 1\{Y_j^{(d)} > q_\tau\}
\]

fold-level Brier score：

\[
BS_{\tau}^{(d)}
=
\frac{1}{n_d}
\sum_j \left(\hat p_{\tau,j}^{(d)} - o_{\tau,j}^{(d)}\right)^2
\]

最後匯出的：

\[
\text{score\_mean}_{\tau}
=
\frac{1}{D}\sum_{d=1}^D BS_{\tau}^{(d)}
\]

\[
\text{score\_se}_{\tau}
=
\frac{\operatorname{sd}(BS_{\tau}^{(1)}, \dots, BS_{\tau}^{(D)})}{\sqrt{D}}
\]

對應程式：

- `code/R/auxfunctions.R` 的 `BrierScore()`
- `refactor_results/01_score_engine.R` 的 `compute_us_all_scores()` 與 `summarize_us_all_scores()`

## 3. MSPE 與 MAPE

MSPE 與 MAPE 使用 posterior predictive mean 作為 point prediction。對每個 validation observation：

\[
\hat \mu_j^{(d)}
=
\frac{1}{m}\sum_{k=1}^m X_{k,j}^{(d)}
\]

fold-level MSPE 定義為：

\[
MSPE^{(d)}
=
\frac{1}{n_d}
\sum_j
\left(\hat \mu_j^{(d)} - Y_j^{(d)}\right)^2
\]

fold-level MAPE 定義為：

\[
MAPE^{(d)}
=
100 \times
\frac{1}{n_d}
\sum_j
\left|
\frac{\hat \mu_j^{(d)} - Y_j^{(d)}}{Y_j^{(d)}}
\right|
\]

其中 \(Y_j^{(d)} = 0\) 的 observation 不計入 MAPE，避免除以零。最後匯出的 `mspe_mean` 與 `mape_mean` 都是先在每個 fold 算完，再跨 folds 取平均；`mape_mean` 單位是百分比。

對應程式：

- `refactor_results/01_score_engine.R` 的 `compute_point_prediction_metrics()`
- `refactor_results/01_score_engine.R` 的 `compute_us_all_scores()` 與 `summarize_us_all_scores()`

## 4. Split Brier score

`comparison_brier_split.csv` 使用的是 same-threshold split。

對每個目標 threshold quantile

- `0.90`
- `0.95`
- `0.98`
- `0.99`
- `0.995`

都先算同一個超標事件：

\[
o_{\tau,j}^{(d)} = \mathbf 1\{Y_j^{(d)} > q_\tau\},\qquad
\hat p_{\tau,j}^{(d)} = P(Y > q_\tau \mid \text{model})
\]

然後把該 fold 的 observation 分成三組：

- `all`: 全部有效 observation
- `below_threshold`: \(Y_j^{(d)} \le q_\tau\)
- `above_threshold`: \(Y_j^{(d)} > q_\tau\)

對每一組 \(b\)，fold-level split-Brier score 定義為：

\[
BS_{\tau,b}^{(d)}
=
\frac{1}{n_{\tau,b}^{(d)}}
\sum_{j \in b}
\left(\hat p_{\tau,j}^{(d)} - o_{\tau,j}^{(d)}\right)^2
\]

其中：

- `all` 同時看 false alarm 和 miss
- `below_threshold` 只保留 \(o_{\tau,j}^{(d)} = 0\) 的 observation，偏向看 false alarm
- `above_threshold` 只保留 \(o_{\tau,j}^{(d)} = 1\) 的 observation，偏向看 miss

匯出欄位意義：

- `score_mean`

\[
\text{score\_mean}_{\tau,b}
=
\frac{1}{D}\sum_{d=1}^D BS_{\tau,b}^{(d)}
\]

- `score_se`

\[
\text{score\_se}_{\tau,b}
=
\frac{\operatorname{sd}(BS_{\tau,b}^{(1)}, \dots, BS_{\tau,b}^{(D)})}{\sqrt{D}}
\]

- `n_obs_total`

\[
n\_obs\_total_{\tau,b}
=
\sum_{d=1}^D n_{\tau,b}^{(d)}
\]

- `n_obs_mean_per_fold`

\[
n\_obs\_mean\_per\_fold_{\tau,b}
=
\frac{1}{D}\sum_{d=1}^D n_{\tau,b}^{(d)}
\]

- `obs_share`

\[
obs\_share_{\tau,b}
=
\frac{n\_obs\_total_{\tau,b}}{n\_obs\_total_{\tau,\text{all}}}
\]

注意：

- `score_mean` 是 fold 平均，不是 pooled observation 直接重算一次
- `n_obs_total` 不是理論上的 `folds × sites × days`，而是只計有限值 observation
- 在目前資料中，`n_obs_total` 變少主要是因為 `Y` 裡有 `NA`

對應程式：

- `refactor_results/01_score_engine.R` 的 `compute_split_brier_scores()`
- `refactor_results/01_score_engine.R` 的 `compute_us_all_scores()`
- `refactor_results/01_score_engine.R` 的 `summarize_us_all_scores()`
- `refactor_results/02_comparison_tables.R` 的 `build_comparison_brier_split_table()`

## 5. Relative score to Gaussian baseline

多數 comparison CSV 都會輸出 `rel_to_gaussian` 或同義欄位。

baseline 固定是 `setting = 1` 的 Gaussian reference model。

對任一 score \(S\)，relative score 定義為：

\[
\text{rel\_to\_gaussian}
=
\frac{S_{\text{model}}}{S_{\text{gaussian baseline}}}
\]

因此：

- `< 1`：比 baseline 好
- `= 1`：和 baseline 一樣
- `> 1`：比 baseline 差

在 paired table 中：

- `baseline_rel_score`：baseline model 的 relative score
- `proposed_rel_score`：proposed model 的 relative score
- `delta_rel_score = proposed_rel_score - baseline_rel_score`

## 6. Classification metrics

`comparison_classification_metrics.csv` 是事件型分類指標表。

事件定義沿用 Brier score：

\[
o_{\tau,j}^{(d)} = \mathbf 1\{Y_j^{(d)} > q_\tau\}
\]

先用 posterior predictive draws 估計超標機率：

\[
\hat p_{\tau,j}^{(d)}
=
\frac{1}{m^\star}\sum_{k=1}^{m^\star}\mathbf 1\{X_{k,j}^{(d)} > q_\tau\}
\]

其中 \(m^\star\) 是實際拿來做 classification 的 posterior predictive draws 數。  
目前實作預設使用 `400` 個等距抽樣的 draws。

再用固定機率 cutoff 把機率預報轉成二元預測：

\[
\hat o_{\tau,j}^{(d)} = \mathbf 1\{\hat p_{\tau,j}^{(d)} \ge 0.5\}
\]

於是每個 fold 內都可以累計：

- `TP`: \(\hat o=1,\, o=1\)
- `TN`: \(\hat o=0,\, o=0\)
- `FP`: \(\hat o=1,\, o=0\)
- `FN`: \(\hat o=0,\, o=1\)

對固定 threshold quantile \(\tau\)，fold-level confusion counts 記為：

\[
TP_\tau^{(d)},\;
TN_\tau^{(d)},\;
FP_\tau^{(d)},\;
FN_\tau^{(d)}
\]

fold-level 效能指標為：

\[
\text{Accuracy}_\tau^{(d)}
=
\frac{TP_\tau^{(d)} + TN_\tau^{(d)}}{TP_\tau^{(d)} + TN_\tau^{(d)} + FP_\tau^{(d)} + FN_\tau^{(d)}}
\]

\[
\text{Precision}_\tau^{(d)}
=
\frac{TP_\tau^{(d)}}{TP_\tau^{(d)} + FP_\tau^{(d)}}
\]

\[
\text{Recall}_\tau^{(d)}
=
\frac{TP_\tau^{(d)}}{TP_\tau^{(d)} + FN_\tau^{(d)}}
\]

\[
\text{Specificity}_\tau^{(d)}
=
\frac{TN_\tau^{(d)}}{TN_\tau^{(d)} + FP_\tau^{(d)}}
\]

\[
F1_\tau^{(d)}
=
\frac{2 \cdot \text{Precision}_\tau^{(d)} \cdot \text{Recall}_\tau^{(d)}}{\text{Precision}_\tau^{(d)} + \text{Recall}_\tau^{(d)}}
\]

若分母為 0，該 fold 的對應指標記為 `NA`。

匯出欄位意義：

- `tp_total`, `tn_total`, `fp_total`, `fn_total`

\[
\sum_{d=1}^D TP_\tau^{(d)},\;
\sum_{d=1}^D TN_\tau^{(d)},\;
\sum_{d=1}^D FP_\tau^{(d)},\;
\sum_{d=1}^D FN_\tau^{(d)}
\]

- `*_mean_per_fold`

\[
\frac{1}{D}\sum_{d=1}^D TP_\tau^{(d)},\quad \dots
\]

- `accuracy_mean`, `precision_mean`, `recall_mean`, `specificity_mean`, `f1_mean`

\[
\frac{1}{D}\sum_{d=1}^D M_\tau^{(d)}
\]

其中 \(M\) 代表各個 fold-level classification metric。

- `*_se`

\[
\frac{\operatorname{sd}(M_\tau^{(1)}, \dots, M_\tau^{(D)})}{\sqrt{D}}
\]

- `actual_positive_share`

\[
\frac{1}{D}\sum_{d=1}^D
\frac{TP_\tau^{(d)} + FN_\tau^{(d)}}{n_\tau^{(d)}}
\]

- `predicted_positive_share`

\[
\frac{1}{D}\sum_{d=1}^D
\frac{TP_\tau^{(d)} + FP_\tau^{(d)}}{n_\tau^{(d)}}
\]

對應程式：

- `refactor_results/01_score_engine.R` 的 `compute_classification_metrics()`
- `refactor_results/01_score_engine.R` 的 `compute_us_all_scores()` 與 `summarize_us_all_scores()`
- `refactor_results/02_comparison_tables.R` 的 `build_comparison_classification_metrics_table()`

## 7. CRPS

`comparison_scalar_metrics.csv` 的 `crps_mean` 是 sample-based CRPS。

先把 posterior draws 攤成矩陣：

\[
\mathbf X \in \mathbb R^{m \times n}
\]

其中每一欄是一個 observation 的 posterior predictive draws。

對單一 observation，CRPS 採用 sample formula：

\[
\mathrm{CRPS}(F, y)
=
\mathbb E |X - y| - \frac{1}{2}\mathbb E |X - X'|
\]

目前實作使用排序後的 sample 等價公式，以向量化方式一次算所有 observation。

fold-level CRPS 是該 fold 內 observation 平均：

\[
\mathrm{CRPS}^{(d)}
=
\frac{1}{n_d}\sum_j \mathrm{CRPS}_j^{(d)}
\]

匯出的：

- `crps_mean`：fold 平均
- `crps_se`：fold-level CRPS 的標準誤
- `crps_rel_to_gaussian`：相對於 `setting = 1` 的比值

注意：

- `summary_draws` 可能小於完整 posterior draws 數；這是為了控制計算量

## 8. Coverage

對每個 coverage level \(1-\alpha\)：

- 下界：posterior predictive 的 \(\alpha/2\) 分位數
- 上界：posterior predictive 的 \(1-\alpha/2\) 分位數

fold-level coverage：

\[
\mathrm{Coverage}_{1-\alpha}^{(d)}
=
\frac{1}{n_d}
\sum_j
\mathbf 1\{
Y_j^{(d)} \in [L_{j,\alpha}^{(d)}, U_{j,\alpha}^{(d)}]
\}
\]

匯出欄位：

- `coverage_80_mean`, `coverage_90_mean`, `coverage_95_mean`
- 對應的 `*_se`
- `*_gap = observed - target`

## 9. PIT 與 calibration summary

對單一 observation，PIT 近似為：

\[
\mathrm{PIT}_j
=
\frac{\#\{k: X_{k,j} \le Y_j\}}{m}
\]

目前程式用 `findInterval()` 依照排序後 draws 估計。

匯出的 summary 指標：

- `pit_mean`

\[
\frac{1}{n}\sum_j \mathrm{PIT}_j
\]

理想值是 `0.5`

- `pit_variance`

\[
\operatorname{Var}(\mathrm{PIT}_1,\dots,\mathrm{PIT}_n)
\]

理想值是 `1/12`

- `pit_ks_stat`
  - PIT empirical CDF 與 `Uniform(0,1)` 的 KS 距離

- `pit_uniformity_mae`

\[
\frac{1}{B}\sum_{b=1}^B |\hat s_b - s_b^\ast|
\]

其中 \(\hat s_b\) 是 PIT histogram bin 的觀測比例，\(s_b^\ast\) 是理論均勻比例

- `pit_uniformity_rmse`

\[
\sqrt{
\frac{1}{B}\sum_{b=1}^B (\hat s_b - s_b^\ast)^2
}
\]

這些都收在 `comparison_uncertainty_summary.csv`。

## 10. Calibration bins

`comparison_calibration_bins.csv` 是 PIT histogram 的 bin-level 匯出。

每個 bin 會輸出：

- `bin_left`, `bin_right`
- `expected_share`
- `observed_share_mean`
- `observed_share_se`
- `share_gap = observed_share_mean - expected_share`
- `count_total`
- `count_mean_per_fold`

用途是直接在 Excel 或其他工具畫 PIT histogram。

## 11. Top-2 table

`comparison_top2.csv` 不是新的分數，而是從 relative score 排名得出的摘要表。

對指定 target levels：

- 如果 `metric = "brier"`，用 `bs.mean.ref.gau`
- 如果 `metric = "quantile"`，用 `qs.mean.ref.gau`

然後對候選 settings 依 relative score 由小到大排序，取前兩名。

## 12. LOO-ELPD / WAIC

這兩欄目前只是 placeholder，會輸出 `NA`。

原因是現有 `fit[[d]]` 只存 posterior predictive draws `fit[[d]]$yp`，沒有 pointwise log-likelihood。

因此目前：

- `loo_elpd = NA`
- `waic = NA`

若未來 result files 額外存每個 observation 的 log-likelihood，才可以補上真正的 LOO-ELPD / WAIC。

## 13. 好壞方向

目前所有已實作的分數欄位，方向如下：

- `Quantile score`: lower is better
- `Brier score`: lower is better
- `Split Brier score`: lower is better
- `CRPS`: lower is better
- `pit_ks_stat`: lower is better
- `pit_uniformity_mae`: lower is better
- `pit_uniformity_rmse`: lower is better
- `coverage_gap`: closer to zero is better
- `pit_mean_gap`: closer to zero is better
- `pit_variance_gap`: closer to zero is better

## 14. 目前最重要的實作細節

- `score_mean` 類欄位：先各 fold 算分數，再跨 folds 取平均
- `score_se` 類欄位：fold-level 分數的標準誤
- `rel_to_gaussian` 類欄位：除以 `setting = 1`
- `n_obs_total`：跨 folds 加總的有效 observation 數
- `obs_share`：band observation 數占同 threshold 全部 observation 的比例
- threshold event 一律採 `Y > threshold`，不是 `Y >= threshold`
