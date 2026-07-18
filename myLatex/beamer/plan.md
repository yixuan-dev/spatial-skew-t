# Defense Slides Plan (page by page)

Thesis: **Extending the Spatio-Temporal Skew-t Process for Threshold Exceedance Prediction: AR(2) Temporal Priors and Multi-Resolution Spline Basis Functions**
Student: Yi-Hsuan Hsieh   Advisor: Prof. Hao-Yun Huang

- File: [slides.tex](slides.tex)   Compile: XeLaTeX → biber → XeLaTeX ×2
- Target length: ~20 minutes (~1–1.5 min per page)
- Figures in `../img/`; bibliography shared with the thesis at `../references/references.bib`
- The slides themselves are **all English**. This plan additionally provides, for each page: a **重點整理（中文）** and a **繁體中文翻譯**. Statistical terms are kept in English throughout.

---

## Title page
- Thesis title, name, advisor, department, date, university logo.
- **重點整理**：一句話點出研究主題與兩項延伸。
- **繁體中文翻譯**：時空 skew-$t$ 過程於 threshold exceedance prediction 之延伸——AR(2) temporal priors 與 multi-resolution spline basis functions。研究生：謝易軒；指導教授：黃灝勻 博士；國立東華大學應用數學系；115 年 7 月。

## Outline
- Five sections: Motivation → Methods → Simulation → Application → Conclusion.
- **重點整理**：讓委員掌握整體結構，每章開頭會再出現一次。
- **繁體中文翻譯**：大綱——動機、方法、模擬、實證、結論。

---

## Section 1: Motivation

### Page 1 — Motivation
- Problem: threshold exceedance of ozone is spatially dependent and heavy-tailed.
- Gap: Gaussian processes cannot jointly model skewness and tail dependence.
- Two extensions: AR(2) temporal priors, multi-resolution spline basis functions.
- Goal block: evaluate predictive gain via CRPS and Brier score.
- **TODO**: add a motivating figure (ozone time series / exceedance illustration).
- **重點整理**：空污的 threshold exceedance 事件具空間相依與 heavy tail；Gaussian process 無法同時刻畫 skewness 與 tail dependence，故以 spatio-temporal skew-$t$ 過程延伸（AR(2) prior、multi-resolution spline basis），並以 CRPS、Brier score 衡量預測增益。
- **繁體中文翻譯**：
  - 空污（如 ozone）的 threshold exceedance 事件具有空間相依與 heavy-tailed 特性。
  - Gaussian process 無法同時刻畫 skewness 與 tail dependence。
  - 本研究以兩項延伸擴充 spatio-temporal skew-$t$ 過程：(1) AR(2) temporal priors，捕捉更豐富的時間動態；(2) multi-resolution spline basis functions，刻畫 non-stationary 空間結構。
  - 目標：以 CRPS、Brier score 量化上述延伸於 threshold exceedance prediction 的預測增益。

### Page 2 — Related Work
- skew-t processes and tail-dependence models; temporal priors; basis expansions and autoFRK.
- **TODO**: fill in citation keys from `references.bib` (`\textcite{...}`).
- **重點整理**：回顧三條主線：skew-$t$ process 與 tail dependence 模型、spatio-temporal 模型的 temporal prior、以及 basis function 展開與 autoFRK。
- **繁體中文翻譯**：相關研究——(1) skew-$t$ process 與 tail dependence 模型；(2) spatio-temporal 模型中的 temporal prior 設定；(3) basis function 展開與 autoFRK。

---

## Section 2: Methods

### Page 3 — Model Setup
- Observation equation: fixed effects + basis expansion + skew-t error.
- $\bB(\bs)$ multi-resolution basis, $\bomega_t$ AR(2) prior, skewness $\blambda$.
- **重點整理**：觀測 $=$ fixed effects $+$ basis expansion $+$ skew-$t$ error；其中 $\bB(\bs)$ 為 multi-resolution spline basis、$\bomega_t$ 服從 AR(2) prior、誤差帶 skewness 參數 $\blambda$。
- **繁體中文翻譯**：
  - 設 $Y(\bs,t)$ 為位置 $\bs$、時間 $t$ 的觀測值，模型為 fixed effects $+$ basis expansion $+$ skew-$t$ error。
  - $\bB(\bs)$：multi-resolution spline basis functions。
  - $\bomega_t$：temporal coefficients，服從 AR(2) prior $\bomega_t=\phi_1\bomega_{t-1}+\phi_2\bomega_{t-2}+\bfe_t$。
  - $\varepsilon(\bs,t)$：skew-$t$ error，帶 skewness 參數 $\blambda$。

### Page 4 — Inference & Evaluation
- Left: Bayesian inference (MCMC, identifiability of $\phi_\tau$). Right: predictive scoring (CRPS, Brier).
- **Key**: $\phi_\tau$ is identifiable (posterior ~4× sharper than prior) — pre-empts the "AR(2) over-parameterization" question.
- **重點整理**：以 MCMC 進行 Bayesian inference 並檢驗 $\phi_\tau$ 的 identifiability；預測以 CRPS 與 threshold exceedance 的 Brier score 評分。
- **繁體中文翻譯**：
  - Bayesian inference：以 MCMC sampling 抽樣，並檢驗 $\phi_\tau$ 的 identifiability（posterior 較 prior 約 4 倍收斂）。
  - Predictive scoring：以 CRPS 與 threshold exceedance 的 Brier score 進行評分。

---

## Section 3: Simulation

### Page 5 — Simulation Design
- Scenario table: Strong vs. Near-unit `(0.15, 0.80)`; baselines; replicates.
- **TODO**: fill the table with the final parameter values.
- **重點整理**：比較 Strong 與 Near-unit $(0.15, 0.80)$ 兩情境；AR(2) 的優勢在 near unit root 情境才顯現。
- **繁體中文翻譯**：模擬設定——比較兩情境：Strong（strong temporal dependence）與 Near-unit $(\phi_1,\phi_2)=(0.15,0.80)$（near unit root）。關鍵：AR(2) 的優勢在 near-unit 情境才顯現。

### Page 6 — Simulation Results
- Core message: AR(2) advantage emerges in near-unit (CRPS ≈ −31%, significant); no detectable difference in strong.
- **TODO**: add CRPS / Brier boxplots or table (`../img/`).
- **重點整理**：Near-unit 情境下 AR(2) 的 CRPS 約降低 $31\%$ 且顯著；Strong 情境測不到差異。（此頁待放結果圖）
- **繁體中文翻譯**：模擬結果——在 near-unit 情境下，AR(2) 的 CRPS 約降低 31% 且具顯著性；在 strong 情境下則測不到差異。（此頁待放 CRPS／Brier 結果圖）

---

## Section 4: Application

### Page 7 — Ozone Data Analysis
- Data source, preprocessing, study region and time span.
- Key finding: after removing autoFRK spatial structure, residuals heavier-tailed ($\nu$: 15.7 → 3.8, near-symmetric).
- **TODO**: add the residual QQ plot.
- **重點整理**：扣除 autoFRK 的空間結構後，residual 更 heavy-tailed（$\nu$ 由 15.7 降至 3.8）；heavy tail 是資料特徵而非 lack of fit。
- **繁體中文翻譯**：ozone 資料分析——說明資料來源、前處理與研究時空範圍。關鍵發現：扣除 autoFRK 的空間結構後，residual 更 heavy-tailed（$\nu$ 由 15.7 降至 3.8、近對稱）；此 heavy tail 是資料特徵，而非 lack of fit，故支持 skew-$t$／heavy-tailed 模型的必要性。

---

## Section 5: Conclusion

### Page 8 — Conclusion & Contributions
1. AR(2) temporal priors give significant predictive gains in the near-unit regime.
2. Multi-resolution spline bases capture non-stationary spatial structure.
3. $\phi_\tau$ is identifiable.
- Future work: higher-dimensional space-time, more threshold settings, real-time prediction.
- **重點整理**：AR(2) prior 在 near unit root 情境提供顯著預測增益；multi-resolution spline basis 刻畫 non-stationary 空間結構；$\phi_\tau$ 為 identifiable，posterior 明顯較 prior 收斂。
- **繁體中文翻譯**：
  - AR(2) temporal priors 在 near-unit 情境提供顯著預測增益。
  - Multi-resolution spline basis 刻畫 non-stationary 空間結構。
  - $\phi_\tau$ 為 identifiable，posterior 明顯較 prior 收斂。
  - 未來工作：更高維時空、更多 threshold 設定與 real-time prediction。

## References
- `\printbibliography`, with `allowframebreaks` auto-paging.

## Thanks page
- "Thank You — Questions & Comments Welcome".
- **繁體中文翻譯**：感謝聆聽，敬請指教。

---

## TODO checklist
- [ ] Pages 2 & 3: add real citations (matching keys in `references.bib`).
- [ ] Page 5: fill the scenario table with final parameter values.
- [ ] Page 6: add simulation result figures (CRPS / Brier).
- [ ] Page 7: add the ozone residual QQ plot.
- [ ] Page 1: consider a motivating figure.
- [ ] Rehearse the full talk; confirm it fits within 20 minutes.
