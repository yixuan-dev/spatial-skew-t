# >>> Insert above: # rm(list = ls())
# ======================================================================
# [00] OPTIONAL LEGACY BATCH ASSEMBLY (COMMENTED BLOCK)
# Purpose: Assemble CV metrics/chains from per-model result files and save
#          consolidated objects (kept as historical reproducibility code).
# ======================================================================


# >>> Insert above: rm(list = ls())  # first active analysis section
# ======================================================================
# [01] INITIALIZE ANALYSIS SESSION
# Purpose: Reset workspace, load helper functions/data/setup/results, and
#          unpack saved score objects for downstream comparison.
# ======================================================================


# >>> Insert above: quant.score.mean <- matrix(NA, 74, length(probs))
# ======================================================================
# [02] COMPUTE MODEL-WISE SUMMARY STATISTICS
# Purpose: Compute mean and standard error of Quantile/Brier scores across
#          CV splits for each model setting.
# ======================================================================


# >>> Insert above: quant.score.mean[c(1:10, 13, 14, 17:19, 21, 23, 25, 26), ]
# ======================================================================
# [03] IDENTIFY BEST MODELS BY THRESHOLD
# Purpose: Report model indices achieving minimum average Quantile/Brier
#          score at each threshold level.
# ======================================================================


# >>> Insert above: bs.mean.ref.gau <- matrix(NA, nrow = 73, ncol = 11)
# ======================================================================
# [04] BUILD RELATIVE-TO-GAUSSIAN PERFORMANCE MATRICES
# Purpose: Normalize each model’s mean score by Gaussian baseline (model 1)
#          to compare relative gains/losses.
# ======================================================================


# >>> Insert above: # spatially plot brier scores for q(0.95): Gaussian, ...
# ======================================================================
# [05] SITE-LEVEL BRIER SCORES FOR SELECTED MODELS (1, 3, 8, 71)
# Purpose: Compute per-site Brier scores at high threshold probability and
#          store spatially aligned performance vectors.
# ======================================================================


# >>> Insert above: windows(width = 10, height = 14)
# ======================================================================
# [06] SPATIAL BRIER MAPS AND EXPORT
# Purpose: Visualize site-level Brier score surfaces and export publication-
#          ready PDF comparisons.
# ======================================================================


# >>> Insert above: load("results/us-all-16.RData")
# ======================================================================
# [07] ADDITIONAL SITE-LEVEL COMPARISON (MODELS 16, 36)
# Purpose: Compute alternative per-site Brier score diagnostics for selected
#          candidate models.
# ======================================================================


# >>> Insert above: # get tau such that q(tau) = 75 for each site
# ======================================================================
# [08] LINK SITE OZONE REGIME TO MODEL PERFORMANCE
# Purpose: Derive marginal site behavior around 75 ppb and examine how
#          relative Brier performance changes across site regimes.
# ======================================================================


# >>> Insert above: fit.np <- npreg(brier.score.site[, 1] ~ ozone.quant.site,
# ======================================================================
# [09] NONPARAMETRIC TREND FIT
# Purpose: Smooth Brier score vs site ozone regime relationship for visual
#          trend interpretation.
# ======================================================================


# >>> Insert above: # find top two for selected quantiles
# ======================================================================
# [10] TOP-MODEL LOOKUP AT KEY THRESHOLDS
# Purpose: Rank and inspect best/second-best model settings for selected
#          threshold quantiles.
# ======================================================================


# >>> Insert above: # three main plots (keep max-stable in all for now)
# ======================================================================
# [11] COMPARATIVE SUMMARY PLOTS: TIME SERIES × THRESHOLD × KNOTS
# Purpose: Build line-based summaries to compare model families under
#          multiple design factors.
# ======================================================================


# >>> Insert above: # another set of plots, 1 time series, 1 no time series
# ======================================================================
# [12] PAPER PANEL FIGURES (BS-OZONE)
# Purpose: Create side-by-side “time-series vs non-time-series” panels and
#          export final manuscript figures.
# ======================================================================


# >>> Insert above: # one knot
# ======================================================================
# [13] EXTENDED EXPLORATORY PANELS (LEGACY FIGURES)
# Purpose: Produce detailed model-by-model exploratory plots across knot and
#          threshold combinations (diagnostic/appendix style).
# ======================================================================


# >>> Insert above: wilcox.results.gau <- matrix(NA, nrow = length(probs), ncol = 4)
# ======================================================================
# [14] STATISTICAL SIGNIFICANCE CHECKS
# Purpose: Run Wilcoxon tests to assess pairwise performance differences in
#          Brier scores across selected settings.
# ======================================================================


# >>> Insert above: library(fields)  # second appearance near heatmap section
# ======================================================================
# [15] HEATMAP-STYLE SCORE SUMMARY
# Purpose: Visualize selected quantile-score blocks as a matrix image for
#          compact comparison across (K, T) combinations.
# ======================================================================


# >>> Insert above: # posterior predictions
# ======================================================================
# [16] POSTERIOR PREDICTION WORKFLOW: SETUP
# Purpose: Reset session, load ozone grid/data, define prediction region, and
#          configure exceedance threshold.
# ======================================================================


# >>> Insert above: load("us-all-pred-1.RData")
# ======================================================================
# [17] POSTERIOR SUMMARIES — MODEL SET 1 (GAUSSIAN)
# Purpose: Compute site-level 95th/99th posterior quantiles and probabilities
#          of at least 1/2/3 exceedance days.
# ======================================================================


# >>> Insert above: # 1 knot - No Time Series - T = 0
# ======================================================================
# [18] POSTERIOR SUMMARIES — MODEL SET 3
# Purpose: Repeat posterior quantile/exceedance summaries for set 3.
# ======================================================================


# >>> Insert above: # Skew-t - No Time series - T = 50
# ======================================================================
# [19] POSTERIOR SUMMARIES — MODEL SET 8
# Purpose: Repeat posterior quantile/exceedance summaries for set 8.
# ======================================================================


# >>> Insert above: # 6 knots - Time series - T = 75
# ======================================================================
# [20] POSTERIOR SUMMARIES — MODEL SET 59
# Purpose: Repeat posterior quantile/exceedance summaries for set 59.
# ======================================================================


# >>> Insert above: # 10 knots - Time series - T = 75
# ======================================================================
# [21] POSTERIOR SUMMARIES — MODEL SET 71
# Purpose: Repeat posterior quantile/exceedance summaries for set 71.
# ======================================================================


# >>> Insert above: save.image(file = "predict-maps.RData")  # first appearance
# ======================================================================
# [22] SAVE / MERGE PREDICTION-MAP OBJECTS
# Purpose: Persist and combine posterior summary objects for unified map
#          rendering.
# ======================================================================


# >>> Insert above: # make the prediction maps
# ======================================================================
# [23] GENERATE AND EXPORT PREDICTION MAPS
# Purpose: Draw exceedance probability maps, posterior quantile maps, and
#          model-difference maps; export figure PDFs.
# ======================================================================


# >>> Insert above: # load("us-all-setup.RData")  # large commented diagnostic block
# ======================================================================
# [24] OPTIONAL MCMC TRACE DIAGNOSTICS (COMMENTED BLOCK)
# Purpose: Archived chain-trace diagnostics for latent processes/parameters
#          under multiple model runs.
# ======================================================================


# >>> Insert above: # get an idea of two sites that are close to one another vs far apart
# ======================================================================
# [25] CLOSE-VS-FAR SITE DEPENDENCE DIAGNOSTICS
# Purpose: Compare nearby vs distant stations using distance filters,
#          rank-style daily position summaries, and paired visual diagnostics.
# ======================================================================