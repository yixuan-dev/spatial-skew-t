# US-all outputs

All generated outputs for the US-all ozone workflow are collected here.

## Folder layout

- `results/` — serialized R objects and score caches (`.RData`)
- `tables/` — CSV summaries and comparison tables
- `plots/` — figures and map exports
- `logs/` — run logs or batch output

## What goes where

- Use the refactored runners in `refactor_results/`
- They write into this folder automatically
- Keep legacy root-level outputs only for historical reference

# 

- `rank_value` = **真正拿來排序、決定 Top 1 / Top 2 的數值**
- `score_value` = **你想看的原始分數本身**

**1. 當 `ranking_basis = rel_to_gaussian`**
像：
- `brier`
- `quantile`
- `crps`
- `brier_split`

這時：

- `rank_value` = 相對 Gaussian baseline 的比值
- `score_value` = 原始 mean score

例子：
如果某列是 `brier__q0.95`，
- `rank_value = 0.9856` 代表它用 `rel_to_gaussian` 來排
- `score_value = 0.02664` 代表它自己的 Brier score mean

也就是：
- 排名看相對表現
- 顯示給你看的主分數是原始分數

**2. 當 `ranking_basis = raw_score`**
像：
- `accuracy`
- `precision`
- `recall`
- `specificity`
- `f1`
- `pit_ks`
- `pit_uniformity_mae`
- `pit_uniformity_rmse`

這時通常：

- `rank_value = score_value`

因為它就是直接用原始分數排名。

例如 `accuracy`：
- `rank_value = 0.9383`
- `score_value = 0.9383`

因為 Top2 就是直接比誰 accuracy 較高。

**3. 當 `ranking_basis = abs_gap_to_target`**
像：
- `coverage`
- `pit_mean`
- `pit_variance`

這時：

- `rank_value` = 距離理想 target 的絕對差
- `score_value` = 原始 observed value

例子：
如果是 `coverage__level0.95`，
- `score_value` 可能是 `0.948`
- `rank_value` 可能是 `|0.948 - 0.95| = 0.002`

所以：
- 排名看「離理想值多近」
- 顯示分數看「實際估到多少」

**最短結論**

你可以直接這樣讀：

- `rank_value`: 排名依據
- `score_value`: 原始分數

所以它們：
- 在 `accuracy/precision/recall/specificity/f1` 這類通常一樣
- 在 `brier/quantile/crps/brier_split` 這類通常不一樣，因為排名用 `rel_to_gaussian`
- 在 `coverage/pit_mean/pit_variance` 也通常不一樣，因為排名用「距離 target 的絕對差」

如果你要，我可以下一步再幫你把 `comparison_top2` 裡這三種情況整理成一張對照表，直接放進 README。
