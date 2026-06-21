# =====================================================================
# plots_map.R
# AQS / CMAQ / DS 三來源「臭氧設計值」地圖
# 設計值 = 4th-highest daily maximum 8-hour average ozone (EPA design value)
#
# 三張圖共用：CONUS 底圖、州界 (borders)、離散分組圖例 (bin_cols)
# 不標經緯度/座標軸；各 panel 短標題 AQS / CMAQ / DS
#
# 輸出到 myLatex/pdf/：
#   map_aqs.pdf, map_cmaq.pdf, map_ds.pdf  (3 個獨立 PDF)
#   map_design_value_vertical.pdf          (縱向 3x1 合併, 單一共用圖例)
#
# 州界風格參考 LaTeX/plots/ozone-10jul-us.pdf
# =====================================================================

library(fields)

# ---------------------------------------------------------------------
# 0. 設定區 (好替換的 seam)
# ---------------------------------------------------------------------
repo <- "d:/Github/spatial-skew-t"
ddir <- file.path(repo, "code/analysis/ozone")
outdir <- file.path(repo, "myLatex/pdf")

ds_file <- file.path(ddir, "US-all/us-all-pred-71.RData") # <- 日後換預測檔只改這行
design_rank <- 4 # 第 4 高 = EPA 設計值
ds_summary <- mean # DS 後驗點估計: mean (預設) 或 median

# --- 離散分組圖例 (共用) ---------------------------------------------
# 9 組顏色: 由使用者提供確切 hex (低 -> 高), 對應下列 9 個區間。
# 收到正式色碼後只需替換這個向量。
bin_cols <- c(
  "#bcbcbc", "#b4e4eb", "#82bcf8", "#7ddc64", "#feff00",
  "#fcbd00", "#ff3a40", "#7c2d2b", "#c228e0"
)
bin_labels <- c(
  "(-Inf,55]", "(55,60]", "(60,65]", "(65,70]", "(70,75]",
  "(75,80]", "(80,85]", "(85,90]", "(90, Inf]"
)
# 內部分界 (繪圖用有限外界取代 -Inf/Inf, 由資料範圍決定)
inner_breaks <- c(55, 60, 65, 70, 75, 80, 85, 90)
legend_title <- "July 2005\n4th Max, Daily max\n8-hour avg\nozone (ppb)"

# 取每列 (一個位置跨 31 天) 的第 design_rank 高值; NA 自動忽略
nth_max <- function(v, k = design_rank) {
  v <- sort(v, decreasing = TRUE) # sort() 預設移除 NA
  if (length(v) < k) NA else v[k]
}

# ---------------------------------------------------------------------
# 1. 載入資料與共用底圖元素 (單位皆為 km)
# ---------------------------------------------------------------------
load(file.path(ddir, "ozone_data.RData")) # x, y, CMAQ, Y, s, borders, nx, ny ...

# 框定在 CMAQ 建模域 (CONUS), 而非 range(borders) (含整個加拿大會把美國壓小)
xlim <- range(x)
ylim <- range(y)

# AQS 站點以較粗網格上色 (參考 ozone-10jul-us.pdf 風格), 才看得見
nx.aqs <- 120
ny.aqs <- 84

add_borders <- function() lines(borders) # 參考 ozone-10jul-us.pdf

# ---------------------------------------------------------------------
# 2. 計算三組 4th-max 設計值
# ---------------------------------------------------------------------

# --- AQS: 監測站觀測值 (剔除缺測 >50% -> 1089 站) ---
S <- cbind(x[s[, 1]], y[s[, 2]])
excl <- which(rowMeans(is.na(Y)) > 0.50)
S <- S[-excl, ]
Yk <- Y[-excl, ]
dv_aqs <- apply(Yk, 1, nth_max)

# --- CMAQ: 全格網 (137241 = nx * ny, 無 NA) ---
dv_cmaq <- apply(CMAQ, 1, nth_max)

# --- DS: downscaler 後驗預測 ---
# 重建預測點 S.p (同 predict-all-71.R), 再 *1000 還原成 km 對齊底圖。
# 換全域預測檔時: 須讓此處 S.p 與該檔產生時相同 (np 對應)。
# 若日後改為全格網, 可把 keep.these 設為全 TRUE, 或在預測檔旁另存 S.p。
S.p <- expand.grid(x / 1000, y / 1000)
keep.these <- (S.p[, 1] > 1.03 & S.p[, 1] < 1.7) &
  (S.p[, 2] > -0.96 & S.p[, 2] < -0.40) # 目前: 東南角
S.p <- as.matrix(S.p[keep.these, ]) * 1000 # -> km
nx.ds <- length(unique(S.p[, 1]))
ny.ds <- length(unique(S.p[, 2]))

load(ds_file) # 提供 y.pred (iters x np x nt)
yp <- y.pred
# 每抽樣每點的 4th-max, 再對抽樣取後驗摘要
dv_ds <- apply(yp, c(1, 2), nth_max) # iters x np
dv_ds <- apply(dv_ds, 2, ds_summary) # np
rm(yp, y.pred)

stopifnot(nrow(S.p) == length(dv_ds))

# ---------------------------------------------------------------------
# 3. 離散分組配色 (共用, 相同值 -> 相同顏色)
# ---------------------------------------------------------------------
all_vals <- c(dv_aqs, dv_cmaq, dv_ds)
lo <- floor(min(all_vals, na.rm = TRUE)) - 1 # 有限外界取代 -Inf
hi <- ceiling(max(all_vals, na.rm = TRUE)) + 1 # 有限外界取代 +Inf
breaks <- c(lo, inner_breaks, hi)
stopifnot(length(bin_cols) == length(breaks) - 1)

# ---------------------------------------------------------------------
# 4. 單一 panel 繪圖工具 (重用; 不畫 panel 內圖例)
# ---------------------------------------------------------------------
panel_aqs <- function() {
  quilt.plot(S[, 1], S[, 2], dv_aqs,
    nx = nx.aqs, ny = ny.aqs,
    breaks = breaks, col = bin_cols, add.legend = FALSE,
    main = "AQS", axes = FALSE, xlab = "", ylab = "",
    xlim = xlim, ylim = ylim
  )
  add_borders()
}

panel_cmaq <- function() {
  quilt.plot(rep(x, times = ny), rep(y, each = nx), dv_cmaq,
    nx = nx, ny = ny,
    breaks = breaks, col = bin_cols, add.legend = FALSE,
    main = "CMAQ", axes = FALSE, xlab = "", ylab = "",
    xlim = xlim, ylim = ylim
  )
  add_borders()
}

panel_ds <- function() {
  quilt.plot(S.p[, 1], S.p[, 2], dv_ds,
    nx = nx.ds, ny = ny.ds,
    breaks = breaks, col = bin_cols, add.legend = FALSE,
    main = "DS", axes = FALSE, xlab = "", ylab = "",
    xlim = xlim, ylim = ylim
  )
  add_borders()
}

# 離散圖例 (在當前繪圖區放置)
draw_legend <- function(pos = "right", cex = 0.9) {
  legend(pos,
    legend = bin_labels, fill = bin_cols, title = legend_title,
    bty = "n", cex = cex, xpd = NA
  )
}

# 以全裝置座標覆蓋, 放單一圖例 (避免與 layout/mfrow 衝突)
overlay_legend <- function(pos = "right", cex = 0.9) {
  par(
    fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0),
    new = TRUE
  )
  plot.new()
  draw_legend(pos, cex)
}

# ---------------------------------------------------------------------
# 5. 輸出
# ---------------------------------------------------------------------
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# --- 3 個獨立 PDF (各自附離散圖例於右側) ---
make_single <- function(file, panel) {
  pdf(file, width = 8, height = 5)
  par(oma = c(0, 0, 0, 6), mar = c(1, 1, 2, 1))
  panel()
  overlay_legend("right")
  dev.off()
}
make_single(file.path(outdir, "map_aqs.pdf"), panel_aqs)
make_single(file.path(outdir, "map_cmaq.pdf"), panel_cmaq)
make_single(file.path(outdir, "map_ds.pdf"), panel_ds)

# --- 縱向 3x1 合併 PDF (單一共用離散圖例) ---
# 註: 用 mfrow 排三列, 再以全裝置座標覆蓋一個圖例, 避免 layout 衝突。
draw_combined <- function() {
  par(oma = c(0, 0, 0, 6)) # 右側預留圖例空間
  par(mfrow = c(3, 1), mar = c(1, 1, 2, 1))
  panel_aqs()
  panel_cmaq()
  panel_ds()
  overlay_legend("right")
}
pdf(file.path(outdir, "map_design_value_vertical.pdf"), width = 6.5, height = 13)
draw_combined()
dev.off()

cat("Done. Wrote 4 PDFs to", outdir, "\n")
cat(sprintf("breaks (finite): lo=%.1f .. hi=%.1f ppb\n", lo, hi))

# --- 驗證用 PNG 快照 (僅供目視, 寫到暫存) ---
qlpng <- Sys.getenv("PLOTS_MAP_PNG", "")
if (nzchar(qlpng)) {
  png(qlpng, width = 760, height = 1500)
  draw_combined()
  dev.off()
  cat("Wrote quicklook PNG:", qlpng, "\n")
}
