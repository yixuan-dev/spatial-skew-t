# =============================================================================
# tzeng2018_scenario1_simulate.R
# -----------------------------------------------------------------------------
# 根據 Tzeng & Huang (2018) "Resolution Adaptive Fixed Rank Kriging"
# Section 4.1 的 Scenario 1（非定態, nonstationary）生成模擬資料。
#
# 模型：
#   z(s, t) = y(s, t) + eps(s, t)
#   y(s, t) = w1(t) * f1(s) + w2(t) * f2(s) + xi(s, t)
#
#   f1(s) = cos(   pi * || s - (0,   1)   || )
#   f2(s) = cos( 2*pi * || s - (3/4, 1/4) || )
#   (w1(t), w2(t))' ~ N(0, diag(25, 9))
#   eps(s, t)       ~ N(0, sigma2_eps),   sigma2_eps = 3
#   xi(s, t)        ~ N(0, sigma2_xi),    sigma2_xi  = 0  (預設關閉；論文未明寫此值)
#
# 設定：ns = 50, nt = 20, D = [0,1]^2, s_i 由 simple random sampling 抽出。
# 輸出：RData (完整物件)。
# =============================================================================

# ---- 0. 設定 ----------------------------------------------------------------
set.seed(1234)

ns <- 144
nt <- 50
sigma2_eps <- 3
sigma2_xi <- 0 # 若要加上 small-scale white-noise，把這裡改成正值
M <- diag(c(25, 9)) # Var((w1, w2))

c1 <- c(0, 1) # f1 中心
c2 <- c(3 / 4, 1 / 4) # f2 中心

# ---- 1. 抽位置 s_1, ..., s_n （[0,1]^2 simple random sampling）-------------
s <- cbind(x = runif(ns), y = runif(ns)) # ns x 2

# ---- 2. 基底函數 ------------------------------------------------------------
f1 <- cos(pi * sqrt((s[, 1] - c1[1])^2 + (s[, 2] - c1[2])^2))
f2 <- cos(2 * pi * sqrt((s[, 1] - c2[1])^2 + (s[, 2] - c2[2])^2))
F_mat <- cbind(f1 = f1, f2 = f2) # ns x 2

# ---- 3. 隨機係數 w_t (nt x 2) -----------------------------------------------
L <- chol(M) # M = L'L
W <- matrix(rnorm(nt * 2), nt, 2) %*% L # nt x 2

# ---- 4. 真值 y_t  與 觀測 z_t -----------------------------------------------
Y <- F_mat %*% t(W) # ns x nt
if (sigma2_xi > 0) {
  Y <- Y + matrix(
    rnorm(ns * nt, sd = sqrt(sigma2_xi)),
    ns, nt
  )
}
E <- matrix(rnorm(ns * nt, sd = sqrt(sigma2_eps)), ns, nt)
Z <- Y + E # ns x nt

# ---- 5. 存 ------------------------------------------------------------
rda_fp <- "non_sta_bf.RData"
save(
  s, F_mat, W, Z, ns, nt, M, c1, c2,
  file = rda_fp
)

cat("Done.\n",
  "  RData : ", rda_fp, "\n",
  sep = ""
)
