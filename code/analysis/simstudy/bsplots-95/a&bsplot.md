## Record the idea of the incorperating mrts in morris model

1. As **extra covariates**

   $ x_{i,t} = [x^{\text{original}}_{i,t},\ b_1(s_i),\ldots,b_K(s_i)] $

   $ \mu_{i,t} = x_{i,t}^\top \beta + \lambda z_{i,t} $

   1. In previous work, only find the a little improvement on MSPE, mean behavior.

2. As a **separate latent block**

3. Replace the matern block

## 觀察bsplot95, bsplot98

當我們要預測**分位數95**以上的超標：

​	當真實資料是Guassian：max-stable 超爛。

​	當真實資料是max-stable：max-stable還最差。

​	當真實資料是Skew−t (K = 1, lambda = 3)：guassian, skew-t 差不多，max-stable還最差。

​	當真實資料是Skew−t (K = 5, lambda = 3)：皆差異不大的還行。

​	當真實資料是Symmetric−t (K = 1)：皆差異不大的差。

​	當真實資料是Symmetric−t (K = 5)：皆差異不大的差。

​	當真實資料是Transform below T：Skew−t, K = 1, T = q(0.0)最好，其他不差。

當我們要預測**分位數98**以上的超標：(較簡單的任務，相較預測**分位數95**以上的超標)

​	當真實資料是Guassian：皆差異不大的還行。

​	當真實資料是max-stable：皆差異不大的還行。

​	當真實資料是Skew−t (K = 1, lambda = 3)：max-stable 超爛。

​	當真實資料是Skew−t (K = 5, lambda = 3)：max-stable還最差。

​	當真實資料是Symmetric−t (K = 1)：max-stable還最差。

​	當真實資料是Symmetric−t (K = 5)：皆差異不大的差。

​	當真實資料是Transform below T：Skew−t, K = 1, T = q(0.0)最好，其他很差。

---

### 總結：

1. **錯誤率**更直觀，定義>0.5即1。todo
2. Guassian其實在Skew−t表現好的情況下也表現不差。
3. Symmetric−t最難預測
4. 分位數95比分位數98難預測

---



