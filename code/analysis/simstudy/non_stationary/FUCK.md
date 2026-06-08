FUCK

我的目的是比較加 mrts covariates 後的表現，那我不應該觀察 brier score relative to gaussian method ，而應該要觀察 brier score relative to original method before adding mrts covariates 。那我們可以做一張折線圖，x軸放由小到大的 mrts_k ，y軸放 relative brier score，不同的 legend 代表不同的 methods。那這張圖可以看出甚麼：在哪些數量的 mrts 擁有最低的 relative brier score，以及 mrts covariates 對該method是否有加分；但不能看出甚麼：即便 a method with mrts covariates which beats the original method ，能否在該 data setting 下打敗 other original methods with no mrts covariates 是看不出來的。

