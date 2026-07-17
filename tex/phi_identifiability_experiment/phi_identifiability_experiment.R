###############################################################################
# FINAL phi-identifiability experiment (verified code path).
###############################################################################
set.seed(20260715)
results_dir <- "d:/Github/spatial-skew-t/code/analysis/simstudy/results"
truth <- list("9"=c(0.80,-0.35),"14"=c(0.80,-0.35),"10"=c(0.12,-0.05),"16"=c(0.15,0.80))
g1_true <- function(p) p[1]/(1-p[2])
in_S <- function(p1,p2) abs(p2)<1 & (p1+p2)<1 & (p2-p1)<1
acf1 <- function(v) acf(v,lag.max=1,plot=FALSE,na.action=na.pass)$acf[2]
thin <- 20

# prior SD restricted to S
pp <- matrix(rnorm(6e6,0,0.5),ncol=2); pp <- pp[in_S(pp[,1],pp[,2]),]
prior_sd <- apply(pp,2,sd)

# finite-sample ACF null: AR(2) at truth, T=50
sim_ar2 <- function(nt,p1,p2){g1<-p1/(1-p2);g2<-p1*g1+p2;se<-sqrt(1-p1*g1-p2*g2)
  L<-chol(matrix(c(1,g1,g1,1),2));x<-numeric(nt);x[1:2]<-as.numeric(t(L)%*%rnorm(2))
  for(t in 3:nt)x[t]<-p1*x[t-1]+p2*x[t-2]+rnorm(1,0,se);x}

rows <- list()
for (st in c(9,14,10,16)) {
  tp <- truth[[as.character(st)]]; g1t <- g1_true(tp)
  pm1<-pm2<-ps1<-ps2<-numeric(0); rho1_ds<-numeric(0)  # per-dataset aggregates
  for (ds in 1:10) {
    f <- file.path(results_dir, sprintf("%d-8-%d.RData",st,ds)); if(!file.exists(f)) next
    e<-new.env(); load(f,envir=e); fit<-e$fit.1; rm(e)
    idx<-seq(1,dim(fit$tau)[1],by=thin); tau<-fit$tau[idx,,,drop=FALSE]
    al<-fit$tau.alpha[idx]/2; be<-fit$tau.beta[idx]/2; nd<-length(idx); K<-dim(tau)[2]
    ts<-array(NA_real_,dim(tau)); for(d in seq_len(nd)) ts[d,,]<-qnorm(pgamma(tau[d,,],shape=al[d],rate=be[d]))
    ts[!is.finite(ts)]<-NA
    rd<-numeric(0); for(d in seq_len(nd)) for(k in seq_len(K)) rd<-c(rd,acf1(ts[d,k,]))
    rho1_ds <- c(rho1_ds, mean(rd,na.rm=TRUE))       # this dataset's mean per-draw rho1
    ph<-fit$phi.tau
    pm1<-c(pm1,mean(ph[,1]));pm2<-c(pm2,mean(ph[,2]));ps1<-c(ps1,sd(ph[,1]));ps2<-c(ps2,sd(ph[,2]))
    rm(fit,tau,ts);gc(verbose=FALSE)
  }
  rho1 <- mean(rho1_ds); a_obs <- rho1/g1t
  rnull <- replicate(5000, acf1(sim_ar2(50,tp[1],tp[2]))); a_null <- mean(rnull)/g1t
  rows[[as.character(st)]] <- data.frame(setting=st, phi_true=sprintf("(%.2f,%.2f)",tp[1],tp[2]),
    phi1_hat=mean(pm1), postSD1=mean(ps1), betaSD1=sd(pm1), ratio1=mean(ps1)/prior_sd[1], z1=(mean(pm1)-tp[1])/mean(ps1),
    phi2_hat=mean(pm2), postSD2=mean(ps2), betaSD2=sd(pm2), ratio2=mean(ps2)/prior_sd[2], z2=(mean(pm2)-tp[2])/mean(ps2),
    g1=g1t, rho1=rho1, a_obs=a_obs, a_null=a_null, a_adj=a_obs/a_null)
}
res<-do.call(rbind,rows)
cat(sprintf("prior SD|_S = (%.3f, %.3f)\n\n", prior_sd[1], prior_sd[2]))
cat("=== (C) phi posterior vs prior vs truth ===\n")
print(format(res[,c("setting","phi_true","phi1_hat","postSD1","ratio1","z1","phi2_hat","postSD2","ratio2","z2")],digits=3),row.names=FALSE)
cat("\n=== (D) attenuation vs finite-sample null ===\n")
print(format(res[,c("setting","phi_true","g1","rho1","a_obs","a_null","a_adj")],digits=3),row.names=FALSE)
saveRDS(res, sub("final_experiment.R","final_result.rds", "C:/Users/Yi-Xuan/AppData/Local/Temp/claude/d--Github-spatial-skew-t/e32df0d9-5c8a-4547-a357-9ec95b24a907/scratchpad/final_experiment.R"))
