set.seed(20260715)
sim_ar2 <- function(nt, phi1, phi2) {
  g1 <- phi1/(1-phi2); g2 <- phi1*g1+phi2
  se <- sqrt(1 - phi1*g1 - phi2*g2)
  L <- chol(matrix(c(1,g1,g1,1),2))
  x <- numeric(nt); x[1:2] <- as.numeric(t(L)%*%rnorm(2))
  for(t in 3:nt) x[t] <- phi1*x[t-1]+phi2*x[t-2]+rnorm(1,0,se)
  x
}
acf1 <- function(v) acf(v,lag.max=1,plot=FALSE)$acf[2]
truth <- list("9"=c(0.80,-0.35),"14"=c(0.80,-0.35),"10"=c(0.12,-0.05),"16"=c(0.15,0.80))
Tt <- 50; nnull <- 5000
cat(sprintf("%-8s %-14s %8s %10s %10s %8s\n","setting","phi_true","gamma1","E[rho1hat]","a_null","sd"))
for(st in c(9,14,10,16)){
  tp <- truth[[as.character(st)]]; g1t <- tp[1]/(1-tp[2])
  r <- replicate(nnull, acf1(sim_ar2(Tt,tp[1],tp[2])))
  cat(sprintf("%-8s (%.2f,%5.2f) %8.4f %10.4f %10.3f %8.3f\n",
      st, tp[1],tp[2], g1t, mean(r), mean(r)/g1t, sd(r)))
}
