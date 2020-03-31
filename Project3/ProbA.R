# setwd("C:/Users/vrekk/OneDrive/Dokumenter/NTNU/Vaar4/Beregningskrevende/Project/Project3")
library(matrixStats)
source("probAhelp.R")
source("probAdata.R")

x0 <- data3A$x
n <- length(x0)
plot(1:n,x0)

betahat <- ARp.beta.est(x0,2)
res_LS <- ARp.resid(x0,betahat$LS)
res_LA <- ARp.resid(x0,betahat$LA)

### a)

B <- 1500
n_res <- length(res_LS)
#Bootstrap x for LS
xb_LS <- matrix(NA,n,B)
for (b in 1:B){
  res <- sample(res_LS,n_res,replace=TRUE)
  x <- ARp.filter(x0[rep(sample(99,1),2)+c(0,1)],betahat$LS,res)
  xb_LS[,b] <- x
}

#Bootstrap x for LA
xb_LA <- matrix(NA,n,B)
for (b in 1:B){
  res <- sample(res_LA,n_res,replace=TRUE)
  x <- ARp.filter(x0[rep(sample(99,1),2)+c(0,1)],betahat$LA,res)
  xb_LA[,b] <- x
}

#Compute beta_LS from bootstrap sample
betahatb_LS <- matrix(NA,2,B)
for (b in 1:B){
  betahatb_LS[,b] <- ARp.beta.est(xb_LS[,b],2)$LS
}
#Estimate variance of beta_LS
var_beta_LS <- rowVars(betahatb_LS)
bias_beta_LS <- rowMeans(betahatb_LS) - betahat$LS

#Compute beta_LA from bootstrap sample
betahatb_LA <- matrix(NA,2,B)
for (b in 1:B){
  betahatb_LA[,b] <- ARp.beta.est(xb_LA[,b],2)$LA
}
#Estimate variance of beta_LA
var_beta_LA <- rowVars(betahatb_LA)
bias_beta_LA <- rowMeans(betahatb_LA) - betahat$LA

### b) 
#########Usikker her!
resb_LS <- matrix(NA,n_res,B)
for (b in 1:B){
  resb_LS[,b] <- ARp.resid(xb_LS[,b],betahatb_LS[,b])
}

resb_LA <- matrix(NA,n_res,B)
for (b in 1:B){
  resb_LA[,b] <- ARp.resid(xb_LA[,b],betahatb_LA[,b])
}

##USIKKER på res her!
x101_LS <- betahat$LS[1]*x0[n] + betahat$LS[2]*x0[n-1] + as.vector(resb_LS)
x101_LA <- betahat$LA[1]*x0[n] + betahat$LA[2]*x0[n-1] + as.vector(resb_LA)
pred_LS <- c(quantile(x101_LS,0.025), quantile(x101_LS,.975))
pred_LA <- c(quantile(x101_LA,0.025), quantile(x101_LA,.975))
