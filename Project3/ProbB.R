bilirubin <- read.table("bilirubin.txt",header=T)

boxplot(log(meas)~pers,data=bilirubin)

fit <- lm(log(meas)~pers,data=bilirubin)
Fval <- summary(fit)$fstatistic[1]

permTest <- function(x0,y0){
  n <- length(x0)
  x <- sample(x0,n)
  return(as.numeric(summary(lm(y0~x))$fstatistic[1]))
}

n_perms <- 999
F_sample <- rep(NA,n_perms)
for (i in 1:n_perms){
  F_sample[i] <- permTest(bilirubin$pers,log(bilirubin$meas))
}
pval <- sum(F_sample > Fval)/n_perms
