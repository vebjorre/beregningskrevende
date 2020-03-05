library(boot)
library(MASS)
library(invgamma)
date <- coal$date #Fetch data
t0 <- date[1] #Start date
t2 <- date[191] #End date
date <- date[-c(1,191)] #Data
plot(date, 1:189)

logtarget <- function(x){ # x1=t1, x2=lam0, x3=lam1, x4=beta, x5=y0
  return (-5*log(x[4]) + (x[5]+1)*log(x[2]) + (189-x[5]+1)*log(x[3]) - (1+x[2]+x[3])/x[4] + x[1]*(x[3]-x[2]) + x[2]*t0 - x[3]*t2 )
}

fullcond <- function(x,r){ #x1=t1, x2=lam0, x3=lam1, x4=beta, x5=y0
  if (r==2){
    return (rgamma(1, shape=x[5]+2, rate=1/x[4]+x[1]-t0))
  }
  else if (r==3){
    return (rgamma(1, shape=189-x[5]+2, rate=1/x[4]+t2-x[1]))
  }
  else if (r==4){
    return (rinvgamma(1, shape=4, rate=1+x[2]+x[3]))
  }
}



mcmc_RW <- function(d, ntimes = 1000){ #x1=t1, x2=lam0, x3=lam1, x4=beta, x5=y0
  x <- matrix(nrow=ntimes, ncol=5)
  x[1,1] <- 1900
  x[1,2:4] <- runif(3,0,1)
  x[1,5] <- sum(date<x[1,1])
  for(i in 2:ntimes){
    x[i,2] <- fullcond(x[i-1,],2)
    x[i,3] <- fullcond(x[i-1,],3)
    x[i,4] <- fullcond(x[i-1,],4)
    temp <- x[i-1,]
    t1 <- rnorm(1, temp[1], d)
    while (t1<t0 || t1>t2){
      t1 <- rnorm(1, temp[1], d)
    }
    temp[1] <- t1
    temp[5] <- sum(date<t1)
    alpha <- min(1,exp(logtarget(temp)- logtarget(x[i-1,])), na.rm=TRUE)
    # cat("alpha:",alpha, "t1:",t1, "x:",x[i-1,], "\n")
    if (runif(1)<alpha && t0<t1 && t1<t2){
      x[i,1] <- t1
      x[i,5] <- sum(date<=t1)
      }
    else{
      x[i,1] <- x[i-1,1]
      x[i,5] <- x[i-1,5]
      }
  }
  return(x)
}

mcmc_block <- function(d, ntimes)
{
  x <- matrix(nrow=ntimes, ncol=5)
  x[1,1] <- runif(1,t0,t2)
  x[1,2:4] <- runif(3,0,1)
  x[1,5] <- sum(date<x[1,1])
  for(i in 2:ntimes)
  {
    t1 <- rnorm(1,mean=x[i-1,1],sd=d)
    while (t0 > t1 || t1 > t2){
      t1 <- rnorm(1,mean=x[i-1,1],sd=d)
    }
    y0 <- sum(date<t1)
    temp <- x[i-1,]
    temp[1] <- t1
    temp[5] <- y0
    lam0 <- fullcond(temp,2)
    lam1 <- fullcond(temp,3)
    temp[2] <- lam0
    temp[3] <- lam1
    alpha <- min(1,exp(logtarget(temp)-logtarget(x[i-1,])), na.rm=TRUE)
    if (runif(1)<alpha){
      x[i,] <- temp
    }
    else{
      x[i,] <- x[i-1,]
    }
    
    temp <- x[i,]
    beta <- rnorm(1,mean=temp[4], sd=d)
    while (beta <= 0){
      beta <- rnorm(1,mean=temp[4], sd=d)
    }
    temp[4] <- beta
    lam0 <- fullcond(temp,2)
    lam1 <- fullcond(temp,3)
    temp[2] <- lam0
    temp[3] <- lam1
    alpha <- min(1,exp(logtarget(temp)-logtarget(x[i-1,])),na.rm=TRUE)
    if (runif(1) < alpha){
      x[i,] <- temp
    }
  }
  return(x)
}

d <- 4
d.b <- .5
ntimes <- 50000



x1 = mcmc_RW(d, ntimes)

par(mfrow=c(2,2))
plot(x1[,1], type='l')
plot(x1[,2], type='l')
plot(x1[,3], type='l')
plot(x1[,4], type='l')

par(mfrow=c(2,2))
truehist(x1[-c(1:1000),1])
truehist(x1[-c(1:1000),2])
truehist(x1[-c(1:1000),3])
truehist(x1[-c(1:1000),4], xlim=c(0,10))

par(mfrow=c(2,2))
acf(x1[,1])
acf(x1[,2])
acf(x1[,3])
acf(x1[,4])

mean(x1[-c(1:1000),1])
mean(x1[-c(1:1000),2])
mean(x1[-c(1:1000),3])
mean(x1[-c(1:1000),4])


# x.b <- mcmc_block(d.b,ntimes)
# 
# par(mfrow=c(2,2))
# plot(x.b[,1], type='l')
# plot(x.b[,2], type='l')
# plot(x.b[,3], type='l')
# plot(x.b[,4], type='l')
# 
# par(mfrow=c(2,2))
# truehist(x.b[-c(1:1000),1])
# truehist(x.b[-c(1:1000),2])
# truehist(x.b[-c(1:1000),3])
# truehist(x.b[-c(1:1000),4])

