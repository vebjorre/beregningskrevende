library(boot)
library(MASS)
date <- coal$date
t0 <- date[1]
t2 <- date[191]
date <- date[-c(1,191)]
plot(date, 1:189)

target <- function(x){
  return (1/x[4]^5 * exp(-(1+x[2]+x[3])/x[4]) * x[2]^(x[5]+1) * x[3]^(189-x[5]+1) * exp(x[1]*(x[3]-x[2]) + x[2]*t0 - x[3]*t2))
}

fullcond <- function(x,r){
  if (r==1){
    return (log(dexp(x[1], rate=1/(x[1]-x[2]))))
  }
  else if (r==2){
    return (log((dgamma(x[2], shape=x[5]+2, scale=1/(1/x[4]+x[1]-t0)))))
  }
  else if (r==3){
    return (log(dgamma(x[3], shape=189-x[5]+2, scale=1/(1/x[4]+t2-x[1]))))
  }
  else if (r==4){
    return (log(1/x[4]^5 * exp(-1/x[4]*(1+x[2]+x[3]))))
  }
}

fullcond.b1 <- function(x){
  return (-5*log(x[4]) + log( x[2]^(x[5]+1) * x[3]^(189-x[5]+1) * exp(-(1+x[2]+x[3])/x[4] + x[1]*(x[3]-x[2]) + x[2]*t0 - x[3]*t2) ))
}

mcmc_RW <- function(d, ntimes = 1000)
{
  x <- matrix(nrow=ntimes, ncol=5)
  x[1,1] <- runif(1,t0,t2)
  x[1,2:4] <- runif(3,0,1)
  x[1,5] <- sum(date<x[1,1])
  for(i in 2:ntimes)
  {
    for (j in 1:4){
      y <- max(0,runif(1,x[i-1,j]-d[j],x[i-1,j]+d[j]))
      temp <- x[i-1,]
      temp[j] <- y
      alpha <- min(1,exp(fullcond(temp,j)- fullcond(x[i-1,],j)), na.rm=TRUE)
      # cat("alpha:",alpha, "j:",j, "y:",y, "x:",x[i-1,], "\n")
      if (runif(1)<alpha){
        if (j==1){
          if (t0 <= y && y <= t2){
            x[i,j] = y
          }
          else{
            x[i,j] = x[i-1,j]
          }
        }
        else{
          x[i,j] <- y
        }
        
      }
      else{
        x[i,j] <- x[i-1,j]
      }
    }
    x[i,5] <- sum(date<=x[i,1])
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
    t1 <- rnorm(1,mean=x[i-1,1],sd=.5)
    y0 <- sum(date<t1)
    lam0 <- max(0,runif(1, x[i-1,2]-d[2], x[i-1,2]+d[2]))
    lam1 <- max(0,runif(1, x[i-1,3]-d[3], x[i-1,3]+d[3]))
    temp <- x[i-1,]
    temp[1] <- t1
    temp[2] <- lam0
    temp[3] <- lam1
    temp[5] <- y0
    alpha <- min(1,exp(fullcond.b1(temp)-fullcond.b1(x[i-1,])), na.rm=TRUE)
    # cat("alpha:",alpha, "x:",x[i-1,], "\n")
    if (runif(1)<alpha){
      x[i,] <- temp
    }
    else{
      x[i,] <- x[i-1,]
    }
    
    temp <- x[i,]
    beta <- max(0,rnorm(1,mean=temp[4], sd=.5))
    temp[4] <- beta
    lam0 <- max(0,runif(1, x[i,2]-d[2], x[i,2]+d[2]))
    lam1 <- max(0,runif(1, x[i,3]-d[3], x[i,3]+d[3]))
    temp[2] <- lam0
    temp[3] <- lam1
    alpha <- min(1, exp(fullcond.b1(temp)-fullcond.b1(x[i,])), na.rm=TRUE)
    # cat("alpha:",alpha, "x:",x[i-1,], "\n")
    if (runif(1)<alpha){
      x[i,] <- temp
    }
  }
  return(x)
}

# mcmc_block2 <- function(d, ntimes)
# {
#   x <- matrix(nrow=ntimes, ncol=5)
#   x[1,1] <- runif(1,t0,t2)
#   x[1,2:4] <- runif(3,0,1)
#   x[1,5] <- sum(date<x[1,1])
#   for(i in 2:ntimes)
#   {
#     temp <- x[i-1,]
#     temp[4] <- max(0,rnorm(1,mean=x[i-1,4], sd=0.1))
#     lam0 <- rgamma(1,shape=x[i-1,5]+2, rate=1/(1/temp[4]+x[i-1,1]-t0))
#     lam1 <- rgamma(1,shape=189-x[i-1,5]+2, rate=1/(1/temp[4]+t2-x[i-1,1]))
#     temp[2] <- lam0
#     temp[3] <- lam1
#     alpha <- min(1,target(temp)/target(x[i-1,]))
#     cat("alpha:",alpha, "temp:", temp, "x:",x[i-1,], "\n")
#     if (runif(1)<alpha){
#       x[i,] <- temp
#     }
#     else{
#       x[i,] <- x[i-1,]
#     }
#     temp <- x[i,]
#     temp[1] <- runif(1,temp[1]-d[1], temp[1]+d[1])
#     alpha <- min(1,target(temp)/target(x[i,]))
#     if (runif(1)<alpha){
#       if (t0 <= temp[1] && temp[1] <= t2){
#         x[i,] = temp
#       }
#     }
#     x[i,5] <- sum(date < x[i,1])
#   }
#   return(x)
# }


d <- c(5,.1,.1,.1)
d.b <- c(1,.1,.1,.1)
ntimes <- 50000


x.b <- mcmc_block(d.b,ntimes)

par(mfrow=c(2,2))
plot(x.b[,1], type='l')
plot(x.b[,2], type='l')
plot(x.b[,3], type='l')
plot(x.b[,4], type='l')

par(mfrow=c(2,2))
truehist(x.b[-c(1:1000),1])
truehist(x.b[-c(1:1000),2])
truehist(x.b[-c(1:1000),3])
truehist(x.b[-c(1:1000),4])


# x1 = mcmc_RW(d, ntimes)
# 
# par(mfrow=c(2,2))
# plot(x1[,1], type='l')
# plot(x1[,2], type='l')
# plot(x1[,3], type='l')
# plot(x1[,4], type='l')
# 
# par(mfrow=c(2,2))
# truehist(x1[-c(1:1000),1])
# truehist(x1[-c(1:1000),2])
# truehist(x1[-c(1:1000),3])
# truehist(x1[-c(1:1000),4])

# par(mfrow=c(2,2))
# acf(x1[,1])
# acf(x1[,2])
# acf(x1[,3])
# acf(x1[,4])
# 
# mean(x1[,1])
# mean(x1[,2])
# mean(x1[,3])
# mean(x1[,4])
