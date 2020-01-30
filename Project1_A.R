library(Matrix)

random.exp <- function(n,lam){
  return( - 1 / lam * log(runif(n)))
}

random.g <- function(n,alpha){
  c <- alpha*exp(1)/(alpha+exp(1))
  u <- runif(n)
  x <- rep(0,n)
  for (i in 1:n){
    if (u[i]<c/alpha){
      x[i] <- (alpha*u[i]/c)^(1/alpha)
    }
    else{
      x[i] <- (-log(1/alpha + exp(-1) - u[i]/c))
    }
  }
  return(x)
}

random.norm <- function(n){
  x1 <- runif(n)*2*pi
  x2 <- random.exp(n,.5)
  return(sqrt(x2)*cos(x1))
}

random.multinorm <- function(n,mu,Sigma){
  x <- random.norm(n)
  A <- chol(Sigma)
  return(mu + t(A)%*%x)
}

random.dirichlet <- function(alpha){
  K <- length(alpha)
  z <- rgamma(K,alpha,1) #n,shape,rate
  v <- sum(z)
  return (z/v)
}

# testing g:
a <- .6
c <- a*exp(1)/(a+exp(1))
test <- random.g(10000, a)
e.mean <- mean(test)
e.var <- var(test)
a.mean <- c/(a+1)+2*c*exp(-1)
a.var <- c*(1/(a+2)+ 5*exp(-1)) - a.mean^2

#testing standard normal:
normal.sample <- random.norm(10000)
mean(normal.sample)
var(normal.sample)

#Testing multinormal
d = 2
N = 10000
sigma <- cbind(c(2,1), c(1,2))
mu <- c(2,1)
multinormal.sample <- matrix(NA,N,d)
for (i in 1:N){
  multinormal.sample[i,] <- random.multinorm(d, mu, sigma)
}
colMeans(multinormal.sample)
var(multinormal.sample)
sigma
hist(multinormal.sample[,1])
hist(multinormal.sample[,2])

K <- 10
N_dir <- 1000
alpha <- rep(0.5,K)
dirichlet.sample <- matrix(NA,N_dir,K)
for (i in 1:N_dir){
  dirichlet.sample[i,] <- random.dirichlet(alpha)
  
}
mean(dirichlet.sample)
var(dirichlet.sample)
alphasum <- sum(alpha)
alpha1 <- alpha/alphasum
alpha1*(1-alpha)/(1+alphasum)

f <- function(theta){
  return ((2+theta)^125*(1-theta)^(18+20)*theta^34)
}
  
logf <- function(theta){
  return (log((2+theta)^125*(1-theta)^(18+20)*theta^34))
}

cf <- function(theta){
  return ((2+theta)^125*(1-theta)^(18+20)*theta^34) / f(15/394+sqrt(53809)/394)
}

theta_new_f <- function(theta){
  return (theta*new_f(theta))
}

new_f <- function(theta){
  return (((2+theta)^125*(1-theta)^(18+20)*theta^34)/as.numeric(integrate(f,0,1)[1]))
}


random.multinomial <- function(n){
  x.out <- rep(NA,n)
  rejections <- 0
  for (i in 1:n){
    finished <- FALSE
    c <- f(15/394+sqrt(53809)/394)
    while(!finished){
      u <- runif(1)
      x <- runif(1)
      alpha <- f(x) / c
      if (u <= alpha){
        finished <- TRUE
      }
      else{
        rejections = rejections + 1
      }
    }
    x.out[i] <- x
  }
  cat(rejections,"rejections out of",(n+rejections),"trials.(",(rejections)/(n+rejections)*100, "% )")
  return (x.out)
}

importance_sampling <- function(n){
  x <- random.multinomial(n)
  mu <- sum(x*(1-x)^4/f(x))/sum((1-x)^4/f(x))
  return(mu)
}

c <- f(15/394+sqrt(53809)/394)

test <- random.multinomial(10000)
mean(test)
var(test)
hist(test, prob=TRUE)
xlin <- seq(0,1,.01)
lines(new_f(xlin))
curve(new_f,0,1)


newc <- integrate(f, 0, 1)
integrate(theta_new_f,0,1)


test2 <- importance_sampling(10000)
test2
