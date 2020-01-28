library(Matrix)

random.exp <- function(n,lam){
  return( - 1 / lam * log(runif(n)))
}

random.g <- function(n,a){
  c <- a*exp(1)/(a+exp(1))
  u <- runif(n)
  x <- rep(0,n)
  for (i in 1:n){
    if (u[i]<c/a){
      x[i] <- (a*u[i]/c)^(1/a)
    }
    else{
      x[i] <- (-log(1/a + exp(-1) - u[i]/c))
    }
  }
  return(x)
}

random.norm <- function(n){
  x1 <- runif(n)*2*pi
  x2 <- random.exp(n,.5)
  return(sqrt(x2)*cos(x1))
}

random.multinorm <- function(n,mu,S){
  x <- random.norm(n)
  A <- chol(S)
  return(mu + A%*%x)
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
