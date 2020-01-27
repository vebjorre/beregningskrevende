

random.exp <- function(n,lam){
  return( - 1 / lam * log(runif(n)))
}

G <- function(x, a){ 
  c = 1/(1/a + exp(-1) - exp(-1000000))
  #usikker ^
  return (c*(1/a + exp(-1) - exp(-x)))
}

random.g <- function(x,a)

x = 1:5000/1000
a = .5
plot(x,G(x,a))

     