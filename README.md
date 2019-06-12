# Introduction
The scripts implement an EM algorithm to fit the mixture model in Li and Shen (2019) to assess the pleiotropic effects of the CCR5delta32 mutation.

# EM algorithm for the mixture model
The following three functions implements the EM algorithm for fitting the 3-component mixture model: 
<img src="http://bit.ly/31cW8QX" align="center" border="0" alt="\pi_{-}N(-\mu_Z,1)+\pi_{0}N(0,1)+\pi_{+}N(\mu_Z,1)" width="311" height="19" /> 
for the Z scores of the CCR5delta32 variant across a series of complex diseases.

``` r
## mixture model 190603 three components

compute.log.lik <- function(X, L, w) {
  L[,1] <- L[,1]*w[1]
  L[,2] <- L[,2]*w[2]
  L[,3] <- L[,3]*w[3]
  return(sum(log(rowSums(L))))
}

mixture.EM <- function(w.init, mu.init, X) {
  w.curr <- w.init
  mu <- mu.init
  L = matrix(NA, nrow = length(X), ncol= 3)
	L[, 1] <- dnorm(X, mean = 0, sd = 1)
	L[, 2] <- dnorm(X, mean = mu, sd = 1)
	L[, 3] <- dnorm(X, mean = -mu, sd = 1)
  L.curr <- L

  # store log-likehoods for each iteration
  log_liks <- c()
  ll       <- compute.log.lik(X, L.curr, w.curr)
  log_liks <- c(log_liks, ll)
  delta.ll <- 1
  
  while(delta.ll > 1e-5) {
    em <- EM.iter(w.curr, L.curr, X)
    mu.curr   <- em$mu.next
    w.curr   <- em$w.next
    L.curr[, 2] <- dnorm(X, mean = mu.curr, sd = 1)
    L.curr[, 3] <- dnorm(X, mean = -mu.curr, sd = 1)
    ll       <- compute.log.lik(X, L.curr, w.curr)
    log_liks <- c(log_liks, ll)
    delta.ll <- log_liks[length(log_liks)]  - log_liks[length(log_liks)-1]
  }
  return(list(pi = w.curr, mu = mu.curr, logliks = log_liks))
}

EM.iter <- function(w.curr, L, X) {
  # E-step
  z_ik <- L
  for(i in 1:ncol(L)) {
    z_ik[,i] <- w.curr[i]*z_ik[,i]
  }
  z_ik <- z_ik / rowSums(z_ik)
  # M-step
  mu.next <- sum(z_ik[,2]*X)/sum(z_ik[,2])
  w.next   <- colSums(z_ik)/sum(z_ik)
  return(list(mu.next = mu.next, w.next = w.next)) 
}
```

# The analysis of CCR5delta32 in UK Biobank
Below, we set up starting values to estimate the parameters for the paper.

``` r
## real data: X is a vector containing 131 Z scores for the associations 
## between the mutation and 131 curated disease phenotypes

estimation <- mixture.EM(w.init = c(.3, .4, .3), mu.init = 1, X)

estimation
```


``` r
## $pi
## [1] 0.00011604534 0.76893798357 0.23094597109

## $mu
## [1] 1.0026908

## $logliks
## ...
```

# Bootstrapping the standard errors
The following code performs bootstrap resampling to estimate the standard errors.

``` r
myfun <- function(X, indices) {
  X <- as.numeric(as.matrix(X[indices,]))
  w.curr <- c(.3, .4, .3)
  mu <- 1
  L = matrix(NA, nrow = length(X), ncol= 3)
	L[, 1] = dnorm(X, mean = 0, sd = 1)
	L[, 2] = dnorm(X, mean = mu, sd = 1)
	L[, 3] = dnorm(X, mean = -mu, sd = 1)
  L.curr <- L
  
  # store log-likehoods for each iteration
  log_liks <- c()
  ll       <- compute.log.lik(X, L.curr, w.curr)
  log_liks <- c(log_liks, ll)
  delta.ll <- 1
  
  while(delta.ll > 1e-5) {
    em <- EM.iter(w.curr, L.curr, X)
    mu.curr   <- em$mu.next
    w.curr   <- em$w.next
    L.curr[, 2] = dnorm(X, mean = mu.curr, sd = 1)
    L.curr[, 3] = dnorm(X, mean = -mu.curr, sd = 1)
    ll       <- compute.log.lik(X, L.curr, w.curr)
    log_liks <- c(log_liks, ll)
    delta.ll <- log_liks[length(log_liks)]  - log_liks[length(log_liks)-1]
  }
  return(c(w.curr, mu.curr))
}


data <- data.frame(x = X)
require(boot)
boot.result <- boot(data = data, statistic = myfun, R = 99)

boot.result
```

``` r
## ORDINARY NONPARAMETRIC BOOTSTRAP


## Call:
## boot(data = data, statistic = myfun, R = 99)


## Bootstrap Statistics :
##          original       bias    std. error
## t1* 0.00011604534  0.090873373 0.129549705
## t2* 0.76893798357 -0.070296231 0.124585479
## t3* 0.23094597109 -0.020577142 0.045616628
## t4* 1.00269083189  0.069282159 0.139459538
```




