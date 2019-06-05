# Introduction
The scripts implement an EM algorithm to fit the mixture model in Li and Shen (2019) to assess the pleiotropic effects of the CCR5Delta32 mutation.

# EM algorithm for the mixture model
The following three functions implements the EM algorithm for fitting the 3-component mixture model: $\pi_{-}N(-\mu_Z,1)+\pi_{0}N(0,1)+\pi_{+}N(\mu_Z,1)$ for the Z scores of the CCR5$\Delta$32 variant across a series of complex diseases.

``` r
## mixture model 190603 three components

compute.log.lik <- function(X, L, w) {
  L[,1] = L[,1]*w[1]
  L[,2] = L[,2]*w[2]
  L[,3] = L[,3]*w[3]
  return(sum(log(rowSums(L))))
}

mixture.EM <- function(w.init, mu.init, X) {
  
  w.curr <- w.init
  mu <- mu.init
  L = matrix(NA, nrow = length(X), ncol= 3)
	L[, 1] = dnorm(X, mean = 0, sd = 1)
	L[, 2] = dnorm(X, mean = mu, sd = 1)
	L[, 3] = dnorm(X, mean = -mu, sd = 1)
  L.curr <- L

  #print(summary(L.curr))
  #print(dim(L.curr))
  # store log-likehoods for each iteration
  log_liks <- c()
  ll       <- compute.log.lik(X, L.curr, w.curr)
  log_liks <- c(log_liks, ll)
  delta.ll <- 1
  
  while(delta.ll > 1e-5) {
    em <- EM.iter(w.curr, L.curr, X)
    mu.curr   <- em$mu.next
    #sigma0.curr   <- em$sigma0.next
    #sigma1.curr   <- em$sigma1.next
    w.curr   <- em$w.next
    #print(mu.curr)
    #print(w.curr)
    L.curr[, 2] = dnorm(X, mean = mu.curr, sd = 1)
    L.curr[, 3] = dnorm(X, mean = -mu.curr, sd = 1)
    ll       <- compute.log.lik(X, L.curr, w.curr)
    log_liks <- c(log_liks, ll)
    delta.ll <- log_liks[length(log_liks)]  - log_liks[length(log_liks)-1]
  }
  return(list(pi = w.curr, mu = mu.curr, logliks = log_liks))
}

EM.iter <- function(w.curr, L, X) {
  
  # E-step: compute E_{Z|X,w0}[I(Z_i = k)]
  z_ik <- L
  for(i in 1:ncol(L)) {
    z_ik[,i] <- w.curr[i]*z_ik[,i]
  }
  z_ik     <- z_ik / rowSums(z_ik)
 
 #print(dim(z_ik))
  #print(length(X))
  # M-step
  mu.next <- sum(z_ik[,2]*X)/sum(z_ik[,2])
  #sigma0.next <- sum(z_ik[,1]*(X - mu.next)**2)/sum(z_ik[,1])
  #sigma1.next <- sum(z_ik[,2]*(X - mu.next)**2)/sum(z_ik[,2])
  w.next   <- colSums(z_ik)/sum(z_ik)
  return(list(mu.next = mu.next, w.next = w.next)) #sigma0.next, sigma1.next
}
```

