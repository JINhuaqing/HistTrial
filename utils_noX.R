## Some utils functions to implement the simulation of 
## Step 2 -4 under a very simple linear model
## when not considering X

#rm(list=ls())# need followings pkgs
# mvtnorm: for mvnorm kernel function
library(mvtnorm)





# truncate function
trunFn <- function(xs, lam.tru){
  rvs <- xs
  rvs[abs(xs) <= lam.tru] = 0
  return(rvs)
  
}



# optimize mu_1 given sigma and phi1
optMu1.noX <- function(data, invsigma2, phi1, Theta1){
  #args
  #    x: current data point, 1xp
  #    data: the datase, n x (2+p): [Y, Z, Xs]
  #    invsigma2: 1/sigma2, sigma2: variance of prior of mu1
  #    phi1: the value of phi1, scalar
  #    Theta1: the [theta1(X1), ..., theta1(Xn)], 1 x n
  
  Zs <- data$Z
  
  Ws <- 1/2/phi1**2 + invsigma2/2
  Ms <- data$Y/2/phi1**2 + invsigma2*Theta1/2
  
  num <- sum(Zs*Ms)
  den <- sum(Zs*Ws)
  num/den
}



# optimize mu_0 given tau2 and phi0
optMu0.noX <- function(data, Tau2, phi0, Theta0){
  #args
  #    data: the dataset, n x (2+p): [Y, Z, Xs]
  #    Tau2: 
  #    phi0: the value of phi0, scalar
  #    Theta0: 
  #    H: cov mat, p x p
  
  sZs <- 1 - data$Z
  
  Ws <- 1/phi0**2 + Tau2
  Ms <- data$Y/phi0**2 + Tau2* Theta0
  
  num <- sum(sZs*Ms)
  den <- sum(sZs*Ws)
  num/den
}


# optimize tau2 given mu0
optTau.noX <- function(data, Mu0, lam, Theta0, invgam2=0, lam.tru=0){
  #args
  #    data: the datase, n x (2+p): [Y, Z, Xs]
  #    Mu0
  #    lam: penalty parameter
  #    lam.tru: threshold to truncate the small value of Tau
  #    Theta0
  #    invgam2: 1/invgam2 is the variance parameter of HN prior
  sZs <- 1 - data$Z
  
  num <- sum(sZs)
  den <- sum(sZs * (Mu0-Theta0)**2) + invgam2
  tau2.td <- num/den
  
  tau2 <- trunFn(tau2.td, lam.tru)
  tau2 <- euclidean_proj_l1ball(rep(1, length(sZs))*tau2, lam*log(length(sZs)))
  return(list(tau2=tau2[1], tau2.td=tau2.td))
}


# can be used for ph0
optPhi0.noX <- function(data, Mu){
  #args
  #    data: the datase, n x (2+p): [Y, Z, Xs]
  #    Mu: 
  
  a <- 0.01
  b <-  0.01
  sZs <- 1 - data$Z
  
  num <- sum(sZs * (data$Y - Mu)**2)/2 + b
  den <- 1 + a + sum(sZs)/2
  res <- sqrt(num/den)
  res
}

# can be used for phi1
optPhi1.noX <- function(data, Mu){
  #args
  #    data: the datase, n x (2+p): [Y, Z, Xs]
  #    Mu:
  
  a <- 0.01
  b <-  0.01
  Zs <- data$Z
  
  num <- sum(Zs * (data$Y - Mu)**2)/2 + b
  den <- 1 + a + sum(Zs)/2
  res <- sqrt(num/den)
  res
}

# function to evaluate estimated mu0
mu0.noX.efn <- function(res){
  rv <- optMu0.noX(res$data, res$tau2, res$phi0, res$theta0)
  rv
}

mu1.noX.efn <- function(res){
  rv <- optMu1.noX(res$data, res$invsigma2, res$phi1, res$theta1)
  rv
}

# function to evaluate estimated tau2(X)
tau2.noX.efn <- function(res){
  rv <- optTau.noX(res$data, res$mu0, res$lam, res$theta0, res$invgam2, res$lam.tru)$tau2
  rv
}

# calculate the posterior variance of the mu0
post.var.mu0.noX.fn <- function(res){
  
  tau2 <- res$tau2
  sZs <- 1 - res$data$Z
  rv <- sum(sZs*(1/res$phi0/res$phi0+tau2))
  return(1/rv)
}

# calculate the posterior mean of the mu0(X)
post.mean.mu0.noX.fn <- function(res){
  tau2 <- res$tau2
  sZs <- 1 - res$data$Z
  Y <- res$data$Y
  theta0 <- res$theta0
  
  den <- sum(sZs*(1/res$phi0/res$phi0+tau2))
  num <- sum(sZs*(Y/res$phi0/res$phi0+tau2*theta0))
  rv <- num/den
  
  return(rv)
}

# function for estimate mu0s, tau2s, and phi0 under info and reference model
mu0.info.est.noX.fn <- function(Theta0, data, lam, phi0=NA, invgam2=0, is.ref=FALSE, maxit=100, lam.tru=0, lastRes=NA){
  #args:
  #  Theta0: Estimated ys based on historical data
  #  data: the dataset, n x (2+p): [Y, Z, Xs]
  #  lam: The penalty parameters
  #  is.ref: if true, reture the results of reference model distribution
  #  maxit: Maximal number of times for iteration
  #  LastRes: When not NA, initial value by results of last run
  
  phi0.tk <- c()
  Mu0.tk <- list()
  Tau2.tk <- list()
  p <- dim(data)[2] - 2
  Xs <- data[, 3:(p+2), drop=FALSE]
  n <- dim(Xs)[1]
  if (is.na(lastRes)[[1]]){
    Tau2 <- 0
  }else{
    Tau2 <- lastRes$tau2
  }
  
  if (!is.ref| is.na(phi0)){
    if (is.na(lastRes)[[1]]){
      phi0 <- 1
    }else{
      phi0 <- lastRes$phi0
    }
    for (i in 1:maxit){
      t0 <- Sys.time()
      Mu0 <- optMu0.noX(data, Tau2, phi0, Theta0)
      Mu0.tk[[i]] <- Mu0
      
      t1 <- Sys.time()
      phi0 <- optPhi0.noX(data, Mu0)
      phi0.tk[i] <- phi0
      
      t2 <- Sys.time()
      Tau2.res  <- optTau.noX(data, Mu0, lam, Theta0, invgam2, lam.tru)
      Tau2 <- Tau2.res$tau2
      Tau2.tk[[i]] <- Tau2.res$tau2
      t3 <- Sys.time()
      ts <- c(t0, t1, t2, t3)
      #print(diff(ts))
      
      if (i > 1){
        err.tau2 <- mean((Tau2.tk[[i]] - Tau2.tk[[i-1]])**2)
        err.mu <- mean((Mu0.tk[[i]] - Mu0.tk[[i-1]])**2)
        err.phi0 <- mean((phi0.tk[[i]] - phi0.tk[[i-1]])**2)
        
        err.all <- max(c(err.tau2, err.mu, err.phi0))
        if (err.all < 1e-5){
          break
        }
      }
    }
  }
  
  if (!is.ref){
    res <- list(tau2=Tau2, mu0=Mu0, phi0=phi0, 
                theta0=Theta0, data=data, lam=lam, 
                invgam2=invgam2, lam.tru=lam.tru, tau2.td=Tau2.res$tau2.td, iterN=i)
  }else{
    Tau20 <- 0
    Mu0.no <- optMu0.noX(data, Tau20, phi0, Theta0)
    res <- list(tau2=Tau20, mu0=Mu0.no, phi0=phi0, theta0=Theta0, data=data, 
                lam=lam, invgam2=invgam2, lam.tru=lam.tru, tau2.td=Tau20)
  }
  
  res
}



# function for estimate mu1s, phi when non-borrowing
mu1.no.est.noX.fn <- function(data){
  #args:
  #  data: the dataset, n x (2+p): [Y, Z, Xs]
  
  p <- dim(data)[2] - 2
  Xs <- data[, 3:(p+2), drop=F]
  n <- dim(Xs)[1]
  phi1 <- 1
  
  Mu1 <- optMu1.noX(data, 0, phi1, rep(0, n))
  phi1 <- optPhi1.noX(data, Mu1)
  
  
  res <- list(mu1=Mu1, phi1=phi1, invsigma2=0, data=data, theta1=rep(0, n))
}

# calculate the posterior variance of the mu1
post.var.mu1.noX.fn <- function(res){
  
  Zs <- res$data$Z
  rv <- sum(Zs*(1/res$phi1/res$phi1+res$invsigma2))
  return(1/rv)
}

# calculate the posterior mean of the mu1(X)
post.mean.mu1.noX.fn <- function(res){
  Zs <- res$data$Z
  Y <- res$data$Y
  theta1 <- res$theta1
  
  den <- sum(Zs*(1/res$phi1/res$phi1+res$invsigma2))
  num <- sum(Zs*(Y/res$phi1/res$phi1+res$invsigma2*theta1))
  rv <- num/den
  
  return(rv)
}


# sampling from posterior distribution of mu0
r.postMu0.noX <- function(res, M){
  # M: num of repitions
  m <- post.mean.mu0.noX.fn(res)
  v <- post.var.mu0.noX.fn(res)
  sps <- rnorm(M, m, sqrt(v))
  list(sps=sps, trts=sps, m=m, v=v)
}

# sampling from posterior distribution of mu1
r.postMu1.noX <- function(res, M){
  # M: num of repitions
  m <- post.mean.mu1.noX.fn(res)
  v <- post.var.mu1.noX.fn(res)
  sps <- rnorm(M, m, sqrt(v))
  list(sps=sps, trts=sps, m=m, v=v)
}


# Reweighted Pocock-Simon design without X
RPS.noX.design <- function(Z, R=1){
  #args:
  # Z: Assignment indicator, n, vec
  #Return:
  # probs: assign probs
  # grp: assigned group
  
  # when grp=1 for n+1 subjects
  Z1 <- c(Z, 1)
  n11 <- sum(Z1)
  n10 <- R*sum(1-Z1)
  gn1 <- sum((n11-n10)**2)/2
  
  # when grp=0 for n+1 subjects
  Z0 <- c(Z, 0)
  n01 <- sum(Z0)
  n00 <- R*sum(1-Z0)
  gn0 <- sum((n01-n00)**2)/2
  
  gs.raw <- c(gn0, gn1)
  gs <- gs.raw/sum(gs.raw)
  
  if(sum(gs==0) == 1){
    probs <- c(gs[2], gs[1])
  }else{
    probs <- allo.probs.fn(gs)
  }
  
  
  grp <- sample.int(2, size=1, prob=probs)
  rv <- list(grp=grp, probs=probs, gs=gs.raw)
  rv
}

# Reweighted biased coin design
RBC.noX.design <- function(Z, R=1){
  #args:
  # Z: Assignment indicator, n, vec
  #Return:
  # ns: num of subs in each grp
  # probs: assign probs
  # grp: assigned group
  
  
  n0 <- R* sum((1-Z))
  n1 <- sum(Z)
  
  ns <- c(n0, n1)
  gs <- ns/sum(ns)
  if (sum(ns)==0){
    probs <- c(0.5, 0.5)
  }else if(sum(gs==0) == 1){
    probs <- c(gs[2], gs[1])
  }else{
    probs <- allo.probs.fn(gs)
  }
  grp <- sample.int(2, size=1, prob=probs)
  rv <- list(grp=grp, probs=probs, ns=ns)
  rv
}

