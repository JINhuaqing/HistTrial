## Some utils functions to implement the simulation of 
## Step 2 -4 under a very simple linear model
## Y_i = b I(Z_i) + betas X_i + eps
#rm(list=ls())# need followings pkgs
# rootSolve: to solve the root
# mvtnorm: for mvnorm kernel function
library(mvtnorm)
library(arrApply)
library(kedd)




# generate the covariates
gen.Data.Xs <- function(n, X.tps){
  # args
  #   n: Sample size for both groups
  #   X.tps: A vector to contain the variable types of beta
  #           "c": continous variable
  #           number: number of level for discrete variable
  Xs.list <- list()
  
  p <- length(X.tps)
  i <- 1
  for (X.tp in X.tps){
    if (substr(tolower(X.tp), 1, 1) == "c"){
        cX <- rnorm(n, i/2, sd=3)
        cX <-  2*exp(cX)/(1+exp(cX)) - 1
    }else {
        nlvl <- as.numeric(X.tp)
        cX <- sample(nlvl, size=n, replace=TRUE) - 1
        # cX <- as.factor(cX)
    }
    Xs.list[[i]] <- cX
    i <- i + 1
  }
  Xs <- do.call(cbind, Xs.list)
  Xs <- as.data.frame(Xs)
  names(Xs) <- paste0("X", 1:p)
  Xs
}

# functions of beta parameters given  coviariates X
# only used for subgroup model
sub.Paras.fn <- function(Xs, betass){
    #args:
    #   Xs: covariates, n x p 
    idxs <- Xs$X1*2 + Xs$X2 + 1
    betMat <- do.call(rbind, betass[idxs])
    return(betMat)
    
}
# Generate the mean of current data given the Xs, Zs, and true parameters
# Vary depends on model
curMean.fn <- function(Xs, Zs, betMat, b){
  #args:
  #  Xs: covariates matrix, n x p, data frame
  #  Zs: Group indicate function, n, vector
  #  betMat: matrix of beta parameters, n x p, matrix
  #  b: treatment effect, scalar
  if (is.null(dim(Xs))){
        fXs <- c(1, as.vector(Xs))
  }else{
        n <- dim(Xs)[1]
        fXs <- cbind(rep(1, n), Xs)
        fXs <- as.matrix(fXs)
  }
  mys <- rowSums(fXs * betMat) + b * Zs
  mys
}

# Reweighted biased coin design
RBC.design <- function(xNew, Xs, Z, hs, R=1){
    #args:
    # xNew: new data point, p x 1, vec
    # Xs: Cumulated covariates, n x p, data frame
    # Z: Assignment indicator, n, vec
    # hs: band widths, p, vec
    #Return:
    # ns: num of subs in each grp
    # probs: assign probs
    # grp: assigned group
    
    p <- length(xNew)
    n <- dim(Xs)[1]
    if (missing(hs)){
        hs <- rep(2.1, p)
    }
    
    xNewMat <- matrix(rep(xNew, n), nrow=n, byrow=TRUE)
    diff <- as.matrix(Xs) - xNewMat
    wss <- normKfn(diff, hs)
    ws <- arrApply(wss, 2, "prod")
    
    n0 <- R* sum((1-Z)*ws)
    n1 <- sum(Z*ws)
    
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

# Reweighted Pocock-Simon design
RPS.design <- function(xNew, Xs, Z, hs, R=1){
    #args:
    # xNew: new data point, p x 1, vec
    # Xs: Cumulated covariates, n x p, data frame
    # Z: Assignment indicator, n, vec
    # hs: band widths, p, vec
    #Return:
    # probs: assign probs
    # grp: assigned group
    
    p <- length(xNew)
    n <- dim(Xs)[1]
    if (missing(hs)){
        hs <- rep(2.1, p)
    }
    xNewMat <- matrix(rep(xNew, n+1), nrow=n+1, byrow=TRUE)
    Xs.ext <- rbind(as.matrix(Xs), xNew)
    diff <- Xs.ext - xNewMat
    wss <- normKfn(diff, hs)
    
    # when grp=1 for n+1 subjects
    Z1 <- c(Z, 1)
    Z1Mat <- matrix(rep(Z1, p), ncol=p, byrow=FALSE)
    n11s <- colSums(Z1Mat * wss)
    n10s <- R*colSums((1-Z1Mat) * wss)
    gn1 <- sum((n11s-n10s)**2)/2
    
    # when grp=0 for n+1 subjects
    Z0 <- c(Z, 0)
    Z0Mat <- matrix(rep(Z0, p), ncol=p, byrow=FALSE)
    n01s <- colSums(Z0Mat * wss)
    n00s <- R*colSums((1-Z0Mat) * wss)
    gn0 <- sum((n01s-n00s)**2)/2
    
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



# Calculate the value of the kernel functions
MKF <- function(x, Xs, H){
  # args:
  #   x: current x, p x 1
  #   Xs: data matrix, n x p
  #   H: cov mat, p x p
  n <- dim(Xs)[1]
  xmat <- matrix(rep(x, n), nrow=n, byrow=TRUE)
  difX <- xmat - Xs
  vs <- dmvnorm(difX, sigma=H)
  vs
}


# kernel function to evaluate multiple points at one time
mMKF <- function(cxs, Xs, H){
  # args:
  #   cxs: current x, m x p 
  #   Xs: data matrix, n x p
  #   H: cov mat, p x p
  #return:
  #   vs: n x m matrix 
  if (is.null(dim(cxs))){
    vs <- MKF(cxs, Xs, H)
  }else{
    m <- dim(cxs)[1]
    p <- dim(Xs)[2]
    n <- dim(Xs)[1]
    cxArr <- array(rep(cxs, n), dim = c(m, p, n))
    XsArr <- replicate(m, Xs, simplify = "array")
    XsArr <- aperm(XsArr, c(3, 2, 1))
    difXArr <- cxArr - XsArr
    difXArr <- aperm(difXArr, c(3, 1, 2))
    difXmat <- matrix(difXArr, ncol=p)
    vs.vec <- dmvnorm(difXmat, sigma=H)
    vs <- matrix(vs.vec, ncol=m)
  }
  
  return(vs)
}

# soft threholding 
softH <- function(xs, lam){
    rvs <- xs
    rvs[abs(xs) <= lam] = 0
    pXs <- rvs[abs(xs) > lam]
    rvs[abs(xs) > lam] =  sign(pXs) * (abs(pXs)-lam)
    return(rvs)
}

# univariate kernel fn with K(0) = 1
normKfn <- function(x, h=1, kt="epanechnikov"){
  #args:
  # x: n x 1 or n x p
  # h: bw, of size p
  

  if (is.null(dim(x))){
    p <- 1
  }else{
    p <- dim(x)[2]
  }
  
  x <- matrix(x, ncol=p)
  
  if (length(h) == 1) {
    h <- rep(h, p)
  }
  
  c <- kernel.fun(0, kernel=kt)$kx
  ys <- lapply(1:p, function(i){kernel.fun(x[, i]/h[i], kernel=kt)$kx/h[i]/c})
  y <- do.call(cbind, ys)
  y
}

# calculate the allo prob given the imbalance measure
allo.probs.fn <- function(gs){
  invs <- 1/gs - 1
  probs <- invs/sum(invs)
  probs
  
}

# optimize mu_1(X) given sigma and phi1
optMu1 <- function(x, data, invsigma2, phi1, Theta1s, H){
  #args
  #    x: current data point, 1xp
  #    data: the datase, n x (2+p): [Y, Z, Xs]
  #    invsigma2: 1/sigma2
  #    phi1: the value of phi1, scalar
  #    Theta1s: the [theta1(X1), ..., theta1(Xn)], 1 x n
  #    H: cov mat, p x p
    
    p <-length(x)
    Xs <- as.matrix(data[, 3:(p+2)])
    Zs <- data$Z
    
    Ws <- 1/2/phi1**2 + invsigma2/2
    Ms <- data$Y/2/phi1**2 + invsigma2*Theta1s/2
    ks <- MKF(x, Xs, H=H)
    
    num <- sum(Zs*ks*Ms)
    den <- sum(Zs*ks*Ws)
    num/den
}

mOptMu1 <- function(cxs, data, invsigma2, phi1, Theta1s, H){
  #args
  #    cxs: current data point, mxp
  #    data: the datase, n x (2+p): [Y, Z, Xs]
  #    invsigma2: 1/sigma2
  #    phi1: the value of phi1, scalar
  #    Theta1s: the [theta1(X1), ..., theta1(Xn)], 1 x n
  #    H: cov mat, p x p
    
    if (is.null(dim(cxs))){
        vs <- optMu1(cxs, data, invsigma2, phi1, Theta1s, H)
    }else{
        m <- dim(cxs)[1]
        p <- dim(cxs)[2]
        Xs <- as.matrix(data[, 3:(p+2)])
        n <- dim(Xs)[1]
        Zs <- data$Z
        ZsMat <- matrix(rep(Zs, m), ncol=m)
        
        mks <- mMKF(cxs, Xs, H=H) # n x m 
        Ws <- 1/2/phi1**2 + invsigma2/2
        WsMat <- matrix(rep(Ws, m*n), ncol=m)
        Ms <- data$Y/2/phi1**2 + invsigma2* Theta1s/2
        MsMat <- matrix(rep(Ms, m), ncol=m)
        
        num <- colSums(ZsMat*mks*MsMat)
        den <- colSums(ZsMat*mks*WsMat)
        vs <- num/den
      
    }
    return(vs)
}


# optimize mu_0(X) given [tau2(X1), ..., tau2(Xn)] and phi0
optMu0 <- function(x, data, Tau2s, phi0, Theta0s, H){
  #args
  #    x: current data point, 1xp
  #    data: the dataset, n x (2+p): [Y, Z, Xs]
  #    Tau2s: the [tau2(X1), ..., tau2(Xn)], 1 x n
  #    phi0: the value of phi0, scalar
  #    Theta0s: the [theta0(X1), ..., theta0(Xn)], 1 x n
  #    H: cov mat, p x p
    
    p <-length(x)
    Xs <- as.matrix(data[, 3:(p+2)])
    sZs <- 1 - data$Z
    
    Ws <- 1/2/phi0**2 + Tau2s/2
    Ms <- data$Y/2/phi0**2 + Tau2s* Theta0s/2
    ks <- MKF(x, Xs, H=H)
    
    num <- sum(sZs*ks*Ms)
    den <- sum(sZs*ks*Ws)
    num/den
}

mOptMu0 <- function(cxs, data, Tau2s, phi0, Theta0s, H){
  #args
  #    cxs: current data point, mxp
  #    data: the datase, n x (2+p): [Y, Z, Xs]
  #    Tau2s: the [tau2(X1), ..., tau2(Xn)], 1 x n
  #    phi0: the value of phi0, scalar
  #    Theta0s: the [theta0(X1), ..., theta0(Xn)], 1 x n
  #    H: cov mat, p x p
    
    if (is.null(dim(cxs))){
        vs <- optMu0(cxs, data, Tau2s, phi0, Theta0s, H)
    }else{
        m <- dim(cxs)[1]
        p <- dim(cxs)[2]
        Xs <- as.matrix(data[, 3:(p+2)])
        sZs <- 1 - data$Z
        sZsMat <- matrix(rep(sZs, m), ncol=m)
        
        mks <- mMKF(cxs, Xs, H=H) # n x m 
        Ws <- 1/2/phi0**2 + Tau2s/2
        WsMat <- matrix(rep(Ws, m), ncol=m)
        Ms <- data$Y/2/phi0**2 + Tau2s* Theta0s/2
        MsMat <- matrix(rep(Ms, m), ncol=m)
        
        num <- colSums(sZsMat*mks*MsMat)
        den <- colSums(sZsMat*mks*WsMat)
        vs <- num/den
      
    }
    return(vs)
}

# optimize tau2(X) given [mu0(X1), ..., mu0(Xn)]
optTau <- function(x, data, Mu0s, lam, Theta0s, H, invgam2=0){
  #args
  #    x: current data point, 1xp
  #    data: the datase, n x (2+p): [Y, Z, Xs]
  #    Mu0s: the [mu0(X1), ..., mu0(Xn)], 1 x n
  #    lam: penalty parameter
  #    Theta0s: the [theta0(X1), ..., theta0(Xn)], 1 x n
  #    H: cov mat, p x p
  #    invgam2: 1/invgam2 is the variance parameter of HN prior
    p <-length(x)
    Xs <- as.matrix(data[, 3:(p+2)])
    sZs <- 1 - data$Z
    ks <- MKF(x, Xs, H=H)
    
    num <- sum(sZs * ks)
    den <- sum(sZs * ks * (Mu0s-Theta0s)**2) + invgam2
    tau2.td <- num/den
    
    tau2 <- euclidean_proj_l1ball(tau2.td, lam)
    return(tau2)
}

mOptTau <- function(cxs, data, Mu0s, lam, Theta0s, H, invgam2=0){
  
  #args
  #    cxs: current data point, mxp
  #    data: the datase, n x (2+p): [Y, Z, Xs]
  #    Mu0s: the [mu0(X1), ..., mu0(Xn)], 1 x n
  #    lam: penalty parameter
  #    Theta0s: the [theta0(X1), ..., theta0(Xn)], 1 x n
  #    H: cov mat, p x p
  #    invgam2: 1/invgam2 is the variance parameter of HN prior
    if (is.null(dim(cxs))){
        tau2 <- optTau(cxs, data, Mu0s, lam, Theta0s, H, invgam2)
    }else{
        p <-dim(cxs)[2]
        m <-dim(cxs)[1]
        
        Xs <- as.matrix(data[, 3:(p+2)])
        sZs <- 1 - data$Z
        sZsMat <- matrix(rep(sZs, m), ncol=m)
        mks <- mMKF(cxs, Xs, H=H)
        sse <- (Mu0s-Theta0s)**2
        sseMat <- matrix(rep(sse, m), ncol=m)
        
        
        num <- colSums(sZsMat * mks)
        den <- colSums(sZsMat * mks * sseMat)  + invgam2
        tau2.td <- num/den # m
        
        m <- length(tau2.td)
        tau2 <-  euclidean_proj_l1ball(tau2.td, lam*log(m))
        #tau2 <-  euclidean_proj_l1ball(tau2.td, lam*sqrt(log(m)/m))
    }
    return(tau2)
}


# can be used for ph0
optPhi0 <- function(data, Mus){
  #args
  #    data: the datase, n x (2+p): [Y, Z, Xs]
  #    Mus: the [mu0(X1), ..., mu0(Xn)], 1 x n

    a <- 0.01
    b <-  0.01
    sZs <- 1 - data$Z

    num <- sum(sZs * (data$Y - Mus)**2)/2 + b
    den <- 1 + a + sum(sZs)/2
    res <- sqrt(num/den)
    res
}

# can be used for phi1
optPhi1 <- function(data, Mus){
  #args
  #    data: the datase, n x (2+p): [Y, Z, Xs]
  #    Mus: the [mu0(X1), ..., mu0(Xn)], 1 x n

    a <- 0.01
    b <-  0.01
    Zs <- data$Z

    num <- sum(Zs * (data$Y - Mus)**2)/2 + b
    den <- 1 + a + sum(Zs)/2
    res <- sqrt(num/den)
    res
}

# function to evaluate estimated mu0(X)
mu0.efn <- function(xs, res){
    rv <- mOptMu0(xs, res$data, res$tau2s, res$phi0, res$theta0s, res$H)
    rv
}

mu1.efn <- function(xs, res){
    rv <- mOptMu1(xs, res$data, res$invsigma2, res$phi1, res$theta1s, res$H)
    rv
}

# function to evaluate estimated tau2(X)
tau2.efn <- function(xs, res){
    rv <- mOptTau(xs, res$data, res$mu0s, res$lam, res$theta0s, res$H, res$invgam2)
    rv
}

# calculate the posterior variance of the mu0(X)
post.var.mu0.fn <- function(cxs, res){
  #args
  # cxs, data points, m x p
  
  if (is.null(dim(cxs))){
      xs <- rbind(cxs, cxs)
  }else{
      xs <- cxs
  }
      
  m <- dim(xs)[1]
  p <- dim(xs)[2]
  Xs <- as.matrix(res$data[, 3:(p+2)])
  tau2Mat <- matrix(rep(res$tau2s, m), ncol=m)
  sZs <- 1 - res$data$Z
  sZsMat <- matrix(rep(sZs, m), ncol=m)
  
  mks <- mMKF(xs, Xs, res$H)
  ws.sum <- colSums(mks)
  ws.sum.mat <- matrix(rep(ws.sum, dim(Xs)[1]), ncol=m, byrow=TRUE)
  norm.mks <- mks/ws.sum.mat
  rv <- colSums(sZsMat*norm.mks*(1/res$phi0/res$phi0+tau2Mat))
  
  if (is.null(dim(cxs))){
      rv <- rv[1]
  }
  
  return(1/rv)
}

# calculate the posterior mean of the mu0(X)
post.mean.mu0.fn <- function(cxs, res){
  #args
  # cxs, data points, m x p
  
  if (is.null(dim(cxs))){
      xs <- rbind(cxs, cxs)
  }else{
      xs <- cxs
  }
      
  m <- dim(xs)[1]
  p <- dim(xs)[2]
  Xs <- as.matrix(res$data[, 3:(p+2)])
  tau2Mat <- matrix(rep(res$tau2s, m), ncol=m)
  sZs <- 1 - res$data$Z
  sZsMat <- matrix(rep(sZs, m), ncol=m)
  Y <- res$data$Y
  YMat <- matrix(rep(Y, m), ncol=m)
  theta0Mat <- matrix(rep(res$theta0s, m), ncol=m)
  
  mks <- mMKF(xs, Xs, res$H)
  ws.sum <- colSums(mks)
  ws.sum.mat <- matrix(rep(ws.sum, dim(Xs)[1]), ncol=m, byrow=TRUE)
  norm.mks <- mks/ws.sum.mat
  den <- colSums(sZsMat*norm.mks*(1/res$phi0/res$phi0+tau2Mat))
  
  num <- colSums(sZsMat*norm.mks*(YMat/res$phi0/res$phi0+tau2Mat*theta0Mat))
  rv <- num/den
  
  if (is.null(dim(cxs))){
      rv <- rv[1]
  }
  
  return(rv)
}

# function for estimate mu0s, tau2s, and phi0 under info and reference model
mu0.info.est.fn <- function(Theta0s, data, H, lam, phi0=NA, invgam2=0, is.ref=FALSE, maxit=100){
  #args:
  #  Theta0s: Estimated ys based on historical data
  #  data: the dataset, n x (2+p): [Y, Z, Xs]
  #  H: cov mat, p x p
  #  lam: The penalty parameters
  #  is.ref: if true, reture the results of reference model distribution
  #  maxit: Maximal number of times for iteration
    
    phi0.tk <- c()
    Mu0s.tk <- list()
    Tau2s.tk <- list()
    p <- dim(data)[2] - 2
    Xs <- data[, 3:(p+2)]
    n <- dim(Xs)[1]
    Tau2s <- rep(0, n)


    if (!is.ref| is.na(phi0)){
        phi0 <- 1
        for (i in 1:maxit){
            Mu0s <- mOptMu0(as.matrix(Xs), data, Tau2s, phi0, Theta0s, H)
            Mu0s.tk[[i]] <- Mu0s

            phi0 <- optPhi0(data, Mu0s)
            phi0.tk[i] <- phi0
            
            Tau2s  <- mOptTau(as.matrix(Xs), data, Mu0s, lam, Theta0s, H, invgam2)
            Tau2s.tk[[i]] <- Tau2s
            
            if (i > 1){
                err.tau2 <- mean((Tau2s.tk[[i]] - Tau2s.tk[[i-1]])**2)
                err.mu <- mean((Mu0s.tk[[i]] - Mu0s.tk[[i-1]])**2)
                err.phi0 <- mean((phi0.tk[[i]] - phi0.tk[[i-1]])**2)
                
                err.all <- max(c(err.tau2, err.mu, err.phi0))
                if (err.all < 1e-5){
                    break
                }
            }
        }
   }

    if (!is.ref){
        res <- list(tau2s=Tau2s, mu0s=Mu0s, phi0=phi0, mu0tk=Mu0s.tk, tau2stk=Tau2s.tk,
                theta0s=Theta0s, H=H, data=data, lam=lam, phi0.tk=phi0.tk, invgam2=invgam2)
    }else{
        Tau20s <- rep(0, n)
        Mu0s.no <- mOptMu0(as.matrix(Xs), data, Tau20s, phi0, Theta0s, H)
        res <- list(tau2s=Tau20s, mu0s=Mu0s.no, phi0=phi0, theta0s=Theta0s, H=H, data=data, lam=lam, invgam2=invgam2)
    }
    
    res
}

# estimate mu0, phi0 when non-borrow
mu0.no.est.fn <- function(data, H){
  #args:
  #  data: the dataset, n x (2+p): [Y, Z, Xs]
  #  H: cov mat, p x p
    
    p <- dim(data)[2] - 2
    Xs <- data[, 3:(p+2)]
    n <- dim(Xs)[1]
    Tau2s <- rep(0, n)
    phi0 <- 1
    Mu0s <- mOptMu0(as.matrix(Xs), data, Tau2s, phi0, Tau2s, H)
    phi0 <- optPhi0(data, Mu0s)

    res <- list(mu0s=Mu0s, phi0=phi0, H=H, data=data, tau2s=Tau2s, theta0s=rep(0, n))
    res
}

# function for estimate mu1s, phi when borrowing
mu1.info.est.fn <- function(Theta1s, data, H, invsigma2=0, maxit=100){
  #args:
  #  Theta1s: Estimated ys based on historical data
  #  data: the dataset, n x (2+p): [Y, Z, Xs]
  #  H: cov mat, p x p
  #  invsigma2: 1/sigma2
  #  maxit: Maximal number of times for iteration
    
    phi1.tk <- c()
    Mu1s.tk <- list()
    p <- dim(data)[2] - 2
    Xs <- data[, 3:(p+2)]
    n <- dim(Xs)[1]

    phi1 <- 1
    for (i in 1:maxit){
        Mu1s <- mOptMu1(as.matrix(Xs), data, invsigma2, phi1, Theta1s, H)
        Mu1s.tk[[i]] <- Mu1s

        phi1 <- optPhi1(data, Mu1s)
        phi1.tk[i] <- phi1
        
        if (i > 1){
            err.mu <- sum((Mu1s.tk[[i]] - Mu1s.tk[[i-1]])**2)
            err.phi1 <- sum((phi1.tk[[i]] - phi1.tk[[i-1]])**2)
            
            err.all <- max(c(err.mu, err.phi1))
            if (err.all < 1e-5){
                break
            }
        }
    }

    res <- list(mu1s=Mu1s, phi1=phi1, invsigma2=invsigma2, 
                   phi1.tk=phi1.tk, Mu1s.tk=Mu1s.tk,
                   theta1s=Theta1s, H=H, data=data)
    res
}

# function for estimate mu1s, phi when non-borrowing
mu1.no.est.fn <- function(data, H){
  #args:
  #  data: the dataset, n x (2+p): [Y, Z, Xs]
  #  H: cov mat, p x p
    
    p <- dim(data)[2] - 2
    Xs <- data[, 3:(p+2)]
    n <- dim(Xs)[1]
    phi1 <- 1

    Mu1s <- mOptMu1(as.matrix(Xs), data, 0, phi1, rep(0, n), H)
    phi1 <- optPhi1(data, Mu1s)


    res <- list(mu1s=Mu1s, phi1=phi1, invsigma2=0, H=H, data=data, theta1s=rep(0, n))
}

# calculate the posterior variance of the mu1(X)
post.var.mu1.fn <- function(cxs, res){
  #args
  # cxs, data points, m x p
  
  if (is.null(dim(cxs))){
      xs <- rbind(cxs, cxs)
  }else{
      xs <- cxs
  }
      
  m <- dim(xs)[1]
  p <- dim(xs)[2]
  Xs <- as.matrix(res$data[, 3:(p+2)])
  Zs <- res$data$Z
  ZsMat <- matrix(rep(Zs, m), ncol=m)
  
  mks <- mMKF(xs, Xs, res$H)
  ws.sum <- colSums(mks)
  ws.sum.mat <- matrix(rep(ws.sum, dim(Xs)[1]), ncol=m, byrow=TRUE)
  norm.mks <- mks/ws.sum.mat
  rv <- colSums(ZsMat*norm.mks*(1/res$phi1/res$phi1+res$invsigma2))
  
  if (is.null(dim(cxs))){
      rv <- rv[1]
  }
  
  return(1/rv)
}

# calculate the posterior mean of the mu1(X)
post.mean.mu1.fn <- function(cxs, res){
  #args
  # cxs, data points, m x p
  
  if (is.null(dim(cxs))){
      xs <- rbind(cxs, cxs)
  }else{
      xs <- cxs
  }
      
  m <- dim(xs)[1]
  p <- dim(xs)[2]
  Xs <- as.matrix(res$data[, 3:(p+2)])
  Zs <-  res$data$Z
  ZsMat <- matrix(rep(Zs, m), ncol=m)
  Ys <- res$data$Y
  YsMat <- matrix(rep(Ys, m), ncol=m)
  theta1Mat <- matrix(rep(res$theta1s, m), ncol=m)
  
  mks <- mMKF(xs, Xs, res$H)
  ws.sum <- colSums(mks)
  ws.sum.mat <- matrix(rep(ws.sum, dim(Xs)[1]), ncol=m, byrow=TRUE)
  norm.mks <- mks/ws.sum.mat
  den <- colSums(ZsMat*norm.mks*(1/res$phi1/res$phi1+res$invsigma2))
  num <- colSums(ZsMat*norm.mks*(YsMat/res$phi1/res$phi1+theta1Mat*res$invsigma2))
  rv <- num/den
  
  if (is.null(dim(cxs))){
      rv <- rv[1]
  }
  
  return(rv)
}

# sampling from posterior distribution of mu0
r.postMu0 <- function(cxs, res, M){
    # M: num of repitions
  ms <- post.mean.mu0.fn(cxs, res)
  vs <- post.var.mu0.fn(cxs, res)
  nSps <- length(ms)
  spss <- lapply(1:M, function(i)rnorm(nSps, ms, sqrt(vs)))
  trts <- colMeans(do.call(cbind, spss))
  list(spss=spss, trts=trts)
}

# sampling from posterior distribution of mu1
r.postMu1 <- function(cxs, res, M){
    # M: num of repitions
  ms <- post.mean.mu1.fn(cxs, res)
  vs <- post.var.mu1.fn(cxs, res)
  nSps <- length(ms)
  spss <- lapply(1:M, function(i)rnorm(nSps, ms, sqrt(vs)))
  trts <- colMeans(do.call(cbind, spss))
  list(spss=spss, trts=trts)
}

