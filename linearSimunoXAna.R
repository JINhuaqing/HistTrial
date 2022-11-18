rm(list=ls())
setwd("/Users/hujin/Library/CloudStorage/OneDrive-UCSF/Documents/ProjectCode/HistTrial")
source("./util_ana.R")




# the F statistis to measure the imbalance between two arms
F.stat <- function(data){
  Xs <- data$Xs
  Z <- data$Z
  n1 <- sum(Z)
  n0 <- length(Z) - n1
  
  tm0 <- n0*(mean(Xs[Z==0]) -mean(Xs))**2
  tm1 <- n1*(mean(Xs[Z==1]) -mean(Xs))**2
  den <- sum((Xs - mean(Xs))**2)
  v <- (tm0+tm1)/den
  return(v)
}

# the allocation rate for each subgrp
subgrp.Z1.rate.fn <- function(data){
    Xs <- data$Xs
    Z <- data$Z
    rates <- c(mean(Z[Xs==0]), mean(Z[Xs==1]))
    rates
}

get.results <- function(post.res){
    # trt eff est
    mtrts <- do.call(rbind, lapply(post.res, function(ix)ix$mtrt))
    
    # z1 rates
    Z1rate.true <-  sapply(post.res, function(x)mean(x$data.true$Z))
    Z1rate.wrong <-  sapply(post.res, function(x)mean(x$data.wrong$Z))
    Z1rates <- cbind(Z1rate.true, Z1rate.wrong)
    
    # sub allo rate
    Z1rate.sub.true <- do.call(rbind, lapply(post.res, function(x)subgrp.Z1.rate.fn(x$data.true)))
    Z1rate.sub.wrong <- do.call(rbind, lapply(post.res, function(x)subgrp.Z1.rate.fn(x$data.wrong)))
    
    # prob > 0
    prb.true <- sapply(post.res, function(x) post.fn(x$sps.trts[, 1], 0))
    prb.wrong <- sapply(post.res, function(x) post.fn(x$sps.trts[, 2], 0))
    prbs <- cbind(prb.true, prb.wrong)
    
    # imblance metric between 2 arms
    f.true <- sapply(post.res, function(x) F.stat(x$data.true))
    f.wrong <- sapply(post.res, function(x) F.stat(x$data.wrong))
    fs <- cbind(f.true, f.wrong)
    
    rv <- list(
      mtrts = mtrts, 
      Z1rates = Z1rates,
      Z1rate.sub.true = Z1rate.sub.true,
      Z1rate.sub.wrong = Z1rate.sub.wrong,
      prbs = prbs,
      fs = fs
    )
    rv
}

H0.fil <- "./results/linearSimu_noX_1grpsame_H0.RData" 
#H0.fil <- "./results/linearSimu_noX_overallsame_H0.RData" 
load(H0.fil)
res.H0 <- get.results(post.res)
c.true <- eps.alig(res.H0$prbs[, 1], 0.1);c.true
c.wrong <- eps.alig(res.H0$prbs[, 2], 0.1);c.wrong

H1.fil <- "./results/linearSimu_noX_1grpsame_b03_H1.RData" 
#H1.fil <- "./results/linearSimu_noX_overallsame_b03_H1.RData" 
load(H1.fil)
res.H1 <- get.results(post.res)

# powers
power.true <- mean(res.H1$prbs[, 1] >  c.true$eps0);power.true
power.wrong <- mean(res.H1$prbs[, 2] >  c.wrong$eps0);power.wrong

# abs err
errs <- abs(res.H1$mtrts - 0.3)
err.true <- CI.fn(errs[, 1]);err.true
err.wrong <- CI.fn(errs[, 2]);err.wrong

# Z1rate
Z1.true <- CI.fn(res.H1$Z1rates[, 1]);Z1.true
Z1.wrong <- CI.fn(res.H1$Z1rates[, 2]);Z1.wrong
colMeans(res.H1$Z1rate.sub.true)
colMeans(res.H1$Z1rate.sub.wrong)


# allo F stat
f.true <- CI.fn(res.H1$fs[, 1]);f.true
f.wrong <- CI.fn(res.H1$fs[, 2]);f.wrong


