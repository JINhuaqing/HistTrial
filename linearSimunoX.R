# this script is to run simulations
# to compare considering X and not considering X
rm(list=ls())
library(magrittr)
library(dplyr)
library(parallel)
#setwd("/root/Documents/HQ/HistTrial")
setwd("/Users/hujin/Library/CloudStorage/OneDrive-UCSF/Documents/ProjectCode/HistTrial")
source("utils.R")
source("utils_noX.R")
source("simplex.R")
seeds <- 1:10000


CI.fn <- function(errs){
  m.v <- mean(errs)
  sd.v <- sd(errs)
  se <- sd.v/sqrt(length(errs))
  low <- m.v - 1.96*se
  up <- m.v + 1.96*se
  rvs <- c(low, m.v, up)
  names(rvs) <- c("lower", "mean", "upper")
  rvs
}

gen.bin.Xs <- function(n, p=0.5){
  res <- runif(n) < p 
  res <- as.numeric(res)
  res
}

# this function is to generate the historical parameters 
# for corrent and wrong models
gen.hist.paras <- function(n.h){
    Xs.h <- gen.bin.Xs(n.h)
    Y.h <- alps[1] + Xs.h*alps[2] + rnorm(n.h)*phi0.h
    fit.h <- lm(Y.h~Xs.h)
    true.model.ests <- fit.h$coefficients
    wrong.model.ests <- mean(Y.h)
    res <- list(true=true.model.ests,
                wrong=wrong.model.ests
                )
    res
}


fun.real <- function(i){
    ts <- c()
    ts <- c(ts, Sys.time())
    set.seed(seeds[i])
    alpss.both <- gen.hist.paras(n.h) # the historical model
    print(i)
    
    # first stage: simply randomization for first n0 sps
    Xs <- gen.bin.Xs(n0)
    idx0 <- sample.int(n0, size=floor(n0/2))
    Zs <- rep(1, n0)
    Zs[idx0] <- 0
    nerrs <- rnorm(n0, sd=phi0)
    
    # considering X
    Zs.true <- Zs
    Ys.true <- b*Zs + betas[1] + betas[2]*Xs + nerrs
    data.true <- cbind(Ys.true, Zs.true, Xs)
    data.true <- as.data.frame(data.true)
    colnames(data.true)[1:2] <- c("Y", "Z")
    
    
    # no considiering X
    Zs.wrong <- Zs
    Ys.wrong <- b*Zs + betas[1] + betas[2]*Xs + nerrs
    data.wrong <- cbind(Ys.wrong, Zs.wrong, Xs)
    data.wrong <- as.data.frame(data.wrong)
    colnames(data.wrong)[1:2] <- c("Y", "Z")
    
    # some vecs for recording
    Rs.true <- rep(1, N) # record R 
    Rs.wrong <- rep(1, N) # record R 
    lam.trus <- rep(-1, N+1) # record lambda
    ts <- c(ts, Sys.time())
    
    lastRes.true <- NA
    lastRes.wrong <- NA
    for (j in (n0+1):N){
        #print(j)
        # new subject
        cx <- gen.bin.Xs(1)
        
        
        Theta0s.true <- alpss.both$true[1] + alpss.both$true[2]*Xs
        Theta0s.wrong <- alpss.both$wrong
        
        # with X
        res.true <- mu0.info.est.fn(Theta0s.true, data.true, H, lam, invgam2=invgam2, lam.tru=lam.tru, lastRes=lastRes.true)
        lastRes.true <- res.true
        res.true.ref <- mu0.info.est.fn(Theta0s.true, data.true, H, lam, phi0=res.true$phi0, invgam2=invgam2, is.ref=TRUE)
        
        var.true.info <- post.var.mu0.fn(cx, res.true)
        var.true.ref <- post.var.mu0.fn(cx, res.true.ref)
        R.true <- var.true.ref/var.true.info
        Rs.true[j] <- R.true
        
        ass.true.res <- RBC.design(cx, data.true[, 3:(p+2), drop=F], data.true$Z, hs, R.true)
        
        
        # without X
        res.wrong <- mu0.info.est.noX.fn(Theta0s.wrong[1], data.wrong, 
                                         lam, invgam2=invgam2, lam.tru=0, lastRes=lastRes.wrong)
        lastRes.wrong <- res.wrong
        res.wrong.ref <- mu0.info.est.noX.fn(Theta0s.wrong[1], data.wrong, lam, phi0=res.wrong$phi0, 
                                             invgam2=invgam2, is.ref=TRUE)
        
        var.wrong.info <- post.var.mu0.noX.fn(res.wrong)
        var.wrong.ref <- post.var.mu0.noX.fn(res.wrong.ref)
        R.wrong <- var.wrong.ref/var.wrong.info
        Rs.wrong[j] <- R.wrong
        ass.wrong.res <- RBC.design(cx, data.wrong[, 3:(p+2), drop=F], data.wrong$Z, hs, R=R.wrong)
        
        Xs <- c(Xs, cx)
        Zs.true <- c(Zs.true, ass.true.res$grp-1)
        Zs.wrong <- c(Zs.wrong, ass.wrong.res$grp-1)
        
        nerr <- rnorm(1, sd=phi0)
        curN <- length(Xs)
        
        # add new to data.true
        Y.true <- b*Zs.true[curN] + betas[1] + betas[2]*Xs[curN] + nerr
        Ys.true <- c(Ys.true, Y.true)
        data.true <- cbind(Ys.true, Zs.true, Xs)
        rownames(data.true) <- NULL
        data.true <- as.data.frame(data.true)
        colnames(data.true)[1:2] <- c("Y", "Z")
        
        # add new to data.wrong
        Y.wrong <- b*Zs.wrong[curN] + betas[1] + betas[2]*Xs[curN] + nerr
        Ys.wrong <- c(Ys.wrong, Y.wrong)
        data.wrong <- cbind(Ys.wrong, Zs.wrong, Xs)
        rownames(data.wrong) <- NULL
        data.wrong <- as.data.frame(data.wrong)
        colnames(data.wrong)[1:2] <- c("Y", "Z")
        
    }
    ts <- c(ts, Sys.time())
    
    Theta0s.true <- alpss.both$true[1] + alpss.both$true[2]*Xs
    Theta0s.wrong <- alpss.both$wrong
    
    # with X
    res.true <- mu0.info.est.fn(Theta0s.true, data.true, H, lam, invgam2=invgam2, lam.tru=0, lastRes=lastRes.true)
    res.true.mu1 <- mu1.no.est.fn(data.true, H)
    testPt <- as.matrix(Xs)
    trt.eff.true <- mean(mu1.efn(testPt, res.true.mu1)) -  mean(mu0.efn(testPt, res.true))
    msps0.true <- r.postMu0(testPt, res.true, M)$trts
    msps1.true <- r.postMu1(testPt, res.true.mu1, M)$trts
    sps.trts.true <- msps1.true - msps0.true
        
    # without X
    res.wrong <- mu0.info.est.noX.fn(Theta0s.wrong[1], data.wrong, 
                                         lam, invgam2=invgam2, lam.tru=0, lastRes=lastRes.wrong)
    res.wrong.mu1 <- mu1.no.est.noX.fn(data.wrong)
    trt.eff.wrong <- mean(mu1.noX.efn(res.wrong.mu1)) -  mean(mu0.noX.efn(res.wrong))
    msps0.wrong <- r.postMu0.noX(res.wrong, M)$trts
    msps1.wrong <- r.postMu1.noX(res.wrong.mu1, M)$trts
    sps.trts.wrong <- msps1.wrong - msps0.wrong

    sps.trts <- cbind(sps.trts.true, sps.trts.wrong)

    
    ts <- c(ts, Sys.time())
    #print(diff(ts))
    rv <- list(mtrt=c(trt.eff.true, trt.eff.wrong),
               sps.trts=sps.trts, 
               data.true=data.true, data.wrong=data.wrong, 
               res.true=res.true, 
               res.wrong=res.wrong, 
               Rs.true=Rs.true, 
               Rs.wrong=Rs.wrong, 
               alpss.both=alpss.both, 
               betas=betas)
    rv
}


# True model parameters
#under H0
b<-1 # the treatment effect

betas <- c(1, 1) # parameters for current model
alps <- c(6, -9) # parameters for historical model
#alps <- c(1, 10) # parameters for historical model
phi0 = phi1 = 0.5 # sd of the noise for the current trial in two grps
phi0.h = 0.5 # sd of the noise for the historical model


# Simu parameters
n.h <- 100
invgam2 <- 0.33
N <- 200 # total sample size
# parameters
lam.tru <- 0.0 # two subgrps can not use lam.q, so I give lam.tru directly.
lam <- 300
hs <- c(1.1) # kernel bw for allocation
H <- diag(1) * 0.1 # kernel bw for reg
p <- 1
# initial dataset
n0 <- 20

# for hypothesis test
# to calculate the prob
M <- 1000


nSimu <- 2000
post.res <- mclapply(1:nSimu, fun.real, mc.cores=35)
H <- diag(c(0.10, 0.10, 9.99, 9.99))
sv.name <- paste0("./results/linearSimu_noX_overallsame_H0.RData")
paras <- list(invgam2=invgam2, 
              b=b, phi0=phi0, phi1=phi1, N=N, lam.tru=lam.tru, 
              hs=hs,  H=H, M=M, lam=lam)
save(post.res, paras, file=sv.name)

