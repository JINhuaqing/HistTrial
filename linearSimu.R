#rm(list=ls())
library(magrittr)
library(dplyr)
#setwd("/home/huaqingj/MyResearch/HistTrial")
setwd("/home/r13user3/MyResearch/HistTrial")
#setwd("C:/Users/JINHU/Documents/ProjectCode/HistTrial")
source("utils.R")
source("simplex.R")
library(parallel)
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

fun.test <- function(i){
    #set.seed(seeds[i])
    alpss <- betass
    for (jj in 1:4){
        alpss[[jj]] <- betass[[jj]] + rnorm(p, sd=xis[jj])
    }
    print(i)
    Xs <- gen.Data.Xs(n0, x.tps)
    idx0 <- sample.int(n0, size=floor(n0/2))
    Zs <- rep(1, n0)
    Zs[idx0] <- 0
    
    
    betMat <- sub.Paras.fn(Xs, betass)
    nerrs <- rnorm(n0, sd=phi0)
    Ys <- curMean.fn(Xs, Zs, betMat, b) + nerrs
    data <- cbind(Ys, Zs, Xs)
    data <- as.data.frame(data)
    colnames(data)[1:2] <- c("Y", "Z")
    
    
    # no borrowing 
    Zs.no <- Zs
    Ys.no <- curMean.fn(Xs, Zs.no, betMat, b) + nerrs
    data.no <- cbind(Ys.no, Zs.no, Xs)
    data.no <- as.data.frame(data.no)
    colnames(data.no)[1:2] <- c("Y", "Z")
    Rs <- rep(1, N) # record R 
    
    for (j in (n0+1):N){
        cx <- unlist(gen.Data.Xs(1, x.tps))
        #print(j)
        
        # H <- diag(c(bw.nrd(data$X1), bw.nrd(data$X2), bw.nrd(data$X3), bw.nrd(data$X4)))
        alpMat <- sub.Paras.fn(Xs, alpss)
        Theta0s <- curMean.fn(Xs, Zs, alpMat, b=0)
        res <- mu0.info.est.fn(Theta0s, data, H, lam, invgam2=invgam2)
        res.ref <- mu0.info.est.fn(Theta0s, data, H, lam, phi0=res$phi0, invgam2=invgam2, is.ref=TRUE)
        
        var.info <- post.var.mu0.fn(cx, res)
        var.ref <- post.var.mu0.fn(cx, res.ref)
        R <- var.ref/var.info
        Rs[j] <- R
        ass.res <- RBC.design(cx, data[, 3:(p+2)], data$Z, hs, R)
        ass.res.no <- RBC.design(cx, data[, 3:(p+2)], data.no$Z, hs, R=1)
        # ass.res <- RBC.design(cx, data[, 3:(p+2)], data$Z, hs, R)
        
        Xs <- rbind(Xs, cx)
        Zs <- c(Zs, ass.res$grp-1)
        Zs.no <- c(Zs.no, ass.res.no$grp-1)
        
        curN <- dim(Xs)[1]
        curBetMat <- sub.Paras.fn(Xs[curN, ], betass)
        nerr <- rnorm(1, sd=phi0)
        Y <- curMean.fn(Xs[curN, ], Zs[curN], curBetMat, b) + nerr
        Ys <- c(Ys, Y)
        data <- cbind(Ys, Zs, Xs)
        data <- as.data.frame(data)
        colnames(data)[1:2] <- c("Y", "Z")
        
        Y.no <- curMean.fn(Xs[curN, ], Zs.no[curN], curBetMat, b) + nerr
        Ys.no <- c(Ys.no, Y.no)
        data.no <- cbind(Ys.no, Zs.no, Xs)
        data.no <- as.data.frame(data.no)
        colnames(data.no)[1:2] <- c("Y", "Z")
        
        
    }
    
    alpMat <- sub.Paras.fn(Xs, alpss)
    Theta0s <- curMean.fn(Xs, Zs, alpMat, b=0)
    res <- mu0.info.est.fn(Theta0s, data, H, lam, invgam2=invgam2)
    res.no <- mu0.no.est.fn(data.no, H)
    
    res.mu1 <- mu1.no.est.fn(data, H)
    res.no.mu1 <- mu1.no.est.fn(data.no, H)
    #zeroPt <- matrix(rep(0, p), nrow=1)
    zeroPt <- as.matrix(Xs)

    #trt.eff <- mean(mu1.efn(zeroPt, res.mu1)) -  mean(mu0.efn(zeroPt, res))
    trt.eff <- mean(mu1.efn(as.matrix(Xs), res.mu1)) -  mean(mu0.efn(as.matrix(Xs), res))
    #trt.eff.no <- mean(mu1.efn(zeroPt, res.no.mu1)) -  mean(mu0.efn(zeroPt, res.no))
    trt.eff.no <- mean(mu1.efn(as.matrix(Xs), res.no.mu1)) -  mean(mu0.efn(as.matrix(Xs), res.no))
    post.prob.trt <- function(i){
        sps0 <- r.postMu0(zeroPt, res)
        sps0.no <- r.postMu0(zeroPt, res.no)
        sps1 <- r.postMu1(zeroPt, res.mu1)
        sps1.no <- r.postMu1(zeroPt, res.no.mu1)
        
        sps.trt <- mean(sps1) - mean(sps0)
        sps.trt.no <- mean(sps1.no) - mean(sps0.no)
        c(sps.trt, sps.trt.no)
    }
    sps.trts <- lapply(1:M, post.prob.trt)
    
    # print(res$tau2s)
    rv <- list(mtrt=c(trt.eff, trt.eff.no), sps.trts=sps.trts, data=data, data.no=data.no, res=res, Rs=Rs)
    rv
}

# parameters for current data 
betass <- list(
              para1=c(2, 1, -1, 3, -2), 
              para2=c(3, 1, 2, 4, 2), 
              para3=c(0, -1, -1, 3, 2), 
              para4=c(2, 0, -1, 4, -2)
                  )


# parameters for historical data
#alpss <-  list(para1=c(2, 1, -1, 3, -2), 
#             para2=c(3, 1, 2, 4, 2), 
#             para3=c(0, -1, -1, 3, 2), 
#             para4=c(2, 0, -1, 4, -2) )
#
# sd of random error on alpha
xis <- c(0, 0, 0, 0)
invgam2 <- 1

b <- 4
phi0 = phi1 = 2
N <- 150 # total sample size
# parameters
lam <- 20
hs <- rep(2.1, 4)
x.tps <- c(2, 2, "c", "c")
x.tps <- c(2, 2, 2, 2)
p <- length(x.tps)

#H <- diag(rep(1, p)/2)
H <- diag(rep(0.1, p))

# initial dataset
n0 <- 20

# for hypothesis test
# to calculate the prob
M <- 1000


nSimu <- 1000
for (N in c(30, 60, 90, 120, 150)){
        paras <- list(invgam2=invgam2, b=b, phi0=phi0, phi1=phi1, N=N, lam=lam, hs=hs, x.tps=x.tps, H=H, M=M, xis=xis)
        post.res <- mclapply(1:nSimu, fun.test, mc.cores=40)
        sv.name <- paste0("./results/Simplex-Linear-AllBinary-diff-0", "-b-", b, "-N-", N, "-lam-", lam, "-phi0-", phi0, "-invgam2-", invgam2, "-nSimu-", nSimu, ".RData")
        #sv.name <- paste0("./results/Simplex-Linear-diff-Sub10", "-b-", b, "-N-", N, "-lam-", lam, "-phi0-", phi0, "-invgam2-", invgam2, "-nSimu-", nSimu, ".RData")
        save(post.res, paras, file=sv.name)
}
