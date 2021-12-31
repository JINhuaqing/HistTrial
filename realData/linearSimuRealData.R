#rm(list=ls())
library(magrittr)
library(dplyr)
#setwd("/home/huaqingj/MyResearch/HistTrial")
#setwd("/home/r13user3/MyResearch/HistTrial")
setwd("C:/Users/JINHU/Documents/ProjectCode/HistTrial")
source("utils.R")
source("simplex.R")
library(parallel)
seeds <- 1:10000

sub.Paras.fn <- function(Xs, betass){
  #args:
  #   Xs: covariates, n x p 
  if (is.null(dim(Xs))){
    
      idxs <- Xs[1] +  Xs[2]*2 + 1
  }else{
      idxs <- Xs[, 1] +  Xs[, 2]*2 + 1
    
  }
  betMat <- do.call(rbind, betass[idxs])
  return(betMat)
  
}

vec2code <- function(vec,fct=10){
    vecF <- vec * fct
    pre <- c()
    if (fct == 10){
        for (ivecF in vecF) {
            if (ivecF < fct){
                pre <- c(pre, "0")
            }else{
                pre <- c(pre, "")
            }
        }
    }else if (fct==100){
        for (ivecF in vecF) {
            if (ivecF >= 100){
                pre <- c(pre, "")
            }else if(ivecF >=10){
                pre <- c(pre, "0")
            }else{
                pre <- c(pre, "00")
            }
        }
    }
    idxCode <- paste0(paste0(pre, vecF), collapse="")
    idxCode
}

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

fun.real <- function(i){
    #set.seed(seeds[i])
    alpss <- betass
    for (jj in 1:4){
       # error is Uniform(1-xi*|beta|, 1+xi*|beta|) * beta
        alpss[[jj]] <- betass[[jj]] *  runif(p+1, min=1-xis[[jj]]*abs(betass[[jj]]), 
                                             max=1+xis[[jj]]*abs(betass[[jj]]))
    }
    print(i)
    Xs <- gen.Real.Xs(n0, fHats)
    
    idx0 <- sample.int(n0, size=floor(n0/2))
    Zs <- rep(1, n0)
    Zs[idx0] <- 0
    
    
    betMat <- sub.Paras.fn(Xs, betass)
    Ys.m <- curMean.fn(Xs, Zs, betMat, b)
    nerrs <- rnorm(n0, sd=phi0)
    Ys <- Ys.m + nerrs
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
    lam.trus <- rep(-1, N+1) # record lambda
    
    for (j in (n0+1):N){
        #print(j)
        cx <- unlist(gen.Real.Xs(1, fHats))
        
        
        # H <- diag(c(bw.nrd(data$X1), bw.nrd(data$X2), bw.nrd(data$X3), bw.nrd(data$X4)))
        alpMat <- sub.Paras.fn(Xs, alpss)
        Theta0s <- curMean.fn(Xs, Zs, alpMat, b=0)
        lam.tru <- lam.sel.fn(data, H=H, invgam2=invgam2, lam=lam, lam.q=lam.q)
        lam.trus[j] <- lam.tru
        res <- mu0.info.est.fn(Theta0s, data, H, lam, invgam2=invgam2, lam.tru=lam.tru)
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
        rownames(data) <- NULL
        data <- as.data.frame(data)
        colnames(data)[1:2] <- c("Y", "Z")
        
        Y.no <- curMean.fn(Xs[curN, ], Zs.no[curN], curBetMat, b) + nerr
        Ys.no <- c(Ys.no, Y.no)
        data.no <- cbind(Ys.no, Zs.no, Xs)
        rownames(data.no) <- NULL
        data.no <- as.data.frame(data.no)
        colnames(data.no)[1:2] <- c("Y", "Z")
        
        
    }
    
    alpMat <- sub.Paras.fn(Xs, alpss)
    Theta0s <- curMean.fn(Xs, Zs, alpMat, b=0)
    lam.tru <- lam.sel.fn(data, H=H, invgam2=invgam2, lam=lam, lam.q=lam.q)
    lam.trus[N+1] <- lam.tru
    res <- mu0.info.est.fn(Theta0s, data, H, lam, invgam2=invgam2, lam.tru=lam.tru)
    res.no <- mu0.no.est.fn(data.no, H)
    
    res.mu1 <- mu1.no.est.fn(data, H)
    res.no.mu1 <- mu1.no.est.fn(data.no, H)
    #testPt <- matrix(rep(0, p), nrow=1)
    testPt <- as.matrix(Xs)

    trt.eff <- mean(mu1.efn(testPt, res.mu1)) -  mean(mu0.efn(testPt, res))
    trt.eff.no <- mean(mu1.efn(testPt, res.no.mu1)) -  mean(mu0.efn(testPt, res.no))

    msps0 <- r.postMu0(testPt, res, M)$trts
    msps0.no <- r.postMu0(testPt, res.no, M)$trts
    msps1 <- r.postMu1(testPt, res.mu1, M)$trts
    msps1.no <- r.postMu1(testPt, res.no.mu1, M)$trts
    sps.trts <- msps1 - msps0
    sps.trts.no <- msps1.no - msps0.no

    sps.trts <- cbind(sps.trts, sps.trts.no)

    
    # print(res$tau2s)
    rv <- list(mtrt=c(trt.eff, trt.eff.no), sps.trts=sps.trts, data=data, data.no=data.no, res=res, Rs=Rs, lam.trus=lam.trus)
    rv
}

# sd of random error on alpha
xi.vs <- c(seq(0, 1, 0.2), 1.3, 1.6, 1.9, 2, 3, 5, 7, 9)
xiss <- list()
for (ii in 1:length(xi.vs)){
    xiss[[ii]] <- c(0, 0, 0, 0) + xi.vs[ii]
}
xis <- c(0, 1, 1, 0)

invgam2 <- 1

phi0 = phi1 = sd(fit.all$residuals)
N <- 150 # total sample size
# parameters
lam.q <- 0.50 
lam <- 100
hs <- c(1.1, 1.1, 2.1, 2.1)
x.tps <- c(2, 2, "c", "c")
idx.tps <- paste0(x.tps, collapse="")
p <- length(x.tps)

H <- diag(c(0.05, 0.05, 0.5, 0.5))
#H <- diag(c(0.5, 0.5, 0.5, 0.5))

# initial dataset
n0 <- 20

# for hypothesis test
# to calculate the prob
M <- 1000

nSimu <- 100
paras <- list(invgam2=invgam2, b=b, phi0=phi0, phi1=phi1, N=N, lam.q=lam.q, hs=hs, x.tps=x.tps, H=H, M=M, xis=xis, lam=lam)
post.res <- lapply(1:nSimu, fun.real)
#post.res <- mclapply(1:nSimu, fun.real, mc.cores=4)
sv.name <- paste0("./results/RealData", "-b-", format(b, scientific=T, digits=1), "-N-", N, "-lam-", lam, "-lamq-", 100*lam.q, "-phi0-", format(phi0, scientific=T, digits=1), "-invgam2-", invgam2, "-H-", vec2code(diag(H), 100), "-h-", vec2code(hs, 100), "-tps-", idx.tps, "-xis-", vec2code(xis, 10), "-nSimu-", nSimu, ".RData")
save(post.res, paras, file=sv.name)


#H <- diag(c(0.5, 0.5, 0.05, 0.05))
#tXs <- matrix(c(0, 0, 1.9, -1.0, 
#                0, 0, 1.0, -1.0, 
#                1, 0, 1.0, 2.0, 
#                1, 1, 1.2, 1.3), ncol=4, byrow=T);tXs
#mMKF(c(0, 0, 1.0, -1.0), tXs, H)
#
#ttXs <- gen.Real.Xs(1000, fHats)
