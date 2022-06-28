# this script is to run real data when N is changed
#rm(list=ls())
library(magrittr)
library(dplyr)
setwd("/root/Documents/HQ/HistTrial")
#setwd("/home/huaqingj/MyResearch/HistTrial")
#setwd("/home/r13user3/MyResearch/HistTrial")
#setwd("C:/Users/JINHU/Documents/ProjectCode/HistTrial")
source("utils.R")
source("simplex.R")
source("./realData/simuParas.R")
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

alpss.fn <- function(curDat){
    curDat$Y01 <-as.vector(scale(curDat$Y0, center=mean.Y0, scale=sd.Y0))
    curDat$Y1 <-as.vector(scale(curDat$Y, center=mean.Y, scale=sd.Y))
    curDat$menyrs1 <- as.vector(scale(curDat$menyrs, center=mean.menyrs, scale=sd.menyrs))
    curFit<- lm(Y1~Z+subGroupId+subGroupId*menyrs1+subGroupId*Y01, data=curDat)
    coeffsH <- curFit$coefficients
    parasH.Vec <- coeffsH[-2][c(1:4, 5, 7:9, 6, 10:12)]
    eparasH <- matrix(parasH.Vec, ncol=4, byrow=T)
    eparasH[1, 2:4] <- eparasH[1, 2:4] + eparasH[1, 1]
    eparasH[2, 2:4] <- eparasH[2, 2:4] + eparasH[2, 1]
    eparasH[3, 2:4] <- eparasH[3, 2:4] + eparasH[3, 1]
    row.names(eparasH) <- c("Intercept", "Menyrs", "dxhhp_0")
    alpssMat <- t(rbind(eparasH[1, ], rep(0, 4), rep(0, 4), eparasH[c(2, 3), ]))
    colnames(alpssMat) <- c("intercept", "falls", "frx", "menyrs", "dxhhp_0")
    alpss <- list(
       para1=alpssMat[1, ],
       para2=alpssMat[2, ],
       para3=alpssMat[3, ],
       para4=alpssMat[4, ]
    )
    return(alpss)
}

fun.real <- function(i){
    rvs <- list()
    ts <- c()
    ts <- c(ts, Sys.time())
    set.seed(seeds[i])
    is.AllSub <- FALSE
    while (!is.AllSub){
        kpIdx <- sample.int(dim(data.Hist)[1], N0)
        curDat <- data.Hist[kpIdx, ]
        is.AllSub <- sum(table(curDat$subGroupId) <= 2) == 0
        #print(table(curDat$subGroupId))
    }
    alpss <- alpss.fn(curDat)
    print(c(N0, i))
    Xs <- gen.Real.Xs(n0, fHats, subGrpDist)
    
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
    ts <- c(ts, Sys.time())
    
    lastRes <- NA
    for (j in (n0+1):N){
        #print(j)
        t0 <- Sys.time()
        cx <- unlist(gen.Real.Xs(1, fHats, subGrpDist))
        
        
        H <- diag(c(0.10, 0.10, bw.nrd(data$menyr), bw.nrd(data$Y0)))

        alpMat <- sub.Paras.fn(Xs, alpss)
        Theta0s <- curMean.fn(Xs, Zs, alpMat, b=0)
        lam.tru <- lam.sel.fn(data, H=H, invgam2=invgam2, lam=lam, lam.q=lam.q)
        lam.trus[j] <- lam.tru
        t1 <- Sys.time()
        res <- mu0.info.est.fn(Theta0s, data, H, lam, invgam2=invgam2, lam.tru=lam.tru, lastRes=lastRes)
        lastRes <- res
        t2 <- Sys.time()
        ts <- c(t0, t1, t2)
        #print(diff(ts))
        #print(c(j, res$iterN))

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
        
        if (j==N){
            alpMat <- sub.Paras.fn(Xs, alpss)
            Theta0s <- curMean.fn(Xs, Zs, alpMat, b=0)
            lam.tru <- lam.sel.fn(data, H=H, invgam2=invgam2, lam=lam, lam.q=lam.q)
            res <- mu0.info.est.fn(Theta0s, data, H, lam, invgam2=invgam2, lam.tru=lam.tru, lastRes=lastRes)
            res.no <- mu0.no.est.fn(data.no, H)
            
            res.mu1 <- mu1.no.est.fn(data, H)
            res.no.mu1 <- mu1.no.est.fn(data.no, H)
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

            
            rvs <- list(mtrt=c(trt.eff, trt.eff.no), sps.trts=sps.trts, data=data, data.no=data.no, res=res, Rs=Rs, lam.trus=lam.trus,alpss=alpss, betass=betass)
        }
        
    }
    rvs
}

#under H0
b<-0

invgam2 <- 0.33

phi0 = phi1 = sd(fit.all$residuals)
N <- 200 # total sample size
# parameters
lam.q <- 0.10
lam <- 300
hs <- c(1.1, 1.1, 1.3, 1.3)
x.tps <- c(2, 2, "c", "c")
idx.tps <- paste0(x.tps, collapse="")
p <- length(x.tps)


# initial dataset
n0 <- 20
# the sample size to estimate the parameters for historical data
#N0s <- c(seq(50, 350, 50))
N0s <- c(30, seq(50, 350, 50))
#N0s <- c(20, seq(50, 350, 50))

# for hypothesis test
# to calculate the prob
M <- 1000

nSimu <- 2000
for (N0 in N0s){

post.res <- mclapply(1:nSimu, fun.real, mc.cores=18)
H <- diag(c(0.10, 0.10, 9.99, 9.99))
folder.name <- paste0("./results/RealDataPureResampleHistSS", "-b-", format(b, scientific=T, digits=2), "-N-", N, "-N0-", N0, "-lam-", lam, "-lamq-", 100*lam.q, "-phi0-", format(phi0, scientific=T, digits=1), "-invgam2-", invgam2, "-H-", vec2code(diag(round(H, 2)), 100), "-h-", vec2code(hs, 100), "-tps-", idx.tps, "-nSimu-", nSimu)
if (!dir.exists(folder.name)){
    dir.create(folder.name)
}
paras <- list(invgam2=invgam2, b=b, phi0=phi0, phi1=phi1, N=N, lam.q=lam.q, hs=hs, x.tps=x.tps, H=H, M=M, lam=lam)
save(paras, file=paste0(folder.name, "/paras.RData"))
for (ix in 1:nSimu){
    curRes <- post.res[[ix]]
    save(curRes, file=paste0(folder.name, "/", ix, ".RData"))
}

}

