#rm(list=ls())
library(magrittr)
library(dplyr)
setwd("/home/huaqingj/MyResearch/HistTrial")
#setwd("/home/r13user3/MyResearch/HistTrial")
#setwd("C:/Users/JINHU/Documents/ProjectCode/HistTrial")
source("utils.R")
source("simplex.R")
library(parallel)
seeds <- 1:10000

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

fun.test <- function(i){
    #set.seed(seeds[i])
    alpss <- betass
    for (jj in 1:4){
        alpss[[jj]] <- betass[[jj]] + rnorm(p+1, mean=xis[jj], sd=ifelse(xis[jj]==0, 0, 1))
        #alpss[[jj]] <- betass[[jj]] + rnorm(p+1, sd=xis[jj])
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
    lam.trus <- rep(-1, N+1) # record lambda
    
    for (j in (n0+1):N){
        #print(j)
        cx <- unlist(gen.Data.Xs(1, x.tps))
        
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
    #print(trt.eff.no)
    #print(mean(data.no$Y))

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
xi.vs <- c(seq(0, 1, 0.2), 1.3, 1.6, 1.9, 2, 3, 5, 7, 9)
xis <- c(0, 3, 3, 0)
xiss <- list()
for (ii in 1:length(xi.vs)){
    xiss[[ii]] <- c(0, 0, 0, 0) + xi.vs[ii]
}

invgam2 <- 50

b <- 2
phi0 = phi1 = 3
N <- 150 # total sample size
# parameters
lam.q <- 0.10 
lam <- 30
hs <- rep(1.1, 4)
#x.tps <- c(2, 2, "c", "c")
x.tps <- c(2, 2, 2, 2)
idx.tps <- paste0(x.tps, collapse="")
p <- length(x.tps)

#H <- diag(rep(1, p)/2)
H <- diag(rep(0.05, p))

# initial dataset
n0 <- 20

# for hypothesis test
# to calculate the prob
M <- 1000

nSimu <- 1000
for (invgam2 in c(100)) {
for (phi0 in c(3)){
#for (xis.idx in 1:length(xiss)) {
        #xis <- xiss[[xis.idx]]
        phi1 = phi0
        paras <- list(invgam2=invgam2, b=b, phi0=phi0, phi1=phi1, N=N, lam.q=lam.q, hs=hs, x.tps=x.tps, H=H, M=M, xis=xis, lam=lam)
        post.res <- mclapply(1:nSimu, fun.test, mc.cores=20)
        sv.name <- paste0("./results/LinearSwNewAlp", "-b-", b, "-N-", N, "-lam-", lam, "-lamq-", 100*lam.q, "-phi0-", phi0, "-invgam2-", invgam2, "-H-", vec2code(diag(H), 100), "-h-", vec2code(hs, 100), "-tps-", idx.tps, "-xis-", vec2code(xis, 10), "-nSimu-", nSimu, ".RData")
        save(post.res, paras, file=sv.name)
#}
}
}
