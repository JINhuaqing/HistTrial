# script for whole procedure
#rm(list=ls())
#library(magrittr)
#library(dplyr)
#setwd("C:/Users/JINHU/Documents/ProjectCode/HistTrial")
#source("utils.R")


# # parameters for current data 
# betass <- list(para1=c(2, 1, -1, 3, -2), 
#               para2=c(3, 1, 2, 4, 2), 
#               para3=c(0, -1, -1, 3, 2), 
#               para4=c(2, 0, -1, 4, -2)
#                   )
# 
# # parameters for historical data
# alpss <-  list(para1=c(2, 1, -1, 3, -2), 
#              para3=c(3, 1, 2, 4, 2), 
#              para2=c(0, -1, -1, 3, 2), # switch para2 and para3
#              para4=c(2, 0, -1, 4, -2) )
# b <- 2
# phi0 = phi1 = 1
# N <- 100 # total sample size
# # parameters
# lam <- 0.05
# hs <- rep(2.1, 4)
# 
# # initial dataset
# n0 <- 10
# x.tps <- c(2, 2, "c", "c")
Xs <- gen.Data.Xs(n0, x.tps)
idx0 <- sample.int(n0, size=floor(n0/2))
Zs <- rep(1, n0)
Zs[idx0] <- 0


betMat <- sub.Paras.fn(Xs, betass)
Ys <- curMean.fn(Xs, Zs, betMat, b) + rnorm(n0, sd=phi0)
data <- cbind(Ys, Zs, Xs)
data <- as.data.frame(data)
colnames(data)[1:2] <- c("Y", "Z")


# no borrowing 
Zs.no <- Zs
Ys.no <- curMean.fn(Xs, Zs.no, betMat, b) + rnorm(n0, sd=phi0)
data.no <- cbind(Ys.no, Zs.no, Xs)
data.no <- as.data.frame(data.no)
colnames(data)[1:2] <- c("Y", "Z")

for (j in (n0+1):N){
    cx <- unlist(gen.Data.Xs(1, x.tps))
    
    # H <- diag(c(bw.nrd(data$X1), bw.nrd(data$X2), bw.nrd(data$X3), bw.nrd(data$X4)))
    alpMat <- sub.Paras.fn(Xs, alpss)
    Theta0s <- curMean.fn(Xs, Zs, alpMat, b=0)
    res <- info.est.fn(Theta0s, data, H, lam)
    res0 <- info.est.fn(Theta0s, data, H, lam, is.borrow=FALSE)
    
    var.info <- post.var.mu0.fn(cx, res)
    var.ref <- post.var.mu0.fn(cx, res0)
    R <- var.ref/var.info
    ass.res <- RPS.design(cx, data[, 3:6], data$Z, hs, R)
    ass.res.no <- RPS.design(cx, data[, 3:6], data.no$Z, hs, R=1)
    # ass.res <- RBC.design(cx, data[, 3:6], data$Z, hs, R)
    
    Xs <- rbind(Xs, cx)
    Zs <- c(Zs, ass.res$grp-1)
    Zs.no <- c(Zs.no, ass.res.no$grp-1)
    
    curN <- dim(Xs)[1]
    curBetMat <- sub.Paras.fn(Xs[curN, ], betass)
    Y <- curMean.fn(Xs[curN, ], Zs[curN], curBetMat, b) + rnorm(1, sd=phi0)
    Ys <- c(Ys, Y)
    data <- cbind(Ys, Zs, Xs)
    data <- as.data.frame(data)
    colnames(data)[1:2] <- c("Y", "Z")
    
    Y.no <- curMean.fn(Xs[curN, ], Zs.no[curN], curBetMat, b) + rnorm(1, sd=phi0)
    Ys.no <- c(Ys.no, Y.no)
    data.no <- cbind(Ys.no, Zs.no, Xs)
    data.no <- as.data.frame(data.no)
    colnames(data.no)[1:2] <- c("Y", "Z")
    
    
    if (j==N){
        #H <- diag(c(bw.nrd(data$X1), bw.nrd(data$X2), bw.nrd(data$X3), bw.nrd(data$X4)))
        alpMat <- sub.Paras.fn(Xs, alpss)
        Theta0s <- curMean.fn(Xs, Zs, alpMat, b=0)
        res <- info.est.fn(Theta0s, data, H, lam)
        res.no <- info.est.fn(Theta0s, data.no, H, lam, is.borrow=FALSE)
    }
}


res.mu1 <- mu1.est.fn(data$Y, data, H)
res.no.mu1 <- mu1.est.fn(data.no$Y, data.no, H)
trt.eff <- mean(mu1.efn(as.matrix(Xs), res.mu1)) -  mean(mu0.efn(as.matrix(Xs), res));trt.eff
trt.eff.no <- mean(mu1.efn(as.matrix(Xs), res.no.mu1)) -  mean(mu0.efn(as.matrix(Xs), res.no));trt.eff.no

