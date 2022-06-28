setwd("C:/Users/JINHU/OneDrive - connect.hku.hk/文档/ProjectCode/HistTrial")
# 1. reduce the size of the dataset

#rm(list=ls())
source("util_ana.R")
#source("./realData/simuParas.R")
library(TeachingDemos)

idxs <- seq(1, 10, 1)
#H1
H1Dir <- "RealDataPureResampleAlp-b-3.5e-01-N-400-lam-300-lamq-10-phi0-3e-01-invgam2-0.33-H-010010999999-h-110110130130-tps-22cc-nSimu-2000"

filsH1 <- dir(paste0("./results/", H1Dir), pattern="[0-9].RData", full.names=1)
load(filsH1[1])

kpIdx <- !sapply(curRes, is.null)
Ns <- which(kpIdx)

OneResFn <- function(curRes){
    Z1s <- c()
    Z1s.no <- c()
    trtEffs <- c()
    trtEffs.no <- c()
    for (N in Ns){
        Z1r <- mean(curRes[[N]]$data$Z[21:N])
        Z1r.no <- mean(curRes[[N]]$data.no$Z[21:N])
        Z1s <- c(Z1s, Z1r)
        Z1s.no <- c(Z1s.no, Z1r.no)
        trtEffs <- c(trtEffs, curRes[[N]]$mtrt[1])
        trtEffs.no <- c(trtEffs.no, curRes[[N]]$mtrt[2])
    }
    list(
        Z1s=Z1s, 
        Z1s.no=Z1s.no, 
         trtEffs=trtEffs, 
         trtEffs.no=trtEffs.no
    )
}

obResFn <- function(filsH, b=0){
    flag <- 1
    Z1ss <- list()
    trtEffss <- list()
    Z1ss.no <- list()
    trtEffss.no <- list()
    for (filH in filsH){
        load(filH)
        res <- OneResFn(curRes)
        Z1ss[[flag]] <- res$Z1s
        Z1ss.no[[flag]] <- res$Z1s.no
        trtEffss[[flag]] <- abs(res$trtEffs-b)
        trtEffss.no[[flag]] <- abs(res$trtEffs.no-b)
        flag <- flag + 1
    }
    
    Lst2MeanFn <- function(curLst){
        curLstMat <- as.matrix(do.call(rbind, curLst))
        ms <- apply(curLstMat, 2, mean) 
        vs <- apply(curLstMat, 2, var) 
        list(ms=ms, vs=vs)
    }
    
    fres <- list(
        Z1=Lst2MeanFn(Z1ss)$ms, 
        Z1.no=Lst2MeanFn(Z1ss.no)$ms, 
        trtEff=Lst2MeanFn(trtEffss),
        trtEff.no=Lst2MeanFn(trtEffss.no)
                )
}

res <- obResFn(filsH1, b)



res$trtEff$Ls <- res$trtEff$ms - 1.96*sqrt(res$trtEff$vs)/sqrt(length(filsH1)*1)
res$trtEff$Us <- res$trtEff$ms + 1.96*sqrt(res$trtEff$vs)/sqrt(length(filsH1)*1)

jpeg("./plots/RealData_AbsBias.jpg", width=6, height=6, unit="in", res=500)
plot(Ns[idxs], res$trtEff$ms[idxs], type="l", col=1, lty=1, main="(a)",
     ylim=c(0.025, 0.12), ylab = "", xlab="N", lwd=3, cex.lab=1.5, cex.main=1.5)

ylab = expression('|' ~ hat(delta) ~ '-' ~ delta ~ '|')
title(ylab=ylab, line=2.2,  cex.lab=1.5)
polygon(c(Ns[idxs], rev(Ns[idxs])), c(res$trtEff$Ls[idxs], rev(res$trtEff$Us[idxs])), border=NA, col=rgb(1, 0, 0, 0.3))
lines(Ns[idxs], res$trtEff.no$ms[idxs], type="l", col="blue", lty=2, lwd=3)
legend("topright", c("BHCA", "KBCD"), col=c(1, "blue"), lty=1:2, lwd=3, cex=1.5)
tmpF <- function(){
    tmpidxs <- 7:10
plot(Ns[tmpidxs], res$trtEff$ms[tmpidxs], type="l", col=1, lty=1, 
      ylab = "", xlab="N", lwd=2, cex.lab=1.5, cex.main=1.5)
polygon(c(Ns[tmpidxs], rev(Ns[tmpidxs])), c(res$trtEff$Ls[tmpidxs], rev(res$trtEff$Us[tmpidxs])), border=NA, col=rgb(1, 0, 0, 0.3))
lines(Ns[tmpidxs], res$trtEff.no$ms[tmpidxs], type="l", col="blue", lty=2, lwd=2)
}

subplot( 
    tmpF(),
    x=grconvertX(c(0.5,1), from='npc'),
    y=grconvertY(c(0.2,0.7), from='npc'),
    type='fig', 
    pars=list( mar=c(1.5,1.5,0,0)+0.1) )
dev.off()



# fit1 <- smooth.spline(trtAbsMs[, 1], Ns)
# fit2 <- smooth.spline(trtAbsMs[, 2], Ns)
# predict(fit1, x=0.045)
# predict(fit2, x=0.045)


jpeg("./plots/RealData_Z1rate.jpg", width=6, height=6, unit="in", res=500)
plot(Ns[idxs], res$Z1[idxs], type="l", ylim=c(0.5, 0.6), col=1, lty=1, main="(c)",
     ylab="Allocation ratio", xlab="N", lwd=3, cex.lab=1.5, cex.main=1.5)
lines(Ns[idxs], res$Z1.no[idxs], type="l", ylim=c(0.5, 0.6), col=c("blue"), lty=2, lwd=3)
legend("topleft", c("BHCA", "KBCD"), col=c(1, "Blue"), lty=1:2, lwd=3, cex=1.5)
dev.off()

# test
# H0
H0Dir <- "RealDataPureResampleAlp-b-0e+00-N-400-lam-300-lamq-10-phi0-3e-01-invgam2-0.33-H-010010999999-h-110110130130-tps-22cc-nSimu-2000"
filsH0 <- dir(paste0("./results/", H0Dir), pattern="[0-9].RData", full.names=1)

post.b.fn2 <- function(trtsMat, dlt0=0){
    #trtsMat <- abs(trtsMat)
    prb <- post.fn(trtsMat[, 1], dlt0)
    prb.no <- post.fn(trtsMat[, 2], dlt0)
    c(prb, prb.no)
}

prbsFn <- function(filsH){
    prbss <- list()
    prbss.no <- list()
    flag <- 1
    for (filH in filsH){
        load(filH)
        prbs <- c()
        prbs.no <- c()
        for (N in Ns){
            probss <- post.b.fn2(curRes[[N]]$sps.trts, dlt0=0)
            prbs <- c(prbs, probss[1])
            prbs.no <- c(prbs.no, probss[2])
        }
        
        prbss[[flag]] <- prbs
        prbss.no[[flag]] <- prbs.no
        flag <- flag+1
    }
    list(prb=as.matrix(do.call(rbind, prbss)), 
         prb.no=as.matrix(do.call(rbind, prbss.no)))
    
}



#  test method: P(dlt>0|D) > dlt0, dlt0 is selected such that size is 0.05 
prbH0 <- prbsFn(filsH0)
epss.no <- list()
epss <- list()
flag <- 1
for (ix in 1:length(Ns)){
    eps0.no <- eps.alig(prbH0$prb.no[, ix], 0.05)
    eps0 <- eps.alig(prbH0$prb[, ix], 0.05)
    epss[[flag]] <- eps0
    epss.no[[flag]] <- eps0.no
    flag <- flag + 1
}
do.call(rbind, epss)
do.call(rbind, epss.no)


# calculate the power
prbH1 <- prbsFn(filsH1)
powers <- c()
powers.no <- c()
for (ix in 1:length(Ns)){
    eps0 <- epss[[ix]][1]
    eps0.no <- epss.no[[ix]][1]
    powers <- c(powers, mean(prbH1$prb[, ix]>eps0))
    powers.no <- c(powers.no, mean(prbH1$prb.no[, ix]>eps0.no))
}

jpeg("./plots/RealData_Power.jpg", width=6, height=6, unit="in", res=500)
par(mfrow=c(1, 1))
plot(Ns[idxs], powers[idxs], type='l', col=1, ylim=c(0.40, 0.98), ylab="Power", xlab="N", 
     lwd=3, main="(b)", cex.main=1.5, cex.lab=1.5)
lines(Ns[idxs], powers.no[idxs], type='l', col="blue", lty=2, lwd=3)
legend("bottomright", c("BHCA", "KBCD"), col=c(1, "blue"), lty=1:2, lwd=3, cex=1.5)
dev.off()


rbind(Ns, powers, powers.no)
