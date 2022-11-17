setwd("C:/Users/JINHU/OneDrive - connect.hku.hk/文档/ProjectCode/HistTrial")

#rm(list=ls())
source("util_ana.R")
#source("./realData/simuParas.R")

#H1
H1Dirs <- dir("./results/", pattern="RealDataPureResampleHistSS-b-3.5e-01*", full.names=1)
H1Dirs <- sort(H1Dirs);
H1Dirs <- H1Dirs[c(5, 8, 1:4, 6, 7)];H1Dirs
H0Dirs <- dir("./results/", pattern="RealDataPureResampleHistSS-b-0e.*", full.names=1)
H0Dirs <- sort(H0Dirs);
H0Dirs <- H0Dirs[c(5, 8, 1:4, 6, 7)];H0Dirs



OneResFn <- function(curRes){
    N <- dim(curRes$data)[1]
    Z1 <- mean(curRes$data$Z[21:N])
    Z1.no <- mean(curRes$data.no$Z[21:N])
    trtEff <- curRes$mtrt[1]
    trtEff.no <- curRes$mtrt[2]
    list(
        Z1=Z1, 
        Z1.no=Z1.no, 
        trtEff=trtEff, 
        trtEff.no=trtEff.no
    )
}

obResFn <- function(filsH, b=0){
    Z1s <- c()
    trtEffs <- c()
    Z1s.no <- c()
    trtEffs.no <- c()
    for (filH in filsH){
        load(filH)
        res <- OneResFn(curRes)
        Z1s <- c(Z1s, res$Z1)
        Z1s.no <- c(Z1s.no, res$Z1.no)
        trtEffs <- c(trtEffs, abs(res$trtEff-b))
        trtEffs.no <- c(trtEffs.no, abs(res$trtEff.no-b))
    }
    
    
    fres <- list(
        Z1=mean(Z1s),
        Z1.sd=sd(Z1s),
        Z1.no=mean(Z1s.no),
        trtEff=mean(trtEffs),
        trtEff.no=mean(trtEffs.no), 
        trtEff.sd =sd(trtEffs),
        trtEff.no.sd=sd(trtEffs.no)
    )
}


N0s <- c(30, seq(50, 350, 50));N0s
tRes <- list()
for (ix in 1:length(H1Dirs)) {
    print(ix)
    filH1s <- dir(H1Dirs[ix], pattern="[0-9].*", full.names = 1)
    res <- obResFn(filH1s, b=b)
    tRes$trtms <- c(tRes$trtms, res$trtEff)
    tRes$trtms.no <- c(tRes$trtms.no, res$trtEff.no)
    tRes$Z1s <- c(tRes$Z1s, res$Z1)
    tRes$Z1sLs <-  c(tRes$Z1sLs, res$Z1- 1.96 * res$Z1.sd/sqrt(length(filH1s)))
    tRes$Z1sUs <-  c(tRes$Z1sUs, res$Z1+ 1.96 * res$Z1.sd/sqrt(length(filH1s)))
    tRes$Z1s.no <- c(tRes$Z1s.no, res$Z1.no)
    tRes$trtLs <-  c(tRes$trtLs, res$trtEff - 1.96 * res$trtEff.sd/sqrt(length(filH1s)))
    tRes$trtUs <-  c(tRes$trtUs, res$trtEff + 1.96 * res$trtEff.sd/sqrt(length(filH1s)))
    
}

jpeg("./plots/RealDataN0_AbsBias.jpg", width=6, height=6, unit="in", res=500)
plot(N0s, tRes$trtms, type="l", col=1, lty=1, main="(a)",
     ylim=c(0.035, 0.05), ylab = "", xlab=expression(N[0]), lwd=3, cex.lab=1.5, cex.main=1.5)
ylab = expression('|' ~ hat(delta) ~ '-' ~ delta ~ '|')
title(ylab=ylab, line=2.2,  cex.lab=1.5)
polygon(c(N0s, rev(N0s)), c(tRes$trtLs, rev(tRes$trtUs)), border=NA, col=rgb(1, 0, 0, 0.3))
abline(h=mean(tRes$trtms.no), col="blue", lty=2, lwd=3)
legend("topright", c("CAHB", "KBCD"), col=c(1, "blue"), lty=1:2, lwd=3, cex=1.5)
dev.off()

jpeg("./plots/RealDataN0_Z1rate.jpg", width=6, height=6, unit="in", res=500)
plot(N0s, tRes$Z1s, type="l", ylim=c(0.562, 0.575), col=1, lty=1, main="(c)", 
     ylab="Allocation ratio", xlab=expression(N[0]), lwd=3, cex.main=1.5, cex.lab=1.5)
#polygon(c(N0s, rev(N0s)), c(tRes$Z1sLs, rev(tRes$Z1sUs)), border=NA, col=rgb(1, 0, 0, 0.3))
abline(h=mean(tRes$Z1s.no), col="blue", lty=2, lwd=3)
#legend("topleft", c("CAHB", "KBCD"), col=c("Red", "Blue"), lty=1:2, lwd=2)
dev.off()



post.b.fn2 <- function(trtsMat, dlt0=0){
    #trtsMat <- abs(trtsMat)
    prb <- post.fn(trtsMat[, 1], dlt0)
    prb.no <- post.fn(trtsMat[, 2], dlt0)
    c(prb, prb.no)
}

prbsFn <- function(filsH){
    prbs <- c()
    prbs.no <- c()
    flag <- 1
    for (filH in filsH){
        load(filH)
        probs <- post.b.fn2(curRes$sps.trts, dlt0=0)
        prb <- probs[1]
        prb.no <- probs[2]
        
        prbs[flag] <- prb
        prbs.no[flag] <- prb.no
        flag <- flag+1
    }
    list(prb=prbs, 
         prb.no=prbs.no)
    
}



#  test method: P(dlt>0|D) > dlt0, dlt0 is selected such that size is 0.05 
epss.no <- list()
epss <- list()
flag <- 1
for (ix in 1:length(H0Dirs)) {
    print(ix)
    filH0s <- dir(H0Dirs[ix], pattern="[0-9].*", full.names = 1)
    prbH0 <- prbsFn(filH0s)
    eps0.no <- eps.alig(prbH0$prb.no, 0.05)
    eps0 <- eps.alig(prbH0$prb, 0.05)
    epss[[flag]] <- eps0
    epss.no[[flag]] <- eps0.no
    flag <- flag + 1
    
}
do.call(rbind, epss)
do.call(rbind, epss.no)


# calculate the power
powers <- c()
powers.no <- c()
for (ix in 1:length(H1Dirs)){
    print(ix)
    filH1s <- dir(H1Dirs[ix], pattern="[0-9].*", full.names = 1)
    prbH1 <- prbsFn(filH1s)
    eps0 <- epss[[ix]][1]
    eps0.no <- epss.no[[ix]][1]
    powers <- c(powers, mean(prbH1$prb>eps0))
    powers.no <- c(powers.no, mean(prbH1$prb.no>eps0.no))
}

jpeg("./plots/RealDataN0_Power.jpg", width=6, height=6, unit="in", res=500)
par(mfrow=c(1, 1))
plot(N0s, powers, type='l', col=1, main="(b)", 
     ylim=c(0.85, 0.95), ylab="Power", xlab=expression(N[0]), lwd=3, cex.lab=1.5, cex.main=1.5)
abline(h=mean(powers.no), col="blue", lty=2, lwd=3)
legend("topleft", c("CAHB", "KBCD"), col=c(1, "blue"), lty=1:2, lwd=3, cex=1.5)
dev.off()

