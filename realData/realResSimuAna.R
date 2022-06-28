setwd("C:/Users/JINHU/OneDrive - connect.hku.hk/文档/ProjectCode/HistTrial")
source("util_ana.R")

ext.sub.info.fn <- function(res){
        subId <- res$data$frx*2 + res$data$falls+ 1
        res.info <- list(Z=res$data$Z, subId=subId, tau2s=res$tau2s)
        Z1ms <- sapply(1:4, function(i) {mean(res.info$Z[21:N][res.info$subId[21:N]==i])})
        tau2ms <-  sapply(1:4, function(i) {mean(res.info$tau2s[res.info$subId==i])})
        rv <- list(Z1ms=Z1ms, tau2ms=tau2ms)
        rv
}   

post.b.fn2 <- function(trtsMat, dlt0=0){
    #trtsMat <- abs(trtsMat)
    prb <- post.fn(trtsMat[, 1], dlt0)
    prb.no <- post.fn(trtsMat[, 2], dlt0)
    c(prb, prb.no)
    
}


N <- 200
is.subgrp <- F

# all results under no borrow, 
load("./results/realResSimuAnaNoBorrow.RData")
load("./results/realResSimuAnaNoBorrowSubGrp.RData")
noResAll <- list(powers=c(noRes$powers, noRes.sub$powers), 
                 Z1s=c(noRes$Z1s, noRes.sub$Z1s), 
                 errs=c(noRes$errs, noRes.sub$errs))

if (is.subgrp){
    fil <-  dir(paste0("./results/RealDataAgeSimuH1Subgrp", N), pattern="RealData.*RData", full.names=T);fil
}else{
    fil <-  dir(paste0("./results/RealDataAgeSimuH1", N), pattern="RealData.*RData", full.names=T);fil
}
fil <- sort(fil);fil
flag <- 1
Hs <- c()
Z1r.list <- list()
trt.Err.list <- list()
Z1Mats <- list()
tauMats <- list()
for (tf in fil){
    cH <- (as.numeric(strsplit(tf, "-")[[1]][23])/100) %% 100/ 10
    Hs <- c(Hs, cH)
    load(tf)
    
    kpidxs <- sapply(post.res, function(res)length(res$res$phi0.tk)) != 100
    post.res <- post.res[kpidxs];length(post.res)/length(kpidxs)
    
    varYs <- sapply(post.res, function(s.res)var(s.res$data$Y))
    Z1rates <- sapply(post.res, function(s.res)mean(s.res$data$Z[21:N]))
    Z1rates.no <- sapply(post.res, function(s.res)mean(s.res$data.no$Z[21:N]))
    Z1RatesCI <- rbind(CI.fn(Z1rates), CI.fn(Z1rates.no))
    rownames(Z1RatesCI) <- c("Borrow", "No Borrow")
    Z1RatesCI
    
    trtMats <- sapply(post.res, function(x)x$mtrt)
    errs <- abs(trtMats - b)
    trt.Errs <- rbind(CI.fn(errs[1, ]), CI.fn(errs[2, ]))
    rownames(trt.Errs) <- c("Borrow", "No Borrow")
    trt.Errs
    
    Z1mss <- lapply(post.res, function(res)ext.sub.info.fn(res$res)$Z1ms)
    Z1Mat <- do.call(rbind, Z1mss)
    
    taumss <- lapply(post.res, function(res)ext.sub.info.fn(res$res)$tau2ms)
    tauMat <- do.call(rbind, taumss)
        
    Z1r.list[[flag]] <- Z1RatesCI
    trt.Err.list[[flag]] <- trt.Errs
    Z1Mats[[flag]] <- apply(Z1Mat, 2, mean)
    tauMats[[flag]] <- apply(tauMat, 2, mean)
    flag <- flag + 1
}

Z1s <- sapply(Z1r.list, function(x)x[1, 2])
Z1s.no <- sapply(Z1r.list, function(x)x[2, 2])
errs.no <- sapply(trt.Err.list, function(x)return(x[2, 2]))
errs.b <- sapply(trt.Err.list, function(x)return(x[1, 2]))
errs.low  <- sapply(trt.Err.list, function(x)return(x[1, 1]))
errs.up <- sapply(trt.Err.list, function(x)return(x[1, 3]))



if (is.subgrp){
    jpeg("./plots/RealDataSimu_Subgroup_AbsBias.jpg", width=6, height=6, unit="in", res=500)
    par(mfrow=c(1, 1))
    plot(Hs, errs.b, type="l", ylim=c(0.035, 0.048), xlab=expression(xi), ylab="",
         lty=1, col=1, lwd=3, main="(a)",  cex.lab=1.5, cex.main=1.5)
    ylab = expression('|' ~ hat(delta) ~ '-' ~ delta ~ '|')
    title(ylab=ylab, line=2.2,  cex.lab=1.5)
    polygon(c(Hs, rev(Hs)), c(errs.low, rev(errs.up)), col=rgb(1, 0, 0, 0.3), border=NA)
    abline(h=mean(noResAll$errs), col="blue", lty=2, lwd=3)
    legend("topright", c("BHCA", "KBCD"), col=c(1, "blue"), lty=c(1, 2), lwd=3, cex=1.5)
    dev.off()
    
    jpeg("./plots/RealDataSimu_Subgroup_rates.jpg", width=6, height=6, unit="in", res=500)
    Z1MatsArr <- do.call(rbind, Z1Mats)
    plot(Hs, Z1MatsArr[, 1], ylim=c(0.47, 0.60), col=4, lty=2, type="l", lwd=3, xlab=expression(xi), 
         ylab="Allocation ratio", main="(c)",  cex.lab=1.5, cex.main=1.5)
    lines(Hs, Z1MatsArr[, 2], col=5, lty=3, lwd=3)
    lines(Hs, Z1MatsArr[, 3], col=6, lty=4, lwd=3)
    lines(Hs, Z1MatsArr[, 4], col=7, lty=5, lwd=3)
    lines(Hs, Z1s, col=1, lty=1, lwd=4)
    #lines(Hs, Z1s.no, col="blue", lty=2, lwd=2)
    abline(h=mean(noResAll$Z1s), col="blue", lty=2, lwd=3)
    #abline(h=0.5, col=6, lty=2, lwd=1)
    #axis(side=2, at=0.5, label="0.5")
    legend("topright", c(paste0("Subgroup ", 1:4), "Overall"), 
           col=c(4:7, 1), lty=c(2:5, 1), lwd=c(3, 3, 3, 3, 4), ncol=2, cex=1.2)
    dev.off()
    
    jpeg("./plots/RealDataSimu_Subgroup_Tau2.jpg", width=6, height=6, unit="in", res=500)
    tauMatsArr <- do.call(rbind, tauMats)
    plot(Hs, tauMatsArr[, 1], ylim=c(0, 7.2), col=1, lty=1, type="l", lwd=2, xlab=expression(xi), 
         ylab=expression(tau^2))
    lines(Hs, tauMatsArr[, 2], col=2, lty=2, lwd=2)
    lines(Hs, tauMatsArr[, 3], col=3, lty=3, lwd=2)
    lines(Hs, tauMatsArr[, 4], col=4, lty=4, lwd=2)
    abline(h=0, col=6, lty=2, lwd=1)
    legend("topright", paste0("Subgroup", 1:4), col=1:4, lty=1:4, lwd=c(2,2,2,2), ncol=1)
    dev.off()

    
}else{
    jpeg("./plots/RealDataSimu_AbsBias.jpg", width=6, height=6, unit="in", res=500)
    par(mfrow=c(1, 1))
    plot(Hs, errs.b, type="l", ylim=c(0.035, 0.048), xlab=expression(xi), ylab="", main="(a)",
         lty=1, col=1, lwd=3, cex.lab=1.5, cex.main=1.5)
    ylab = expression('|' ~ hat(delta) ~ '-' ~ delta ~ '|')
    title(ylab=ylab, line=2.2,  cex.lab=1.5)
    polygon(c(Hs, rev(Hs)), c(errs.low, rev(errs.up)), col=rgb(1, 0, 0, 0.3), border=NA)
    abline(h=mean(noResAll$errs), col="blue", lty=2, lwd=3)
    legend("topright", c("BHCA", "KBCD"), col=c(1, "blue"), lty=c(1, 2), cex=1.5, lwd=3)
    dev.off()
    
    jpeg("./plots/RealDataSimu_rates.jpg", width=6, height=6, unit="in", res=500)
    Z1MatsArr <- do.call(rbind, Z1Mats)
    plot(Hs, Z1MatsArr[, 1], ylim=c(0.47, 0.60), col=4, lty=2, type="l", lwd=3, xlab=expression(xi), 
         ylab="Allocation ratio", main="(c)",  cex.lab=1.5, cex.main=1.5)
    lines(Hs, Z1MatsArr[, 2], col=5, lty=3, lwd=3)
    lines(Hs, Z1MatsArr[, 3], col=6, lty=4, lwd=3)
    lines(Hs, Z1MatsArr[, 4], col=7, lty=5, lwd=3)
    lines(Hs, Z1s, col=1, lty=1, lwd=4)
    #lines(Hs, Z1s.no, col="blue", lty=2, lwd=2)
    abline(h=mean(noResAll$Z1s), col="blue", lty=2, lwd=3)
    #abline(h=0.5, col=6, lty=2, lwd=1)
    #axis(side=2, at=0.5, label="0.5")
    legend("topright", c(paste0("Subgroup ", 1:4), "Overall"), 
           col=c(4:7, 1), lty=c(2:5, 1), lwd=c(3, 3, 3, 3, 4), ncol=2, cex=1.2)
    dev.off()
}

# powers 
if (is.subgrp){
    filH0 <-  dir(paste0("./results/RealDataAgeSimuH0Subgrp", N), pattern="RealData.*RData", full.names=T);filH0
}else{
    filH0 <-  dir(paste0("./results/RealDataAgeSimuH0", N), pattern="RealData.*RData", full.names=T);filH0
}
filH0 <- sort(filH0);filH0
Hs <- c()
eps0s <- c()
eps0s.no <- c()
for (curFil in filH0){
    cH <- (as.numeric(strsplit(curFil, "-")[[1]][22])/100) %% 100/ 10
    Hs <- c(Hs, cH)
    load(curFil)
    probMat <- sapply(post.res, function(res)post.b.fn2(res$sps.trts, dlt0=0))
    eps0 <- eps.alig(probMat[1, ], 0.05)
    eps0s <- rbind(eps0s, unlist(eps0))
    eps0.no <- eps.alig(probMat[2, ], 0.05)
    eps0s.no <- rbind(eps0s.no, unlist(eps0.no))
}



flag <- 1
powers <- c()
powers.no <- c()
for (tf in fil){
    load(tf)
    probMat <- sapply(post.res, function(res)post.b.fn2(res$sps.trts, dlt0=0))
    eps0 <- eps0s[flag, 1]
    eps0.no <- eps0s.no[flag, 1]
    powers <- c(powers, mean(probMat[1, ]>eps0))
    powers.no <- c(powers.no, mean(probMat[2, ]>eps0.no))
    flag <- flag + 1
}

if (is.subgrp){
jpeg("./plots/RealDataSimu_Subgroup_Power.jpg", width=6, height=6, unit="in", res=500)
}else{
jpeg("./plots/RealDataSimu_Power.jpg", width=6, height=6, unit="in", res=500)
}
plot(Hs, powers, type="l", ylim=c(0.80, 0.95), xlab=expression(xi), ylab="Power", main="(b)",
     lty=1, col=1, lwd=3, cex.main=1.5, cex.lab=1.5)
abline(h=mean(noResAll$powers), col="blue", lty=2, lwd=3)
legend("topright", c("BHCA", "KBCD"), col=c(1, "blue"), lty=c(1, 2), lwd=3, cex=1.5)
dev.off()

if (is.subgrp){
    noRes.sub <- list(powers=powers.no, Z1s=Z1s.no, errs=errs.no)
    save(noRes.sub, file="./results/realResSimuAnaNoBorrowSubGrp.RData")
    
}else{
    noRes <- list(powers=powers.no, Z1s=Z1s.no, errs=errs.no)
    save(noRes, file="./results/realResSimuAnaNoBorrow.RData")
}