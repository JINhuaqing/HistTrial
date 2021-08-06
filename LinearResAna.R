setwd("C:/Users/JINHU/Documents/ProjectCode/HistTrial")
library(magrittr)
source("util_ana.R")




#1. vary phi0
{

typ <- "Simplex-Linear-diff-0"
N <- 150
lam <- 20
b <- 4
invgam2 <- 1
nSimu  <- 1000
#allocation ratio vs Var(Y|X)/Var(Y)
phi0s <- c(1, 3, 5, 7, 9)
Z1Ms.CIs <- list()
Z1noMs <- c()
var.rateMs <- c()
varMs <- c()
varXMs <- c()
biass <- list()
i <- 1
for (phi0 in phi0s){
    f.name <- paste0("./results/Phi0s/", typ, "-b-", b, "-N-", N, 
    "-lam-",  lam, "-phi0-", phi0, "-invgam2-", invgam2, "-nSimu-", nSimu, ".RData")
    load(f.name)
    kpidxs <- sapply(post.res, function(res)length(res$res$phi0.tk)) != 100
    post.res <- post.res[kpidxs]
    varYs <- sapply(post.res, function(s.res)var(s.res$data$Y))
    varMs <- c(varMs, mean(varYs))
    Z1rates <- sapply(post.res, function(s.res)mean(s.res$data$Z))
    Z1rates.no <- sapply(post.res, function(s.res)mean(s.res$data.no$Z))
    varYXs <- phi0*2
    varXMs <- c(varXMs, mean(varYXs))
    var.rates <- varYXs/varYs
    Z1Ms.CIs[[i]] <- CI.fn(Z1rates)
    Z1noMs <- c(Z1noMs, mean(Z1rates.no))
    var.rateMs <- c(var.rateMs, mean(var.rates))
    
    trtMats <- sapply(post.res, function(x)x$mtrt)
    errs <- abs(trtMats - b)
    biass[[i]] <- rbind(CI.fn(errs[1, ]), CI.fn(errs[2, ]))
    
    i <- i+1
    
}

# plot about allocation rate
Z1MsMat <- do.call(rbind, Z1Ms.CIs)
png(filename = paste0("./plots/Z1rateVsphi0_", typ, "_phi0_", phi0, ".png"), width=360, height=480)
plot(phi0s, Z1MsMat[, 2], type="l", xlab="Var(Y|X, Z)", ylab="Percentage of treatment group", ylim=c(0.45, 0.95))
polygon(c(phi0s, rev(phi0s)), c(Z1MsMat[, 1], rev(Z1MsMat[, 3])),  col= rgb(0, 0, 0, 0.1), border=NA)
lines(phi0s, Z1noMs, lwd=2, col=2, lty=2)
legend("topleft", c("Borrow", "no Borrow"), col=c(1, 2), lty=c(1, 2))
dev.off()

png(filename = paste0("./plots/Z1rateVsrv", typ, "_phi0_", phi0, ".png"), width=360, height=480)
plot(var.rateMs, Z1MsMat[, 2], type="p", xlab="Var(Y|X, Z)/Var(Y)", ylab="Percentage of treatment group")
dev.off()

# plot about bias
png(filename = paste0("./plots/biasvsphi0", typ, "_phi0_", phi0, ".png"), width=360, height=480)
biassMat.b <- sapply(biass, function(x)x[1, ])
biassMat.no <- sapply(biass, function(x)x[2, ])
plot(phi0s, biassMat.b[2, ], type="l", lwd=2, col=1, xlab="Var(Y|X, Z)", ylab='Absolute bias', lty=1, ylim=c(0.15, 1.2))
polygon(c(phi0s, rev(phi0s)), c(biassMat.b[1, ], rev(biassMat.b[3, ])),  col= rgb(0, 0, 0, 0.1), border=NA)
lines(phi0s, biassMat.no[2, ], lwd=2, col=2, lty=2)
#polygon(c(phi0s, rev(phi0s)), c(biassMat.no[1, ], rev(biassMat.no[3, ])),  col= rgb(1, 0, 0,0.1), border=NA)
legend("topleft", c("Borrow", "no Borrow"), col=c(1, 2), lty=c(1, 2))
dev.off()
}


#2. Varying b 
{
# for hypothesis test
dlt0 <- 0.5

typ <- "Simplex-Linear-diff-0"
N <- 150
lam <- 20
b <- 0
phi0 <- 5
invgam2 <- 1
nSimu <- 2000
f.name <- paste0("./results/bs/", typ, "-b-", b, "-N-", N, 
    "-lam-",  lam, "-phi0-", phi0, "-invgam2-", invgam2, "-nSimu-", nSimu, ".RData")
load(f.name)
kpidxs <- sapply(post.res, function(res)length(res$res$phi0.tk)) != 100
post.res <- post.res[kpidxs]

# eps0 and eps0.no are tuned such that the size is 0.05 
probMat <- sapply(post.res, post.b.fn, dlt0=dlt0)
eps.alig(probMat[1, ])
eps.alig(probMat[2, ])
eps0 <- eps.alig(probMat[1, ])$eps0
eps0.no <- eps.alig(probMat[2, ])$eps0
mean(probMat[1, ]>eps0)
mean(probMat[2, ]>eps0.no)

bs <- seq(0, 4, 1)
#bs <- seq(0, 4, 0.5)
i <- 1
Z1s <- list()
powers <- list()
biass <- list()
for (b in bs){
    f.name <- paste0("./results/bs/", typ, "-b-", b, "-N-", N, 
    "-lam-",  lam, "-phi0-", phi0, "-invgam2-", invgam2, "-nSimu-", nSimu, ".RData")
    load(f.name)
    kpidxs <- sapply(post.res, function(res)length(res$res$phi0.tk)) != 100
    post.res <- post.res[kpidxs] 
    
    probMat <- sapply(post.res, post.b.fn, dlt0=dlt0)
    power <- mean(probMat[1, ] > eps0)
    power.no <- mean(probMat[2, ] > eps0.no)
    powers[[i]] <- c(power, power.no)
    
    
    trtMats <- sapply(post.res, function(x)x$mtrt)
    rowMeans(trtMats)
    c(var(trtMats[1, ]), var(trtMats[2, ]))
    errs <- abs(trtMats - b)
    biass[[i]] <- rbind(CI.fn(errs[1, ]), CI.fn(errs[2, ]))
    
    Z1.rates <- sapply(post.res, function(i)mean(i$data$Z))
    Z1.rates.no <- sapply(post.res, function(i)mean(i$data.no$Z))
    Z1s[[i]] <-  rbind(CI.fn(Z1.rates), CI.fn(Z1.rates.no))
    i <- i+1
}

powersMat <- do.call(rbind, powers)
biassMat.b <- sapply(biass, function(x)x[1, ])
biassMat.no <- sapply(biass, function(x)x[2, ])
Z1s.b <- sapply(Z1s, function(x)x[1, ])
Z1s.no <- sapply(Z1s, function(x)x[2, ])

res.b <- data.frame(b=bs, 
                    Z1rateLow=Z1s.b[1, ], 
                    Z1rateMean=Z1s.b[2, ], 
                    Z1rateUp=Z1s.b[3, ], 
                    power=powersMat[, 1], 
                    biasLow=biassMat.b[1, ],
                    biasMean=biassMat.b[2, ],
                    biasUp=biassMat.b[3, ])
res.no <- data.frame(b=bs, 
                    Z1rateLow=Z1s.no[1, ], 
                    Z1rateMean=Z1s.no[2, ], 
                    Z1rateUp=Z1s.no[3, ], 
                    power=powersMat[, 2], 
                    biasLow=biassMat.no[1, ],
                    biasMean=biassMat.no[2, ],
                    biasUp=biassMat.no[3, ])


# plot the results
png(filename = paste0("./plots/bias_", typ, "_phi0_", phi0, ".png"), width=360, height=480)
plot(bs, res.b$biasMean, type="l", lwd=2, col=1, xlab="b", ylab='Absolute bias', lty=1, ylim=c(0.0, 0.25))
#plot(bs, res.b$biasMean, type="l", lwd=2, col=1, xlab="b", ylab='Absolute bias', lty=1, ylim=c(0.6, 0.85))
polygon(c(bs, rev(bs)), c(res.b$biasLow, rev(res.b$biasUp)),  col= rgb(0, 0, 0,0.1), border=NA)
lines(bs, res.no$biasMean, lwd=2, col=2, lty=2)
#polygon(c(bs, rev(bs)), c(res.no$biasLow, rev(res.no$biasUp)),  col= rgb(1, 0, 0,0.1), border=NA)
legend("bottomleft", c("Borrow", "no Borrow"), col=c(1, 2), lty=c(1, 2))
dev.off()


png(filename = paste0("./plots/power_", typ, "_phi0_", phi0, ".png"), width=360, height=480)
plot(bs, res.b$power, type="l", lwd=2, col=1, xlab="b", ylab='Power', lty=1, ylim=c(0, 1))
lines(bs, res.no$power, lwd=2, col=2, lty=2)
legend("bottomright", c("Borrow", "no Borrow"), col=c(1, 2), lty=c(1, 2))
dev.off()

png(filename = paste0("./plots/Z1rate_", typ, "_phi0_", phi0, ".png"), width=360, height=480)
plot(bs, res.b$Z1rateMean, type="l", lwd=2, col=1, xlab="b", ylab='Percentage of treatment group', lty=1, ylim=c(0.45, 0.95))
polygon(c(bs, rev(bs)), c(res.b$Z1rateLow, rev(res.b$Z1rateUp)),  col= rgb(0, 0, 0,0.1), border=NA)
lines(bs, res.no$Z1rateMean, lwd=2, col=2, lty=2)
legend("topleft", c("Borrow", "no Borrow"), col=c(1, 2), lty=c(1, 2))
dev.off()
}

#3. vary N
{
typ <- "Simplex-Linear-diff-0"
lam <- 20
b <- 4
phi0 <- 2
invgam2 <- 1
nSimu  <- 1000
#allocation ratio vs Var(Y|X)/Var(Y)
Ns<- c(50, 75, 100, 125, 150, 175)
Ns<- seq(30, 200, 10)
Z1noMs <- c()
Z1Ms.CIs <- list()
var.rateMs <- c()
varMs <- c()
varXMs <- c()
biass <- list()
i <- 1
for (N in Ns){
    f.name <- paste0("./results/Ns/", typ, "-b-", b, "-N-", N, "-lam-",  lam, "-phi0-", phi0, ".RData")
    load(f.name)
    kpidxs <- sapply(post.res, function(res)length(res$res$phi0.tk)) != 100
    post.res <- post.res[kpidxs]
    varYs <- sapply(post.res, function(s.res)var(s.res$data$Y))
    varMs <- c(varMs, mean(varYs))
    Z1rates <- sapply(post.res, function(s.res)mean(s.res$data$Z))
    Z1rates.no <- sapply(post.res, function(s.res)mean(s.res$data.no$Z))
    varYXs <- phi0*2
    varXMs <- c(varXMs, mean(varYXs))
    var.rates <- varYXs/varYs
    Z1Ms.CIs[[i]] <- CI.fn(Z1rates)
    Z1noMs <- c(Z1noMs, mean(Z1rates.no))
    var.rateMs <- c(var.rateMs, mean(var.rates))
    
    trtMats <- sapply(post.res, function(x)x$mtrt)
    print(rowMeans(trtMats))
    errs <- abs(trtMats - b)
    biass[[i]] <- rbind(CI.fn(errs[1, ]), CI.fn(errs[2, ]))
    
    i <- i+1
    
}

# plot about allocation rate
Z1MsMat <- do.call(rbind, Z1Ms.CIs)
png(filename = paste0("./plots/Z1rateVsN_", typ, "_phi0_", phi0, ".png"), width=360, height=480)
plot(Ns, Z1MsMat[, 2], type="l", xlab="N", ylab="Percentage of treatment group", ylim=c(0.45, 0.75))
polygon(c(Ns, rev(Ns)), c(Z1MsMat[, 1], rev(Z1MsMat[, 3])),  col= rgb(0, 0, 0, 0.1), border=NA)
lines(Ns, Z1noMs, lwd=2, col=2, lty=2)
legend("topright", c("Borrow", "no Borrow"), col=c(1, 2), lty=c(1, 2))
dev.off()

png(filename = paste0("./plots/Z1rateVsRRN_", typ, "_phi0_", phi0, ".png"), width=360, height=480)
plot(var.rateMs, Z1MsMat[, 2], type="p", xlab="Var(Y|X, Z)/Var(Y)", ylab="Percentage of treatment group")
dev.off()

# plot about bias
png(filename = paste0("./plots/biasVsN_", typ, "_phi0_", phi0, ".png"), width=360, height=480)
biassMat.b <- sapply(biass, function(x)x[1, ])
biassMat.no <- sapply(biass, function(x)x[2, ])
plot(Ns, biassMat.b[2, ], type="l", lwd=2, col=1, xlab="N", ylab='Absolute bias', lty=1, ylim=c(0.2, 1.0))
polygon(c(Ns, rev(Ns)), c(biassMat.b[1, ], rev(biassMat.b[3, ])),  col= rgb(0, 0, 0, 0.1), border=NA)
lines(Ns, biassMat.no[2, ], lwd=2, col=2, lty=2)
legend("topleft", c("Borrow", "no Borrow"), col=c(1, 2), lty=c(1, 2))
dev.off()
}

#4. Subgroup analysis
{
    ext.sub.info.fn <- function(res){
        subId <- res$data$X1*2 + res$data$X2 + 1
        res.info <- list(Z=res$data$Z, subId=subId, tau2s=res$tau2s)
        Z1ms <- sapply(1:4, function(i) {mean(res.info$Z[res.info$subId==i])})
        tau2ms <-  sapply(1:4, function(i) {mean(res.info$tau2s[res.info$subId==i])})
        rv <- list(Z1ms=Z1ms, tau2ms=tau2ms)
        rv
    }    

    
typ <- "Simplex-Linear-diff-Sub10"
N <- 150
lam <- 20
b <- 4
invgam2 <- 1
nSimu  <- 1000
phi0 <- 1
f.name <- paste0("./results/SubGroups/", typ, "-b-", b, "-N-", N, 
    "-lam-",  lam, "-phi0-", phi0, "-invgam2-", invgam2, "-nSimu-", nSimu, ".RData")
load(f.name)

# remove the non-converging results
kpidxs <- sapply(post.res, function(res)length(res$res$phi0.tk)) != 100
post.res <- post.res[kpidxs]
mean(sapply(post.res, function(res)length(res$res$phi0.tk)) != 100)

png(filename = paste0("./plots/tausubgroup", typ, "_phi0_", phi0, ".png"), width=360, height=480)
tau2mss <- lapply(post.res, function(res)ext.sub.info.fn(res$res)$tau2ms)
tau2Mat <- do.call(rbind, tau2mss)
boxplot(tau2Mat, ylab="Tau2", xlab="Subgroup Index")
dev.off()

png(filename = paste0("./plots/Z1ratesubgroup", typ, "_phi0_", phi0, ".png"), width=360, height=480)
Z1mss <- lapply(post.res, function(res)ext.sub.info.fn(res$res)$Z1ms)
Z1Mat <- do.call(rbind, Z1mss)
boxplot(Z1Mat, ylab="Percentage of treatment group", xlab="Subgroup Index")
dev.off()

# trt effect bias
trtMats <- sapply(post.res, function(x)x$mtrt)
errs <- abs(trtMats - b)
trtRes <- rbind(CI.fn(errs[1, ]), CI.fn(errs[2, ]))
rownames(trtRes) <- c("Borrow", "Non-borrow")
trtRes
}


# 5, N, all discrete
{
typ <- "Simplex-Linear-AllBinary-diff-0"
lam <- 20
b <- 4
phi0 <- 2
invgam2 <- 1
nSimu  <- 1000
#allocation ratio vs Var(Y|X)/Var(Y)
Ns<- seq(30, 150, 30)
Z1noMs <- c()
Z1Ms.CIs <- list()
var.rateMs <- c()
varMs <- c()
varXMs <- c()
biass <- list()
i <- 1
for (N in Ns){
    f.name <- paste0("./results/", typ, "-b-", b, "-N-", N,  "-lam-",  lam, "-phi0-", phi0, "-invgam2-", invgam2, "-nSimu-", nSimu, ".RData")
    load(f.name)
    kpidxs <- sapply(post.res, function(res)length(res$res$phi0.tk)) != 100
    post.res <- post.res[kpidxs]
    varYs <- sapply(post.res, function(s.res)var(s.res$data$Y))
    varMs <- c(varMs, mean(varYs))
    Z1rates <- sapply(post.res, function(s.res)mean(s.res$data$Z))
    Z1rates.no <- sapply(post.res, function(s.res)mean(s.res$data.no$Z))
    varYXs <- phi0*2
    varXMs <- c(varXMs, mean(varYXs))
    var.rates <- varYXs/varYs
    Z1Ms.CIs[[i]] <- CI.fn(Z1rates)
    Z1noMs <- c(Z1noMs, mean(Z1rates.no))
    var.rateMs <- c(var.rateMs, mean(var.rates))
    
    trtMats <- sapply(post.res, function(x)x$mtrt)
    print(rowMeans(trtMats))
    errs <- abs(trtMats - b)
    biass[[i]] <- rbind(CI.fn(errs[1, ]), CI.fn(errs[2, ]))
    
    i <- i+1
    
}

# plot about allocation rate
Z1MsMat <- do.call(rbind, Z1Ms.CIs)
#png(filename = paste0("./plots/Z1rateVsN_", typ, "_phi0_", phi0, ".png"), width=360, height=480)
plot(Ns, Z1MsMat[, 2], type="l", xlab="N", ylab="Percentage of treatment group", ylim=c(0.45, 0.75))
polygon(c(Ns, rev(Ns)), c(Z1MsMat[, 1], rev(Z1MsMat[, 3])),  col= rgb(0, 0, 0, 0.1), border=NA)
lines(Ns, Z1noMs, lwd=2, col=2, lty=2)
legend("topright", c("Borrow", "no Borrow"), col=c(1, 2), lty=c(1, 2))
#dev.off()

png(filename = paste0("./plots/Z1rateVsRRN_", typ, "_phi0_", phi0, ".png"), width=360, height=480)
plot(var.rateMs, Z1MsMat[, 2], type="p", xlab="Var(Y|X, Z)/Var(Y)", ylab="Percentage of treatment group")
dev.off()

# plot about bias
png(filename = paste0("./plots/biasVsN_", typ, "_phi0_", phi0, ".png"), width=360, height=480)
biassMat.b <- sapply(biass, function(x)x[1, ])
biassMat.no <- sapply(biass, function(x)x[2, ])
plot(Ns, biassMat.b[2, ], type="l", lwd=2, col=1, xlab="N", ylab='Absolute bias', lty=1, ylim=c(0.2, 1.0))
polygon(c(Ns, rev(Ns)), c(biassMat.b[1, ], rev(biassMat.b[3, ])),  col= rgb(0, 0, 0, 0.1), border=NA)
lines(Ns, biassMat.no[2, ], lwd=2, col=2, lty=2)
legend("topleft", c("Borrow", "no Borrow"), col=c(1, 2), lty=c(1, 2))
dev.off()
}


#6. Varying b, all discrete
{
    # for hypothesis test
post.b.fn2 <- function(res, dlt0=1){
    sps.trts <- res$sps.trts
    trtsMat <- sps.trts
    #trtsMat <- abs(sps.trts)
    prb <- post.fn(trtsMat[, 1], dlt0)
    prb.no <- post.fn(trtsMat[, 2], dlt0)
    c(prb, prb.no)
}
    
    dlt0 <- 0.5
    
    typ <- "Linear"
    diagH <- rep(0.05, 4)
    x.tps <- c(2, 2, 2, 2)
    hs <- rep(1.1, 4)
    xis <- c(0, 0, 0, 0)
    idx.tps <- paste0(x.tps, collapse = "")
    N <- 150
    lam <- 1e6
    b <- 0
    phi0 <- 2
    invgam2 <- 1
    nSimu <- 1500
    f.name <- paste0("./results/", typ, "-b-", b, "-N-", N, "-lam-", lam, "-phi0-", phi0, "-invgam2-", invgam2, 
                     "-H-", vec2code(diagH, 100), "-h-", vec2code(hs, 100), "-tps-", idx.tps, "-xis-", vec2code(xis, 10), "-nSimu-", nSimu, ".RData")
    #f.name <- paste0("./results/", typ, "-b-", b, "-N-", N, 
    #                 "-lam-",  lam, "-phi0-", phi0, "-invgam2-", invgam2, "-nSimu-", nSimu, ".RData")
    load(f.name)
    kpidxs <- sapply(post.res, function(res)length(res$res$phi0.tk)) != 100
    post.res <- post.res[kpidxs];length(post.res)
    
    # eps0 and eps0.no are tuned such that the size is 0.05 
    probMat <- sapply(post.res, post.b.fn2, dlt0=dlt0)
    eps.alig(probMat[1, ])
    eps.alig(probMat[2, ])
    eps0 <- eps.alig(probMat[1, ])$eps0
    eps0.no <- eps.alig(probMat[2, ])$eps0
    mean(probMat[1, ]>eps0)
    mean(probMat[2, ]>eps0.no)
    #eps0 <- eps0.no
    
    bs <- c(0, 0.5, 1, 2, 3, 4)
    i <- 1
    Z1s <- list()
    powers <- list()
    biass <- list()
    for (b in bs){
        f.name <- paste0("./results/", typ, "-b-", b, "-N-", N, "-lam-", lam, "-phi0-", phi0, "-invgam2-", invgam2, 
                     "-H-", vec2code(diagH, 100), "-h-", vec2code(hs, 100), "-tps-", idx.tps, "-xis-", vec2code(xis, 10), "-nSimu-", nSimu, ".RData")
        #f.name <- paste0("./results/", typ, "-b-", b, "-N-", N, 
        #                 "-lam-",  lam, "-phi0-", phi0, "-invgam2-", invgam2, "-nSimu-", nSimu, ".RData")
        load(f.name)
        kpidxs <- sapply(post.res, function(res)length(res$res$phi0.tk)) != 100
        post.res <- post.res[kpidxs] 
        
        probMat <- sapply(post.res, post.b.fn2, dlt0=dlt0)
        power <- mean(probMat[1, ] > eps0)
        power.no <- mean(probMat[2, ] > eps0.no)
        powers[[i]] <- c(power, power.no)
        
        
        trtMats <- sapply(post.res, function(x)x$mtrt)
        print(rowMeans(trtMats))
        c(var(trtMats[1, ]), var(trtMats[2, ]))
        errs <- abs(trtMats - b)
        biass[[i]] <- rbind(CI.fn(errs[1, ]), CI.fn(errs[2, ]))
        
        Z1.rates <- sapply(post.res, function(i)mean(i$data$Z))
        Z1.rates.no <- sapply(post.res, function(i)mean(i$data.no$Z))
        Z1s[[i]] <-  rbind(CI.fn(Z1.rates), CI.fn(Z1.rates.no))
        i <- i+1
    }
    
    powersMat <- do.call(rbind, powers)
    biassMat.b <- sapply(biass, function(x)x[1, ])
    biassMat.no <- sapply(biass, function(x)x[2, ])
    Z1s.b <- sapply(Z1s, function(x)x[1, ])
    Z1s.no <- sapply(Z1s, function(x)x[2, ])
    
    res.b <- data.frame(b=bs, 
                        Z1rateLow=Z1s.b[1, ], 
                        Z1rateMean=Z1s.b[2, ], 
                        Z1rateUp=Z1s.b[3, ], 
                        power=powersMat[, 1], 
                        biasLow=biassMat.b[1, ],
                        biasMean=biassMat.b[2, ],
                        biasUp=biassMat.b[3, ])
    res.no <- data.frame(b=bs, 
                         Z1rateLow=Z1s.no[1, ], 
                         Z1rateMean=Z1s.no[2, ], 
                         Z1rateUp=Z1s.no[3, ], 
                         power=powersMat[, 2], 
                         biasLow=biassMat.no[1, ],
                         biasMean=biassMat.no[2, ],
                         biasUp=biassMat.no[3, ])
    
    
    # plot the results
    #png(filename = paste0("./plots/bias_", typ, "_phi0_", phi0, ".png"), width=360, height=480)
    plot(bs, res.b$biasMean, type="l", lwd=2, col=1, xlab="b", ylab='Absolute bias', lty=1, ylim=c(0.1, 0.3))
    #plot(bs, res.b$biasMean, type="l", lwd=2, col=1, xlab="b", ylab='Absolute bias', lty=1, ylim=c(0.6, 0.85))
    polygon(c(bs, rev(bs)), c(res.b$biasLow, rev(res.b$biasUp)),  col= rgb(0, 0, 0,0.1), border=NA)
    lines(bs, res.no$biasMean, lwd=2, col=2, lty=2)
    #polygon(c(bs, rev(bs)), c(res.no$biasLow, rev(res.no$biasUp)),  col= rgb(1, 0, 0,0.1), border=NA)
    legend("bottomleft", c("Borrow", "no Borrow"), col=c(1, 2), lty=c(1, 2))
    #dev.off()
    
    
    #png(filename = paste0("./plots/power_", typ, "_phi0_", phi0, ".png"), width=360, height=480)
    plot(bs, res.b$power, type="l", lwd=2, col=1, xlab="b", ylab='Power', lty=1, ylim=c(0, 1))
    lines(bs, res.no$power, lwd=2, col=2, lty=2)
    legend("bottomright", c("Borrow", "no Borrow"), col=c(1, 2), lty=c(1, 2))
    #dev.off()
    
    #png(filename = paste0("./plots/Z1rate_", typ, "_phi0_", phi0, ".png"), width=360, height=480)
    plot(bs, res.b$Z1rateMean, type="l", lwd=2, col=1, xlab="b", ylab='Percentage of treatment group', lty=1, ylim=c(0.45, 0.95))
    polygon(c(bs, rev(bs)), c(res.b$Z1rateLow, rev(res.b$Z1rateUp)),  col= rgb(0, 0, 0,0.1), border=NA)
    lines(bs, res.no$Z1rateMean, lwd=2, col=2, lty=2)
    legend("topleft", c("Borrow", "no Borrow"), col=c(1, 2), lty=c(1, 2))
    #dev.off()
}

#7. Subgroup Varying b, all discrete
{
    {
    # for hypothesis test
post.b.fn2 <- function(res, dlt0=1){
    sps.trts <- res$sps.trts
    trtsMat <- sps.trts
    #trtsMat <- abs(sps.trts)
    prb <- post.fn(trtsMat[, 1], dlt0)
    prb.no <- post.fn(trtsMat[, 2], dlt0)
    c(prb, prb.no)
}
    
    
    dlt0 <- 0.5
    
    typ <- "Linear"
    diagH <- rep(0.05, 4)
    x.tps <- c(2, 2, 2, 2)
    hs <- rep(1.1, 4)
    xis <- c(0, 2, 2, 0)
    idx.tps <- paste0(x.tps, collapse = "")
    N <- 150
    lam <- 1e6
    b <- 0
    phi0 <- 2
    invgam2 <- 1
    nSimu <- 1500
    f.name <- paste0("./results/", typ, "-b-", b, "-N-", N, "-lam-", lam, "-phi0-", phi0, "-invgam2-", invgam2, 
                     "-H-", vec2code(diagH, 100), "-h-", vec2code(hs, 100), "-tps-", idx.tps, "-xis-", vec2code(xis, 10), "-nSimu-", nSimu, ".RData")
    #f.name <- paste0("./results/", typ, "-b-", b, "-N-", N, 
    #                 "-lam-",  lam, "-phi0-", phi0, "-invgam2-", invgam2, "-nSimu-", nSimu, ".RData")
    load(f.name)
    kpidxs <- sapply(post.res, function(res)length(res$res$phi0.tk)) != 100
    post.res <- post.res[kpidxs];length(post.res)
    
    # eps0 and eps0.no are tuned such that the size is 0.05 
    probMat <- sapply(post.res, post.b.fn2, dlt0=dlt0)
    eps.alig(probMat[1, ])
    eps.alig(probMat[2, ])
    eps0 <- eps.alig(probMat[1, ])$eps0
    eps0.no <- eps.alig(probMat[2, ])$eps0
    mean(probMat[1, ]>eps0)
    mean(probMat[2, ]>eps0.no)
    #eps0 <- eps0.no
    
    bs <- c(0, 0.5, 1, 2, 3, 4)
    i <- 1
    Z1s <- list()
    powers <- list()
    biass <- list()
    for (b in bs){
        f.name <- paste0("./results/", typ, "-b-", b, "-N-", N, "-lam-", lam, "-phi0-", phi0, "-invgam2-", invgam2, 
                     "-H-", vec2code(diagH, 100), "-h-", vec2code(hs, 100), "-tps-", idx.tps, "-xis-", vec2code(xis, 10), "-nSimu-", nSimu, ".RData")
        #f.name <- paste0("./results/", typ, "-b-", b, "-N-", N, 
        #                 "-lam-",  lam, "-phi0-", phi0, "-invgam2-", invgam2, "-nSimu-", nSimu, ".RData")
        load(f.name)
        kpidxs <- sapply(post.res, function(res)length(res$res$phi0.tk)) != 100
        post.res <- post.res[kpidxs] 
        
        probMat <- sapply(post.res, post.b.fn2, dlt0=dlt0)
        power <- mean(probMat[1, ] > eps0)
        power.no <- mean(probMat[2, ] > eps0.no)
        powers[[i]] <- c(power, power.no)
        
        
        trtMats <- sapply(post.res, function(x)x$mtrt)
        print(rowMeans(trtMats))
        c(var(trtMats[1, ]), var(trtMats[2, ]))
        errs <- abs(trtMats - b)
        biass[[i]] <- rbind(CI.fn(errs[1, ]), CI.fn(errs[2, ]))
        
        Z1.rates <- sapply(post.res, function(i)mean(i$data$Z))
        Z1.rates.no <- sapply(post.res, function(i)mean(i$data.no$Z))
        Z1s[[i]] <-  rbind(CI.fn(Z1.rates), CI.fn(Z1.rates.no))
        i <- i+1
    }
    
    powersMat <- do.call(rbind, powers)
    biassMat.b <- sapply(biass, function(x)x[1, ])
    biassMat.no <- sapply(biass, function(x)x[2, ])
    Z1s.b <- sapply(Z1s, function(x)x[1, ])
    Z1s.no <- sapply(Z1s, function(x)x[2, ])
    
    res.b <- data.frame(b=bs, 
                        Z1rateLow=Z1s.b[1, ], 
                        Z1rateMean=Z1s.b[2, ], 
                        Z1rateUp=Z1s.b[3, ], 
                        power=powersMat[, 1], 
                        biasLow=biassMat.b[1, ],
                        biasMean=biassMat.b[2, ],
                        biasUp=biassMat.b[3, ])
    res.no <- data.frame(b=bs, 
                         Z1rateLow=Z1s.no[1, ], 
                         Z1rateMean=Z1s.no[2, ], 
                         Z1rateUp=Z1s.no[3, ], 
                         power=powersMat[, 2], 
                         biasLow=biassMat.no[1, ],
                         biasMean=biassMat.no[2, ],
                         biasUp=biassMat.no[3, ])
    
    
    # plot the results
    #png(filename = paste0("./plots/bias_", typ, "_phi0_", phi0, ".png"), width=360, height=480)
    plot(bs, res.b$biasMean, type="l", lwd=2, col=1, xlab="b", ylab='Absolute bias', lty=1, ylim=c(0.2, 0.7))
    #plot(bs, res.b$biasMean, type="l", lwd=2, col=1, xlab="b", ylab='Absolute bias', lty=1, ylim=c(0.6, 0.85))
    polygon(c(bs, rev(bs)), c(res.b$biasLow, rev(res.b$biasUp)),  col= rgb(0, 0, 0,0.1), border=NA)
    lines(bs, res.no$biasMean, lwd=2, col=2, lty=2)
    #polygon(c(bs, rev(bs)), c(res.no$biasLow, rev(res.no$biasUp)),  col= rgb(1, 0, 0,0.1), border=NA)
    legend("bottomleft", c("Borrow", "no Borrow"), col=c(1, 2), lty=c(1, 2))
    #dev.off()
    
    
    #png(filename = paste0("./plots/power_", typ, "_phi0_", phi0, ".png"), width=360, height=480)
    plot(bs, res.b$power, type="l", lwd=2, col=1, xlab="b", ylab='Power', lty=1, ylim=c(0, 1))
    lines(bs, res.no$power, lwd=2, col=2, lty=2)
    legend("bottomright", c("Borrow", "no Borrow"), col=c(1, 2), lty=c(1, 2))
    #dev.off()
    
    #png(filename = paste0("./plots/Z1rate_", typ, "_phi0_", phi0, ".png"), width=360, height=480)
    plot(bs, res.b$Z1rateMean, type="l", lwd=2, col=1, xlab="b", ylab='Percentage of treatment group', lty=1, ylim=c(0.45, 0.95))
    polygon(c(bs, rev(bs)), c(res.b$Z1rateLow, rev(res.b$Z1rateUp)),  col= rgb(0, 0, 0,0.1), border=NA)
    lines(bs, res.no$Z1rateMean, lwd=2, col=2, lty=2)
    legend("topleft", c("Borrow", "no Borrow"), col=c(1, 2), lty=c(1, 2))
    #dev.off()
    }
    
    {  
    ext.sub.info.fn <- function(res){
        subId <- res$data$X1*2 + res$data$X2 + 1
        res.info <- list(Z=res$data$Z, subId=subId, tau2s=res$tau2s)
        Z1ms <- sapply(1:4, function(i) {mean(res.info$Z[res.info$subId==i])})
        tau2ms <-  sapply(1:4, function(i) {mean(res.info$tau2s[res.info$subId==i])})
        rv <- list(Z1ms=Z1ms, tau2ms=tau2ms)
        rv
    }    
    
    b <- 2
    f.name <- paste0("./results/", typ, "-b-", b, "-N-", N, "-lam-", lam, "-phi0-", phi0, "-invgam2-", invgam2, 
                     "-H-", vec2code(diagH, 100), "-h-", vec2code(hs, 100), "-tps-", idx.tps, "-xis-", vec2code(xis, 10), "-nSimu-", nSimu, ".RData")
    #f.name <- paste0("./results/", typ, "-b-", b, "-N-", N, 
    #                 "-lam-",  lam, "-phi0-", phi0, "-invgam2-", invgam2, "-nSimu-", nSimu, ".RData")
    load(f.name)
    kpidxs <- sapply(post.res, function(res)length(res$res$phi0.tk)) != 100
    post.res <- post.res[kpidxs];length(post.res)
    
    
    #png(filename = paste0("./plots/tausubgroup", typ, "_phi0_", phi0, ".png"), width=360, height=480)
    tau2mss <- lapply(post.res, function(res)ext.sub.info.fn(res$res)$tau2ms)
    tau2Mat <- do.call(rbind, tau2mss)
    colMeans(tau2Mat)
    boxplot(tau2Mat, ylab="Tau2", xlab="Subgroup Index")
    #dev.off()
    
    #png(filename = paste0("./plots/Z1ratesubgroup", typ, "_phi0_", phi0, ".png"), width=360, height=480)
    Z1mss <- lapply(post.res, function(res)ext.sub.info.fn(res$res)$Z1ms)
    Z1Mat <- do.call(rbind, Z1mss)
    colMeans(Z1Mat)
    boxplot(Z1Mat, ylab="Percentage of treatment group", xlab="Subgroup Index")
    #dev.off()
    
    # trt effect bias
    trtMats <- sapply(post.res, function(x)x$mtrt)
    errs <- abs(trtMats - b)
    trtRes <- rbind(CI.fn(errs[1, ]), CI.fn(errs[2, ]))
    rownames(trtRes) <- c("Borrow", "Non-borrow")
    trtRes
    }    
    
}