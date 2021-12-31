setwd("C:/Users/JINHU/Documents/ProjectCode/HistTrial")
library(magrittr)
source("util_ana.R")




#1. vary phi0
{

typ <- "LinearSw"
diagH <- rep(0.05, 4)
x.tps <- c(2, 2, 2, 2)
hs <- rep(1.1, 4)
xis <- c(0, 0, 0, 0) 
idx.tps <- paste0(x.tps, collapse = "")
lam.q <- 0.10
lam <- 30
b <- 4
phi0 <- 5
invgam2 <- 10
nSimu  <- 1000
N <- 150
phi0s <- c(1, 2, 3, 4, 5, 6, 8, 10)
Z1Ms.CIs <- list()
Z1noMs <- c()
var.rateMs <- c()
varMs <- c()
varXMs <- c()
biass <- list()
i <- 1
for (phi0 in phi0s){
    f.name <- paste0("./results/", typ, "-b-", b, "-N-", N, "-lam-", lam, "-lamq-", 100*lam.q, "-phi0-", phi0, "-invgam2-", invgam2, 
                     "-H-", vec2code(diagH, 100), "-h-", vec2code(hs, 100), "-tps-", idx.tps, "-xis-", vec2code(xis, 10), "-nSimu-", nSimu, ".RData")
    load(f.name)
    kpidxs <- sapply(post.res, function(res)length(res$res$phi0.tk)) != 100
    post.res <- post.res[kpidxs]
    varYs <- sapply(post.res, function(s.res)var(s.res$data$Y))
    varMs <- c(varMs, mean(varYs))
    Z1rates <- sapply(post.res, function(s.res)mean(s.res$data$Z[21:N]))
    Z1rates.no <- sapply(post.res, function(s.res)mean(s.res$data.no$Z[21:N]))
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
png(filename = paste0("./plots/Z1rateVsphi0_", typ, ".png"), width=360, height=480)
plot(phi0s, Z1MsMat[, 2], type="l", xlab="Var(Y|X, Z)", ylab="Percentage of treatment group", ylim=c(0.45, 0.95))
polygon(c(phi0s, rev(phi0s)), c(Z1MsMat[, 1], rev(Z1MsMat[, 3])),  col= rgb(0, 0, 0, 0.1), border=NA)
lines(phi0s, Z1noMs, lwd=2, col=2, lty=2)
legend("topleft", c("Borrow", "no Borrow"), col=c(1, 2), lty=c(1, 2))
dev.off()

Z1MsMat[, 2] %>% round(2)

#png(filename = paste0("./plots/Z1rateVsrv", typ, ".png"), width=360, height=480)
plot(var.rateMs, Z1MsMat[, 2], type="p", xlab="Var(Y|X, Z)/Var(Y)", ylab="Percentage of treatment group")
#dev.off()

# plot about bias
png(filename = paste0("./plots/biasVsphi0_", typ, ".png"), width=360, height=480)
biassMat.b <- sapply(biass, function(x)x[1, ])
biassMat.no <- sapply(biass, function(x)x[2, ])
plot(phi0s, biassMat.b[2, ], type="l", lwd=2, col=1, xlab="Var(Y|X, Z)", ylab='Absolute bias', lty=1, ylim=c(0.15, 1.3))
polygon(c(phi0s, rev(phi0s)), c(biassMat.b[1, ], rev(biassMat.b[3, ])),  col= rgb(0, 0, 0, 0.1), border=NA)
lines(phi0s, biassMat.no[2, ], lwd=2, col=2, lty=2)
#polygon(c(phi0s, rev(phi0s)), c(biassMat.no[1, ], rev(biassMat.no[3, ])),  col= rgb(1, 0, 0,0.1), border=NA)
legend("topleft", c("Borrow", "no Borrow"), col=c(1, 2), lty=c(1, 2))
dev.off()

rbind(biassMat.b[3, ], biassMat.no[2, ], Ns)
diff <- biassMat.no[2, ] - biassMat.b[2, ]
plot(phi0s, diff, type="b")
}

#2. vary N, all discrete
{
typ <- "LinearSw"
diagH <- rep(0.05, 4)
x.tps <- c(2, 2, 2, 2)
hs <- rep(1.1, 4)
xis <- c(0, 0, 0, 0) 
idx.tps <- paste0(x.tps, collapse = "")
lam.q <- 0.10
lam <- 30
b <- 4
phi0 <- 5
invgam2 <- 10
nSimu  <- 1000
#allocation ratio vs Var(Y|X)/Var(Y)
Ns<- c(seq(30, 100, 10), seq(125, 200, 25), seq(250, 400, 50))
Z1noMs <- c()
Z1Ms.CIs <- list()
var.rateMs <- c()
varMs <- c()
varXMs <- c()
biass <- list()
i <- 1
for (N in Ns){
    f.name <- paste0("./results/", typ, "-b-", b, "-N-", N, "-lam-", lam, "-lamq-", 100*lam.q, "-phi0-", phi0, "-invgam2-", invgam2, 
                     "-H-", vec2code(diagH, 100), "-h-", vec2code(hs, 100), "-tps-", idx.tps, "-xis-", vec2code(xis, 10), "-nSimu-", nSimu, ".RData")
    load(f.name)
    kpidxs <- sapply(post.res, function(res)length(res$res$phi0.tk)) != 100
    post.res <- post.res[kpidxs];print(length(post.res))
    varYs <- sapply(post.res, function(s.res)var(s.res$data$Y))
    varMs <- c(varMs, mean(varYs))
    Z1rates <- sapply(post.res, function(s.res)mean(s.res$data$Z[21:N]))
    Z1rates.no <- sapply(post.res, function(s.res)mean(s.res$data.no$Z[21:N]))
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
plot(Ns, Z1MsMat[, 2], type="l", xlab="N", ylab="Percentage of treatment group", ylim=c(0.45, 0.95))
polygon(c(Ns, rev(Ns)), c(Z1MsMat[, 1], rev(Z1MsMat[, 3])),  col= rgb(0, 0, 0, 0.1), border=NA)
lines(Ns, Z1noMs, lwd=2, col=2, lty=2)
legend("topright", c("Borrow", "no Borrow"), col=c(1, 2), lty=c(1, 2))
dev.off()

#png(filename = paste0("./plots/Z1rateVsRRN_", typ, "_phi0_", phi0, ".png"), width=360, height=480)
plot(var.rateMs, Z1MsMat[, 2], type="p", xlab="Var(Y|X, Z)/Var(Y)", ylab="Percentage of treatment group")
#dev.off()

# plot about bias
png(filename = paste0("./plots/biasVsN_", typ, "_phi0_", phi0, ".png"), width=360, height=480)
biassMat.b <- sapply(biass, function(x)x[1, ])
biassMat.no <- sapply(biass, function(x)x[2, ])
plot(Ns, biassMat.b[2, ], type="l", lwd=2, col=1, xlab="N", ylab='Absolute bias', lty=1, ylim=c(0.3, 1.5))
polygon(c(Ns, rev(Ns)), c(biassMat.b[1, ], rev(biassMat.b[3, ])),  col= rgb(0, 0, 0, 0.1), border=NA)
lines(Ns, biassMat.no[2, ], lwd=2, col=2, lty=2)
legend("topright", c("Borrow", "no Borrow"), col=c(1, 2), lty=c(1, 2))
dev.off()

rbind(biassMat.b[3, ], biassMat.no[2, ], Ns)
diff <- biassMat.no[2, ] - biassMat.b[2, ]
plot(Ns, diff, type="b")
}


#3. Varying b, all discrete
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
    
    typ <- "LinearSw"
    diagH <- rep(0.05, 4)
    x.tps <- c(2, 2, 2, 2)
    hs <- rep(1.1, 4)
    xis <- c(0, 0, 0, 0)
    idx.tps <- paste0(x.tps, collapse = "")
    N <- 150
    lam.q <- 0.10
    lam <- 30
    b <- 0
    phi0 <- 3
    invgam2 <- 10
    nSimu <- 1500
    #f.name <- paste0("./results/", typ, "-b-", b, "-N-", N, "-lamq-", 100*lam.q, "-phi0-", phi0, "-invgam2-", invgam2, 
    #                 "-H-", vec2code(diagH, 100), "-h-", vec2code(hs, 100), "-tps-", idx.tps, "-xis-", vec2code(xis, 10), "-nSimu-", nSimu, ".RData")
    f.name <- paste0("./results/", typ, "-b-", b, "-N-", N, "-lam-", lam, "-lamq-", 100*lam.q, "-phi0-", phi0, "-invgam2-", invgam2, 
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
    eps0.no <- 0.999
    mean(probMat[1, ]>eps0)
    mean(probMat[2, ]>eps0.no)
    #eps0 <- eps0.no
    
    bs <- c(0, 0.5, 1:3)
    i <- 1
    Z1s <- list()
    powers <- list()
    biass <- list()
    for (b in bs){
        f.name <- paste0("./results/", typ, "-b-", b, "-N-", N, "-lam-", lam, "-lamq-", 100*lam.q, "-phi0-", phi0, "-invgam2-", invgam2, 
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
        
        Z1.rates <- sapply(post.res, function(i)mean(i$data$Z[21:N]))
        Z1.rates.no <- sapply(post.res, function(i)mean(i$data.no$Z[21:N]))
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
    png(filename = paste0("./plots/inConsistentbias_", typ, "_phi0_", phi0, ".png"), width=360, height=480)
    plot(bs, res.b$biasMean, type="l", lwd=2, col=1, xlab="b", ylab='Absolute bias', lty=1, ylim=c(0.2, 0.4))
    #plot(bs, res.b$biasMean, type="l", lwd=2, col=1, xlab="b", ylab='Absolute bias', lty=1, ylim=c(0.6, 0.85))
    polygon(c(bs, rev(bs)), c(res.b$biasLow, rev(res.b$biasUp)),  col= rgb(0, 0, 0,0.1), border=NA)
    lines(bs, res.no$biasMean, lwd=2, col=2, lty=2)
    #polygon(c(bs, rev(bs)), c(res.no$biasLow, rev(res.no$biasUp)),  col= rgb(1, 0, 0,0.1), border=NA)
    legend("bottomleft", c("Borrow", "no Borrow"), col=c(1, 2), lty=c(1, 2))
    dev.off()
    
    
    png(filename = paste0("./plots/inConsistentpower_", typ, "_phi0_", phi0, ".png"), width=360, height=480)
    plot(bs, res.b$power, type="l", lwd=2, col=1, xlab="b", ylab='Power', lty=1, ylim=c(0, 1))
    lines(bs, res.no$power, lwd=2, col=2, lty=2)
    legend("bottomright", c("Borrow", "no Borrow"), col=c(1, 2), lty=c(1, 2))
    dev.off()
    
    png(filename = paste0("./plots/inConsistentZ1rate_", typ, "_phi0_", phi0, ".png"), width=360, height=480)
    plot(bs, res.b$Z1rateMean, type="l", lwd=2, col=1, xlab="b", ylab='Percentage of treatment group', lty=1, ylim=c(0.45, 0.95))
    polygon(c(bs, rev(bs)), c(res.b$Z1rateLow, rev(res.b$Z1rateUp)),  col= rgb(0, 0, 0,0.1), border=NA)
    lines(bs, res.no$Z1rateMean, lwd=2, col=2, lty=2)
    legend("topleft", c("Borrow", "no Borrow"), col=c(1, 2), lty=c(1, 2))
    dev.off()
}

#4. Subgroup Varying b, all discrete
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
    
    typ <- "LinearSw"
    diagH <- rep(0.05, 4)
    x.tps <- c(2, 2, 2, 2)
    hs <- rep(1.1, 4)
    xis <- c(0, 3, 3, 0)
    idx.tps <- paste0(x.tps, collapse = "")
    N <- 150
    lam.q <- 0.10
    lam <- 30
    b <- 0
    phi0 <- 3
    invgam2 <- 10
    nSimu <- 1500
    f.name <- paste0("./results/", typ, "-b-", b, "-N-", N, "-lam-", lam, "-lamq-", 100*lam.q, "-phi0-", phi0, "-invgam2-", invgam2, 
                     "-H-", vec2code(diagH, 100), "-h-", vec2code(hs, 100), "-tps-", idx.tps, "-xis-", vec2code(xis, 10), "-nSimu-", nSimu, ".RData")
    #f.name <- paste0("./results/", typ, "-b-", b, "-N-", N, "-lam-", lam, "-phi0-", phi0, "-invgam2-", invgam2, 
    #                 "-H-", vec2code(diagH, 100), "-h-", vec2code(hs, 100), "-tps-", idx.tps, "-xis-", vec2code(xis, 10), "-nSimu-", nSimu, ".RData")
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
    
    bs <- c(0, 0.5, 1:3)
    i <- 1
    Z1s <- list()
    powers <- list()
    biass <- list()
    for (b in bs){
        f.name <- paste0("./results/", typ, "-b-", b, "-N-", N, "-lam-", lam, "-lamq-", 100*lam.q, "-phi0-", phi0, "-invgam2-", invgam2, 
                     "-H-", vec2code(diagH, 100), "-h-", vec2code(hs, 100), "-tps-", idx.tps, "-xis-", vec2code(xis, 10), "-nSimu-", nSimu, ".RData")
        #f.name <- paste0("./results/", typ, "-b-", b, "-N-", N, "-lam-", lam, "-phi0-", phi0, "-invgam2-", invgam2, 
        #             "-H-", vec2code(diagH, 100), "-h-", vec2code(hs, 100), "-tps-", idx.tps, "-xis-", vec2code(xis, 10), "-nSimu-", nSimu, ".RData")
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
        
        Z1.rates <- sapply(post.res, function(i)mean(i$data$Z[21:N]))
        Z1.rates.no <- sapply(post.res, function(i)mean(i$data.no$Z[21:N]))
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
    png(filename = paste0("./plots/Subgroupbias_", typ, "_phi0_", phi0, ".png"), width=360, height=480)
    plot(bs, res.b$biasMean, type="l", lwd=2, col=1, xlab="b", ylab='Absolute bias', lty=1, ylim=c(0.2, 0.3))
    #plot(bs, res.b$biasMean, type="l", lwd=2, col=1, xlab="b", ylab='Absolute bias', lty=1, ylim=c(0.6, 0.85))
    polygon(c(bs, rev(bs)), c(res.b$biasLow, rev(res.b$biasUp)),  col= rgb(0, 0, 0,0.1), border=NA)
    lines(bs, res.no$biasMean, lwd=2, col=2, lty=2)
    #polygon(c(bs, rev(bs)), c(res.no$biasLow, rev(res.no$biasUp)),  col= rgb(1, 0, 0,0.1), border=NA)
    legend("bottomleft", c("Borrow", "no Borrow"), col=c(1, 2), lty=c(1, 2))
    dev.off()
    
    
    png(filename = paste0("./plots/Subgrouppower_", typ, "_phi0_", phi0, ".png"), width=360, height=480)
    plot(bs, res.b$power, type="l", lwd=2, col=1, xlab="b", ylab='Power', lty=1, ylim=c(0, 1))
    lines(bs, res.no$power, lwd=2, col=2, lty=2)
    legend("bottomright", c("Borrow", "no Borrow"), col=c(1, 2), lty=c(1, 2))
    dev.off()
    
    png(filename = paste0("./plots/SubgroupZ1rate_", typ, "_phi0_", phi0, ".png"), width=360, height=480)
    plot(bs, res.b$Z1rateMean, type="l", lwd=2, col=1, xlab="b", ylab='Percentage of treatment group', lty=1, ylim=c(0.45, 0.95))
    polygon(c(bs, rev(bs)), c(res.b$Z1rateLow, rev(res.b$Z1rateUp)),  col= rgb(0, 0, 0,0.1), border=NA)
    lines(bs, res.no$Z1rateMean, lwd=2, col=2, lty=2)
    legend("topleft", c("Borrow", "no Borrow"), col=c(1, 2), lty=c(1, 2))
    dev.off()
    }
    
    {  
    ext.sub.info.fn <- function(res){
        subId <- res$data$X1*2 + res$data$X2 + 1
        res.info <- list(Z=res$data$Z, subId=subId, tau2s=res$tau2s)
        Z1ms <- sapply(1:4, function(i) {mean(res.info$Z[21:N][res.info$subId[21:N]==i])})
        tau2ms <-  sapply(1:4, function(i) {mean(res.info$tau2s[res.info$subId==i])})
        rv <- list(Z1ms=Z1ms, tau2ms=tau2ms)
        rv
    }    
    
    b <- 2
    f.name <- paste0("./results/", typ, "-b-", b, "-N-", N, "-lam-", lam, "-lamq-", 100*lam.q, "-phi0-", phi0, "-invgam2-", invgam2, 
                     "-H-", vec2code(diagH, 100), "-h-", vec2code(hs, 100), "-tps-", idx.tps, "-xis-", vec2code(xis, 10), "-nSimu-", nSimu, ".RData")
    #f.name <- paste0("./results/", typ, "-b-", b, "-N-", N, "-lam-", lam, "-phi0-", phi0, "-invgam2-", invgam2, 
    #                 "-H-", vec2code(diagH, 100), "-h-", vec2code(hs, 100), "-tps-", idx.tps, "-xis-", vec2code(xis, 10), "-nSimu-", nSimu, ".RData")
    #f.name <- paste0("./results/", typ, "-b-", b, "-N-", N, 
    #                 "-lam-",  lam, "-phi0-", phi0, "-invgam2-", invgam2, "-nSimu-", nSimu, ".RData")
    load(f.name)
    kpidxs <- sapply(post.res, function(res)length(res$res$phi0.tk)) != 100
    post.res <- post.res[kpidxs];length(post.res)
    
    
    png(filename = paste0("./plots/SubgroupTau2Boxplot_", typ, "_phi0_", phi0, ".png"), width=360, height=480)
    tau2mss <- lapply(post.res, function(res)ext.sub.info.fn(res$res)$tau2ms)
    tau2Mat <- do.call(rbind, tau2mss)
    colMeans(tau2Mat)
    median(tau2Mat[, 1])
    median(tau2Mat[, 2])
    boxplot(tau2Mat, ylab="Tau2", xlab="Subgroup Index")
    dev.off()
    
    png(filename = paste0("./plots/SubgroupZ1rateBoxplot_", typ, "_phi0_", phi0, ".png"), width=360, height=480)
    Z1mss <- lapply(post.res, function(res)ext.sub.info.fn(res$res)$Z1ms)
    Z1Mat <- do.call(rbind, Z1mss)
    colMeans(Z1Mat)
    boxplot(Z1Mat, ylab="Percentage of treatment group", xlab="Subgroup Index")
    abline(h=0.5, col=2, lty=2)
    axis(side=2, at=0.5, label="0.5")
    dev.off()
    
    # trt effect bias
    trtMats <- sapply(post.res, function(x)x$mtrt)
    errs <- abs(trtMats - b)
    trtRes <- rbind(CI.fn(errs[1, ]), CI.fn(errs[2, ]))
    rownames(trtRes) <- c("Borrow", "Non-borrow")
    trtRes
    }    
    
}

#5. Varying xis
{
    typ <- "LinearSwNewAlp"
    diagH <- rep(0.05, 4)
    x.tps <- c(2, 2, 2, 2)
    hs <- rep(1.1, 4)
    idx.tps <- paste0(x.tps, collapse = "")
    lam.q <- 0.10
    lam <- 30
    b <- 2
    phi0 <- 3
    invgam2 <- 100
    nSimu  <- 2000
    #allocation ratio vs Var(Y|X)/Var(Y)
    N <- 150
    Z1s <- c()
    bias.b <- c()
    bias.no <- c()
    mse.b <- c()
    mse.no <- c()
    Z1s.CIs <- c()
    xi.vs <- c(seq(0, 1, 0.2), 1.3, 1.6, 1.9, 2, 3, 5, 7, 9)
    #xi.vs <- c(0:3, 5, 7, 9)
    xiss <- list()
    for (ii in 1:length(xi.vs)){
        xiss[[ii]] <- c(0, 0, 0, 0) + xi.vs[ii]
        
    }
    for (idx in 1:length(xiss)){
        xis <- xiss[[idx]]
        f.name <- paste0("./results/", typ, "-b-", b, "-N-", N, "-lam-", lam, "-lamq-", 100*lam.q, "-phi0-", phi0, "-invgam2-", invgam2, 
                         "-H-", vec2code(diagH, 100), "-h-", vec2code(hs, 100), "-tps-", idx.tps, "-xis-", vec2code(xis, 10), "-nSimu-", nSimu, ".RData")
        load(f.name)
        kpidxs <- sapply(post.res, function(res)length(res$res$phi0.tk)) != 100
        post.res <- post.res[kpidxs];print(length(post.res))
        Z1rates <- sapply(post.res, function(s.res)mean(s.res$data$Z[21:N]))
        Z1rates.no <- sapply(post.res, function(s.res)mean(s.res$data.no$Z[21:N]))
        CI.fn(Z1rates)
        Z1s <- c(Z1s, CI.fn(Z1rates)[2])
        Z1s.CIs <- rbind(Z1s.CIs, CI.fn(Z1rates))
        CI.fn(Z1rates.no)
        # trt effect bias
        trtMats <- sapply(post.res, function(x)x$mtrt)
        errs <- abs(trtMats - b)
        errs2 <- abs(trtMats - b)**2
        trtRes <- rbind(CI.fn(errs[1, ]), CI.fn(errs[2, ]))
        trtMSEs <- rbind(CI.fn(errs2[1, ]), CI.fn(errs2[2, ]))
        rownames(trtRes) <- c("Borrow", "Non-borrow")
        trtRes
        bias.b <- rbind(bias.b, trtRes[1, ])
        bias.no <- rbind(bias.no, trtRes[2, ])
        mse.b <- rbind(mse.b, trtMSEs[1, ])
        mse.no <- rbind(mse.no, trtMSEs[2, ])
    }
    
    png(filename = paste0("./plots/newAlp_s.png"), width=720, height=480)
    par(mfrow=c(1, 2))
    plot(xi.vs, bias.b[, 2], type="l", lwd=2, col=1, xlab="s", ylab='Absolute bias', lty=1, 
         ylim=c(min(bias.no, bias.b)-0.1, max(bias.no, bias.b)+0.1))
    polygon(c(xi.vs, rev(xi.vs)), c(bias.b[, 1], rev(bias.b[, 3])),  col= rgb(0, 0, 0,0.1), border=NA)
    lines(xi.vs, bias.no[, 2], lwd=2, col=2, lty=2)
    #polygon(c(xi.vs, rev(xi.vs)), c(bias.no[, 1], rev(bias.no[, 3])),  col= rgb(1, 0, 0,0.1), border=NA)
    legend("topright", c("Borrow", "no Borrow"), col=c(1, 2), lty=c(1, 2))
    
    plot(xi.vs, Z1s, type="l", ylab="Allocation rate", xlab="s")
    abline(h=0.5, col=2, lty=2)
    #axis(side=2, at=0.5, label="0.5") 
    dev.off()
    
    # plot(xi.vs, mse.b[, 2], type="l", lwd=2, col=1, xlab="xis", ylab='MSE', lty=1, 
    #      ylim=c(min(mse.no, mse.b)-0.1, max(mse.no, mse.b)+0.1))
    # polygon(c(xi.vs, rev(xi.vs)), c(mse.b[, 1], rev(mse.b[, 3])),  col= rgb(0, 0, 0,0.1), border=NA)
    # lines(xi.vs, mse.no[, 2], lwd=2, col=2, lty=2)
    # #polygon(c(xi.vs, rev(xi.vs)), c(mse.no[, 1], rev(mse.no[, 3])),  col= rgb(1, 0, 0,0.1), border=NA)
    # legend("topright", c("Borrow", "no Borrow"), col=c(1, 2), lty=c(1, 2))
}

#t(round(Z1s.CIs*100, 2))

#test 
{
typ <- "LinearSwNewAlp"
diagH <- rep(0.05, 4)
x.tps <- c(2, 2, 2, 2)
hs <- rep(1.1, 4)
xis <- c(0, 3, 3, 0)
idx.tps <- paste0(x.tps, collapse = "")
lam.q <- 0.10
lam <- 30
b <- 2
phi0 <- 3
invgam2 <- 100
nSimu  <- 1000
#allocation ratio vs Var(Y|X)/Var(Y)
N <- 150
Z1s <- c()
bias.b <- c()
bias.no <- c()
invgam2s <- c(0.5, 1, 2, 4, 15, 30, 50, 100)
invgam2s <- c(30, 50, 100)
#for (invgam2 in invgam2s){
f.name <- paste0("./results/", typ, "-b-", b, "-N-", N, "-lam-", lam, "-lamq-", 100*lam.q, "-phi0-", phi0, "-invgam2-", invgam2, 
                     "-H-", vec2code(diagH, 100), "-h-", vec2code(hs, 100), "-tps-", idx.tps, "-xis-", vec2code(xis, 10), "-nSimu-", nSimu, ".RData")
load(f.name)
kpidxs <- sapply(post.res, function(res)length(res$res$phi0.tk)) != 100
post.res <- post.res[kpidxs];print(length(post.res))
Z1rates <- sapply(post.res, function(s.res)mean(s.res$data$Z[21:N]))
Z1rates.no <- sapply(post.res, function(s.res)mean(s.res$data.no$Z[21:N]))
CI.fn(Z1rates)
Z1s <- c(Z1s, CI.fn(Z1rates)[2])
CI.fn(Z1rates.no)
# trt effect bias
trtMats <- sapply(post.res, function(x)x$mtrt)
errs <- abs(trtMats - b)
trtRes <- rbind(CI.fn(errs[1, ]), CI.fn(errs[2, ]))
rownames(trtRes) <- c("Borrow", "Non-borrow")
trtRes
bias.b <- c(bias.b, trtRes[1, 2])
bias.no <- c(bias.no, trtRes[2, 2])
#}

# par(mfrow=c(1, 2))
# plot(log(invgam2s), Z1s, type="l")
# plot(log(invgam2s), bias.b, type="l")
# lines(log(invgam2s), bias.no, type="l", lty=2, col=2)

ext.sub.info.fn <- function(res){
        subId <- res$data$X1*2 + res$data$X2 + 1
        res.info <- list(Z=res$data$Z, subId=subId, tau2s=res$tau2s)
        Z1ms <- sapply(1:4, function(i) {mean(res.info$Z[21:N][res.info$subId[21:N]==i])})
        tau2ms <-  sapply(1:4, function(i) {mean(res.info$tau2s[res.info$subId==i])})
        rv <- list(Z1ms=Z1ms, tau2ms=tau2ms)
        rv
    }    
    
png(filename = paste0("./plots/newAlp_subgroup.png"), width=720, height=480)
par(mfrow=c(1, 2))
tau2mss <- lapply(post.res, function(res)ext.sub.info.fn(res$res)$tau2ms)
tau2Mat <- do.call(rbind, tau2mss)
colMeans(tau2Mat)
median(tau2Mat[, 1])
median(tau2Mat[, 2])
boxplot(tau2Mat, ylab="Tau2", xlab="Subgroup Index")
    
Z1mss <- lapply(post.res, function(res)ext.sub.info.fn(res$res)$Z1ms)
Z1Mat <- do.call(rbind, Z1mss)
colMeans(Z1Mat)
boxplot(Z1Mat, ylab="Percentage of treatment group", xlab="Subgroup Index")
abline(h=0.5, col=2, lty=2)
axis(side=2, at=0.5, label="0.5")
dev.off()
    
}
