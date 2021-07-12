setwd("C:/Users/JINHU/Documents/ProjectCode/HistTrial")
library(magrittr)
source("util_ana.R")



# for hypothesis test
dlt0 <- 1

typ <- "Simplex-Linear-diff-0"
N <- 100
lam <- 1
b <- 3
phi0 <- 1

#allocation ratio vs Var(Y|X)/Var(Y)
#for (phi0.it in c(1, 2, 3)){
    f.name <- paste0("./results/", typ, "-b-", b, "-N-", N, "-lam-",  lam, "-phi0-", phi0, ".RData")
    load(f.name)
#    varYs <- sapply(post.res, function(s.res)var(s.res$data$Y))
    Z1rates <- sapply(post.res, function(s.res)mean(s.res$data$Z))
    post.res[[3]]$res0$tau2s
#    print(mean(Z1rates))
#    
#}

#dir("./results/", f.name)
f.name <- paste0("./results/", typ, "-b-", b, "-N-", N, "-lam-",  lam, "-phi0-", phi0, ".RData")

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
    f.name <- paste0("./results/", typ, "-b-", b, "-N-", N, "-lam-",  lam, "-phi0-", phi0, ".RData")
    load(f.name)
    
    probMat <- sapply(post.res, post.b.fn, dlt0=dlt0)
    power <- mean(probMat[1, ] > eps0)
    power.no <- mean(probMat[2, ] > eps0.no)
    powers[[i]] <- c(power, power.no)
    
    trtMats <- sapply(post.res, function(x)x$mtrt)
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
png(filename = paste0("./plots/bias_", typ, ".png"), width=360, height=480)
plot(bs, res.b$biasMean, type="l", lwd=2, col=1, xlab="b", ylab='Absolute bias', lty=1, ylim=c(0.25, 0.55))
polygon(c(bs, rev(bs)), c(res.b$biasLow, rev(res.b$biasUp)),  col= rgb(0, 0, 0,0.1), border=NA)
lines(bs, res.no$biasMean, lwd=2, col=2, lty=2)
polygon(c(bs, rev(bs)), c(res.no$biasLow, rev(res.no$biasUp)),  col= rgb(1, 0, 0,0.1), border=NA)
legend("bottomright", c("Borrow", "no Borrow"), col=c(1, 2), lty=c(1, 2))
dev.off()


png(filename = paste0("./plots/power_", typ, ".png"), width=360, height=480)
plot(bs, res.b$power, type="l", lwd=2, col=1, xlab="b", ylab='Power', lty=1, ylim=c(0, 1))
lines(bs, res.no$power, lwd=2, col=2, lty=2)
legend("bottomright", c("Borrow", "no Borrow"), col=c(1, 2), lty=c(1, 2))
dev.off()

png(filename = paste0("./plots/Z1rate_", typ, ".png"), width=360, height=480)
plot(bs, res.b$Z1rateMean, type="l", lwd=2, col=1, xlab="b", ylab='Percentage of treatment group', lty=1, ylim=c(0.4, 0.6))
polygon(c(bs, rev(bs)), c(res.b$Z1rateLow, rev(res.b$Z1rateUp)),  col= rgb(0, 0, 0,0.1), border=NA)
lines(bs, res.no$Z1rateMean, lwd=2, col=2, lty=2)
legend("bottomright", c("Borrow", "no Borrow"), col=c(1, 2), lty=c(1, 2))
dev.off()
