library(magrittr)
source("util_ana.R")


# under H0
fil <-  dir("./results/", pattern="RealDataRSS-b-0.*lam-*", full.names=T);fil
load(fil)
post.res[[1]]
fadsf

kpidxs <- sapply(post.res, function(res)length(res$res$phi0.tk)) != 100
post.res <- post.res[kpidxs]
varYs <- sapply(post.res, function(s.res)var(s.res$data$Y))
Z1rates <- sapply(post.res, function(s.res)mean(s.res$data$Z[21:N]))
Z1rates.no <- sapply(post.res, function(s.res)mean(s.res$data.no$Z[21:N]))
mean(Z1rates)
mean(Z1rates.no)
Z1RatesCI <- rbind(CI.fn(Z1rates), CI.fn(Z1rates.no))
rownames(Z1RatesCI) <- c("Borrow", "No Borrow")
Z1RatesCI
trtMats <- sapply(post.res, function(x)x$mtrt)
errs <- abs(trtMats - b)
trt.Errs <- rbind(CI.fn(errs[1, ]), CI.fn(errs[2, ]))
rownames(trt.Errs) <- c("Borrow", "No Borrow")
trt.Errs
ext.sub.info.fn <- function(res){
        subId <- res$data$frx*2 + res$data$falls+ 1
        res.info <- list(Z=res$data$Z, subId=subId, tau2s=res$tau2s)
        Z1ms <- sapply(1:4, function(i) {mean(res.info$Z[21:N][res.info$subId[21:N]==i])})
        tau2ms <-  sapply(1:4, function(i) {mean(res.info$tau2s[res.info$subId==i])})
        rv <- list(Z1ms=Z1ms, tau2ms=tau2ms)
        rv
    }   
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

tres <- post.res[[2]]
tres$res$phi0
mean(tres$data$Z)
tres$Rs
tres$res$tau2.tds
tres$res$tau2s
tres$res$lam.tru

data <- tres$data
ass.res <- RBC.design(unlist(data[4, 3:(p+2)]), data[, 3:(p+2)], data$Z, hs, R=1)
xNew <- unlist(data[4, 3:(p+2)])
Xs <- data[, 3:(p+2)]
Z <- data$Z
ass.res
