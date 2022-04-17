library(magrittr)
source("util_ana.R")

N <- 150

Hs <- c()
Z1r.list <- list()
trt.Err.list <- list()
Z1Mats <- list()
fil <-  dir("./results/", pattern="RealData-b-3e-01-N-150-lam-200-lamq-15-phi0-3e-01-invgam2-0.33-H-010010999999-h-110110130130-tps-22cc-xis-.*-nSimu-.*", full.names=T);fil
fil <- sort(fil)
flag <- 1
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
        
    Z1r.list[[flag]] <- Z1RatesCI
    trt.Err.list[[flag]] <- trt.Errs
    Z1Mats[[flag]] <- colMeans(Z1Mat)
    flag <- flag + 1
}

Z1s <- sapply(Z1r.list, function(x)x[1, 2])
errs.no <- sapply(trt.Err.list, function(x)return(x[2, 2]))
errs.b <- sapply(trt.Err.list, function(x)return(x[1, 2]))
errs.low  <- sapply(trt.Err.list, function(x)return(x[1, 1]))
errs.up <- sapply(trt.Err.list, function(x)return(x[1, 3]))


par(mfrow=c(1, 2))
plot(Hs, errs.b, type="l", ylim=c(0.04, 0.055), xlab=expression(xi), ylab="Error", 
     lty=1, col="red", lwd=2)
polygon(c(Hs, rev(Hs)), c(errs.low, rev(errs.up)), col=rgb(1, 0, 0, 0.3), border=NA)
lines(Hs, errs.no, col="blue", lty=2, lwd=2)
legend("topright", c("Borrow", "No Borrow"), col=c("red", "blue"), lty=c(1, 2))

Z1MatsArr <- do.call(rbind, Z1Mats)
plot(Hs, Z1MatsArr[, 1], ylim=c(0.45, 0.57), col=1, lty=1, type="l", lwd=2, xlab=expression(xi), 
     ylab="Allocate rates")
lines(Hs, Z1MatsArr[, 2], col=2, lty=2, lwd=2)
lines(Hs, Z1MatsArr[, 3], col=3, lty=3, lwd=2)
lines(Hs, Z1MatsArr[, 4], col=4, lty=4, lwd=2)
lines(Hs, Z1s, col=5, lty=5, lwd=4)
abline(h=0.5, col=6, lty=2, lwd=1)
#axis(side=2, at=0.5, label="0.5")
legend("bottomleft", c(paste0("SubGroup", 1:4), "Overall"), col=1:5, lty=1:5, lwd=c(2,2,2,2, 4), ncol=1)




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
round(colMeans(tau2Mat), 2)
median(tau2Mat[, 1]) 
median(tau2Mat[, 2])
boxplot(tau2Mat, ylab="Tau2", xlab="Subgroup Index")

Z1mss <- lapply(post.res, function(res)ext.sub.info.fn(res$res)$Z1ms)
Z1Mat <- do.call(rbind, Z1mss)
colMeans(Z1Mat)
round(colMeans(Z1Mat), 3)
boxplot(Z1Mat, ylab="Percentage of treatment group", xlab="Subgroup Index")
abline(h=0.5, col=2, lty=2)
axis(side=2, at=0.5, label="0.5")

tres <- post.res[[2]]
tres$alpss
tres$betass
tres$res$phi0
mean(tres$data$Z)
tres$Rs
tres$lam.trus
tres$res$tau2.tds
tres$res$tau2s
tres$res$lam.tru
data <- tres$data
subId <- data$frx*2 + data$falls + 1
data$subId <- subId
data$tau2s <- tres$res$tau2s
data$tau2.tds <- tres$res$tau2.tds
data$Rs <- tres$Rs
group_by(data, subId) %>% summarise(
    Z=mean(Z), 
    Rs=mean(Rs), 
    tau2.tds = mean(tau2.tds),
    tau2s=mean(tau2s))
plot(subId, tres$Rs)



ass.res <- RBC.design(unlist(data[4, 3:(p+2)]), data[, 3:(p+2)], data$Z, hs, R=1)
xNew <- unlist(data[4, 3:(p+2)])
Xs <- data[, 3:(p+2)]
Z <- data$Z
ass.res
