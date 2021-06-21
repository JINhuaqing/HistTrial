# This script to test my function
# RPS and RBC designs are correct
# Following the procedure of (Jiang et al, Stat sinica, 2018)
source("utils.R")

gen.X.cts <- function(p){
    #args:
    #   p: num of covariances
    epss <- sapply(1:p, function(k){rnorm(1, k/2, 5)})
    x <- 2*exp(epss)/(1+exp(epss)) - 1
    names(x) <- paste0("X", 1:p)
    x.b <- sign(x)
    x.b[x.b==-1] <- 0
        
    return(list(x=x, x.b=x.b))
}


F.score.fn <- function(Xs, Z){
    n <- length(Z)
    n0 <- n - sum(Z)
    n1 <- sum(Z)
    p <- dim(Xs)[2]
    mMat <- matrix(rep(colMeans(Xs), n), nrow=n, byrow=TRUE)
    ZMat <- matrix(rep(Z, p), nrow=n, byrow=FALSE)
    
    SSTs <- colSums((Xs-mMat)**2)
    SSBs <- n0 * (colSums((1-ZMat)*Xs)/n0- colSums(Xs)/n)**2 + n1 * (colSums(ZMat*Xs)/n1- colSums(Xs)/n)**2
    Fs <- SSBs/(SSTs-SSBs) * (n-2)
    list(Fm = mean(Fs), Fs=Fs)
    
}

p <- 5

nRep <- 1000

nMax <- 50

XsBs <- list()
XsFs <- list()
Zs <- list()
for (j in 1:nRep){
    print(j)
    dat <- gen.X.cts(p)
    Xsb <- matrix(dat$x.b, nrow=1)
    XsF <- matrix(dat$x, nrow=1)
    Z.PS <- rbinom(1, 1, 0.5)
    Z.WPS <- rbinom(1, 1, 0.5)
    Z.BC <- rbinom(1, 1, 0.5)
    Z.WBC <- rbinom(1, 1, 0.5)
    for (i in 2:nMax){
        cdat <- gen.X.cts(p)
        
        ass.PS.res <- RPS.design(cdat$x.b, Xsb, Z.PS, hs=rep(0.1, p))
        ass.WPS.res <- RPS.design(cdat$x, XsF, Z.WPS, hs=rep(2.1, p))
        ass.BC.res <- RBC.design(cdat$x.b, Xsb, Z.BC, hs=rep(0.1, p))
        ass.WBC.res <- RBC.design(cdat$x, XsF, Z.WBC, hs=rep(2.1, p))
        
        Z.PS <- c(Z.PS, ass.PS.res$grp-1)
        Z.WPS <- c(Z.WPS, ass.WPS.res$grp-1)
        Z.BC <- c(Z.BC, ass.BC.res$grp-1)
        Z.WBC <- c(Z.WBC, ass.WBC.res$grp-1)
        Xsb <- rbind(Xsb, cdat$x.b)
        XsF <- rbind(XsF, cdat$x)
    }
    Zs[[j]] <- list(PS=Z.PS, BC=Z.BC, WPS=Z.WPS, WBC=Z.WBC)
    XsBs[[j]] <- Xsb
    XsFs[[j]] <- XsF
}


nams <- c("PS", "WPS", "BC", "WBC")
m1s <- c()

for (nam in nams){
    ZsMat <- do.call(rbind, lapply(Zs, function(z)z[[nam]]))
    m1 <- mean(abs(2*rowSums(ZsMat) - 50))
    m1s <- c(m1s, m1)
}
names(m1s) <- nams



Fsss <- list()
for (nam in nams){
    Fss <- list()
    for (j in 1:nRep){
        Z <- Zs[[j]][[nam]]
        XsF <- XsFs[[j]]
        Fss[[j]] <- F.score.fn(XsF, Z)$Fs
    
    }
    Fsss[[nam]] <- Fss
    
}

rZs <- lapply(1:nRep, function(i)rbinom(nMax, 1, 0.5))
FssR <- list()
for (j in 1:nRep){
    rZ <- rZs[[j]]
    XsF <- XsFs[[j]]
    FssR[[j]] <- F.score.fn(XsF, rZ)$Fs
}

FsRMat <- do.call(rbind, FssR)
j <- 5
plot(density(FsRMat[, j]), lty=1, lwd=2, col=1, ylim=c(0, 6), xlim=c(0, 5))
for (k in 1:4){
    nam <- nams[k]
    Fss <- Fsss[[nam]]
    FsMat <- do.call(rbind, Fss)
    print(mean(FsMat))
    lines(density(FsMat[, j]), lty=2, lwd=2, col=k+1)
}
legend("topright", nams, lty=2:5, col=2:5)





