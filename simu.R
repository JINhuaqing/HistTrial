#rm(list=ls())
library(magrittr)
library(dplyr)
setwd("C:/Users/JINHU/Documents/ProjectCode/HistTrial")
source("utils.R")



# parameters for current data 
betass <- list(para1=c(2, 1, -1, 3, -2), 
              para2=c(3, 1, 2, 4, 2), 
              para3=c(0, -1, -1, 3, 2), 
              para4=c(2, 0, -1, 4, -2)
                  )

# parameters for historical data
alpss <-  list(para1=c(2, 1, -1, 3, -2), 
             para3=c(3, 1, 2, 4, 2), 
             para2=c(0, -1, -1, 3, 2), # switch para2 and para3
             para4=c(2, 0, -1, 4, -2) )
b <- 2
phi0 = phi1 = 1
N <- 200 # total sample size
# parameters
lam <- 0.1
hs <- rep(2.1, 4)

tXs <- gen.Data.Xs(1000, x.tps)
H <- diag(c(bw.nrd(tXs[, 1]), bw.nrd(tXs[, 2]), bw.nrd(tXs[, 3]), bw.nrd(tXs[, 4])))

# initial dataset
n0 <- 10
x.tps <- c(2, 2, "c", "c")

nSimu <- 100
trt.effs <- list()
for (i in 1:nSimu){
    print(c(i, nSimu))
    source("whole_procedure.R")
    trt.effs[[i]] <- c(trt.eff, trt.eff.no)
}

trt.effs.mat <- do.call(rbind, trt.effs)
errs <- abs(trt.effs.mat - b)
plot(errs[, 1], errs[, 2], xlim=c(0, 2), ylim=c(0, 2), ylab="Borrow", xlab="no Borrow")
abline(a=0, b=1, col=2, lwd=2)
m.err <- colMeans(errs)
names(m.err) <- c("Borrow", "no Borrow")
m.err
err.se1 <- sd(errs[, 1])/sqrt(length(errs[, 1]))
err.se2 <- sd(errs[, 2])/sqrt(length(errs[, 2]))
err.CI <- c(m.err[1]-1.96*err.se1, m.err[1], m.err[1]+1.96*err.se2)
err.CI.no <- c(m.err[2]-1.96*err.se1, m.err[2], m.err[2]+1.96*err.se2)
CIs <- rbind(err.CI, err.CI.no)
colnames(CIs) <- c("Low", "Mean", "Up")
rownames(CIs) <- c("Borrow", "No Borrow")
CIs
