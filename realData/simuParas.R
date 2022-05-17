#rm(list=ls())
#setwd("/home/huaqingj/MyResearch/HistTrial/")
setwd("/root/Documents/HQ/HistTrial")
#setwd("C:/Users/JINHU/Documents/ProjectCode/HistTrial/")

load("./realData/dat.merge.Rdata")
library(dplyr)

## 1. Data prepocessing
# 1.1 remove the SOF study and only keep race ==1. 
RCT.data <- filter(dat.merge, STUDY!="SOF" & race==1)
#summary(RCT.data)
#str(RCT.data)
RCT.data$STUDY <- droplevels(RCT.data$STUDY)
RCT.data$race <- droplevels(RCT.data$race)
#str(RCT.data)
#summary(RCT.data)

# group_by(RCT.data, STUDY) %>% 
# summarise(
#     isNA=sum(is.na(VFXINC_36)), 
#     isNoTNA=sum(!is.na(VFXINC_36))
#                                         )
# group_by(RCT.data, STUDY) %>% 
# summarise(
#     isNA=sum(is.na(VFXINC_24)), 
#     isNOTNA=sum(!is.na(VFXINC_24))
#                                         )
# 
# group_by(RCT.data, STUDY) %>% 
# summarise(
#     isNA=sum(is.na(dxhhp_24)), 
#     isNOTNA=sum(!is.na(dxhhp_24))
#                                         )
# group_by(RCT.data, STUDY) %>% 
# summarise(
#     is0=sum(nfrxvert==0, na.rm = NA), 
#     isNot0=sum(nfrxvert!=0, na.rm=NA))


#plot(RCT.data$dxhhp_24[!is.na(RCT.data$dxhhp_24)])
#plot(RCT.data$VFXINC_36[!is.na(RCT.data$VFXINC_36)])

# so I use dxhhp_24 as a linear model
# Treatment variable: TRTN (Z)
# Covariates
# 1. X1: dxhhp: continous
# 2. X2: falls: 1 or 0
# 3. X3: frxvert: 1 or 0 
# 4. X4: menyrs: decimal
# 5. nfrxvert: interger # I not use it 


# 1.2 Remove obs with NA covariates # I do not remove
#data.CL <- filter(RCT.data, (!is.na(falls)) & (!is.na(dxhhp_0)) & (!is.na(frxvert)) & (!is.na(menyrs)))
#data.CL <- filter(RCT.data, (!is.na(falls)) & (!is.na(nfrxvert)) & (!is.na(frxvert)) & (!is.na(menyrs)))
data.CL <- RCT.data
data.CL <- transmute(data.CL, age=age, falls=falls,  frx=frxvert, menyrs=menyrs, nfrx=nfrxvert, Y0=dxhhp_0, 
                     study=STUDY,  Y=(dxhhp_24-mean(dxhhp_24, na.rm=T))/sd(dxhhp_24, na.rm=T), Z=TRTN)

#summary(data.CL)


## 2. Fit naive regression model
# 2.1 divide into 4 subgroups by frx and falls
# We assume for different frx and falls, trt effect is the same
# while effects of age and menyrs will change 
# So model is 
# dxhhp_24 ~ beta0(X2, X3) + b *Z + beta1(X2, X3) X1 + beta4(X2, X3) X4 
data.CL$subGroupId <- as.factor(data.CL$frx*2 + data.CL$falls)

# I pool two historical datasets together 
# because, frx == 0 in FIT.clin trial while frx==1 in FIT.vert trial
data.Hist <- filter(data.CL, !study=="ZOL")
data.Cur <- filter(data.CL, study=="ZOL")

mean.menyrs<- mean(data.Cur$menyrs, na.rm=T)
sd.menyrs<- sd(data.Cur$menyrs, na.rm=T)
mean.Y0 <- mean(data.Cur$Y0, na.rm=T)
sd.Y0 <- sd(data.Cur$Y0, na.rm=T)

data.Cur$menyrs1 <- (data.Cur$menyrs - mean.menyrs)/sd.menyrs
data.Cur$Y01<- (data.Cur$Y0- mean.Y0)/sd.Y0

fit.all <- lm(Y~Z+subGroupId+subGroupId*menyrs1+subGroupId*Y01, data=data.Cur)
summary(fit.all)

coeffs <- fit.all$coefficients;coeffs
b <- coeffs[2] # trt effect
paras.Vec <- coeffs[-2][c(1:4, 5, 7:9, 6, 10:12)];paras.Vec
eparas <- matrix(paras.Vec, ncol=4, byrow=T)
eparas[1, 2:4] <- eparas[1, 2:4] + eparas[1, 1]
eparas[2, 2:4] <- eparas[2, 2:4] + eparas[2, 1]
eparas[3, 2:4] <- eparas[3, 2:4] + eparas[3, 1]
row.names(eparas) <- c("Intercept", "Menyrs", "dxhhp_0")
eparas






## historical 


hmean.menyrs<- mean(data.Hist$menyrs, na.rm=T)
hsd.menyrs<- sd(data.Hist$menyrs, na.rm=T)
hmean.Y0 <- mean(data.Hist$Y0, na.rm=T)
hsd.Y0 <- sd(data.Hist$Y0, na.rm=T)
data.Hist$Y01 <-(data.Hist$Y0- hmean.Y0)/hsd.Y0
data.Hist$menyrs1 <- (data.Hist$menyrs - hmean.menyrs)/hsd.menyrs

fitH.all <- lm(Y~Z+subGroupId+subGroupId*menyrs1+subGroupId*Y01, data=data.Hist)
summary(fitH.all)
coeffsH <- fitH.all$coefficients
bH <- coeffsH[2] # historical trt effect
parasH.Vec <- coeffsH[-2][c(1:4, 5, 7:9, 6, 10:12)];parasH.Vec
eparasH <- matrix(parasH.Vec, ncol=4, byrow=T)
eparasH[1, 2:4] <- eparasH[1, 2:4] + eparasH[1, 1]
eparasH[2, 2:4] <- eparasH[2, 2:4] + eparasH[2, 1]
eparasH[3, 2:4] <- eparasH[3, 2:4] + eparasH[3, 1]
row.names(eparasH) <- c("Intercept", "Menyrs", "dxhhp_0")
eparasH


# 2.2 The true parameters of the current trial 
betassMat <- t(rbind(eparas[1, ], rep(0, 4), rep(0, 4), eparas[c(2, 3), ]))
colnames(betassMat) <- c("intercept", "falls", "frx", "menyrs", "dxhhp_0")
betassMat
betass <- list(
   para1=betassMat[1, ],
   para2=betassMat[2, ],
   para3=betassMat[3, ],
   para4=betassMat[4, ]
)
betass # the true parameters
b # the true trt effect

# For historical parameters, alpss, I can generate them based on the current model
# or use estimated data
alpssMat <- t(rbind(eparasH[1, ], rep(0, 4), rep(0, 4), eparasH[c(2, 3), ]))
colnames(alpssMat) <- c("intercept", "falls", "frx", "menyrs", "dxhhp_0")
alpss <- list(
   para1=alpssMat[1, ],
   para2=alpssMat[2, ],
   para3=alpssMat[3, ],
   para4=alpssMat[4, ]
)
alpss # the true parameters

  

# 3. Obtain the distribution of X
library(ks)
fHats <- list()
for (i in 1:4){
    X0 <- filter(data.Cur, data.Cur$subGroupId==i-1) %>% transmute(menyrs=menyrs, Y0=Y0)
    nonNaIdxs <- rowSums(is.na(X0)) == 0
    fHat <- kde(as.matrix(X0[nonNaIdxs, ]))
    fHats[[i]] <- fHat
}

# function to generate data from real data
gen.Real.Xs <- function(n, fHats){
    maps <- list()
    maps[[1]] <- c(0, 0)
    maps[[2]] <- c(1, 0)
    maps[[3]] <- c(0, 1)
    maps[[4]] <- c(1, 1)
    
    subIds <- sample.int(4, n, replace=T) - 1
    sps <- c()
    for (i in 1:4){
        fHat <- fHats[[i]]
        nSps <- sum(subIds==i-1)
        if (nSps != 0){
            pSps <- matrix(rep(maps[[i]], nSps), ncol=2, byrow=T)
            curSps <- rkde(nSps, fHat)
            aSps <- cbind(pSps, curSps)
            sps <- rbind(sps, aSps)
        }
    }
    
    colnames(sps) <- c("falls", "frx", "menyrs", "Y0")
    sps[, 3] <- (sps[, 3] - mean.menyrs)/sd.menyrs
    sps[, 4] <- (sps[, 4] - mean.Y0)/sd.Y0
    rownames(sps) <- NULL
    sps
}


if (FALSE){
# 4. Check the true Ys and estimated Ys

sub.Paras.fn <- function(Xs, betass){
  #args:
  #   Xs: covariates, n x p 
  if (is.null(dim(Xs))){
    
      idxs <- Xs[1] +  Xs[2]*2 + 1
  }else{
      idxs <- Xs[, 1] +  Xs[, 2]*2 + 1
    
  }
  betMat <- do.call(rbind, betass[idxs])
  return(betMat)
  
}

curMean.fn <- function(Xs, Zs, betMat, b){
    #args:
    #  Xs: covariates matrix, n x p, data frame
    #  Zs: Group indicate function, n, vector
    #  betMat: matrix of beta parameters, n x p, matrix
    #  b: treatment effect, scalar
    if (is.null(dim(Xs))){
        fXs <- c(1, as.vector(Xs))
    }else{
        n <- dim(Xs)[1]
        fXs <- cbind(rep(1, n), Xs)
        fXs <- as.matrix(fXs)
    }
    mys <- rowSums(fXs * betMat) + b * Zs
    mys
}

ntest <- 1000
Xs <- gen.Real.Xs(ntest, fHats)

idx0 <- sample.int(ntest, size=floor(ntest/2))
Zs <- rep(1, ntest)
Zs[idx0] <- 0


betMat <- sub.Paras.fn(Xs, betass)
Ys.m <- curMean.fn(Xs, Zs, betMat, b)
nerrs <- rnorm(ntest, sd=sd(fit.all$residuals))
Ys <- Ys.m + nerrs

# # Two Densities is similar
plot(density(data.Cur$Y[data.Cur$Z==0], na.rm=T), main="Control Group")
lines(density(Ys[Zs==0], na.rm=T), col="red")
legend("topright", legend=c("True data", "Generated data"), col=c("black", "red"), lty=c(1, 1))
 
plot(density(data.Cur$Y[data.Cur$Z==1], na.rm=T), main="Treatment Group")
lines(density(Ys[Zs==1], na.rm=T), col="red")
legend("topright", legend=c("True data", "Generated data"), col=c("black", "red"), lty=c(1, 1))

}
