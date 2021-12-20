setwd("C:/Users/JINHU/Documents/ProjectCode/HistTrial/realData")

load("dat.merge.RData")

library(dplyr)

## 1. Data prepocessing
# 1.1 remove the SOF study and only keep race ==1. 
RCT.data <- filter(dat.merge, STUDY!="SOF" & race==1)
summary(RCT.data)
str(RCT.data)
RCT.data$STUDY <- droplevels(RCT.data$STUDY)
RCT.data$race <- droplevels(RCT.data$race)
str(RCT.data)
summary(RCT.data)

group_by(RCT.data, STUDY) %>% 
summarise(
    isNA=sum(is.na(VFXINC_36)), 
    isNoTNA=sum(!is.na(VFXINC_36))
                                        )
group_by(RCT.data, STUDY) %>% 
summarise(
    isNA=sum(is.na(VFXINC_24)), 
    isNOTNA=sum(!is.na(VFXINC_24))
                                        )

group_by(RCT.data, STUDY) %>% 
summarise(
    isNA=sum(is.na(dxhhp_24)), 
    isNOTNA=sum(!is.na(dxhhp_24))
                                        )
group_by(RCT.data, STUDY) %>% 
summarise(
    is0=sum(nfrxvert==0, na.rm = NA), 
    isNot0=sum(nfrxvert!=0, na.rm=NA))


plot(RCT.data$dxhhp_24[!is.na(RCT.data$dxhhp_24)])
plot(RCT.data$VFXINC_36[!is.na(RCT.data$VFXINC_36)])

# so I use logit(dxhhp_24) as a linear model
# Treatment variable: TRTN
# Covariates
# 1. age: interger
# 2. falls: 1 or 0
# 3. nfrxvert: interger
# 4. frxvert: 1 or 0 
# 5. menyrs: decimal


# 1.2 Remove obs with NA covariates
data.CL <- filter(RCT.data, (!is.na(falls)) & (!is.na(nfrxvert)) & (!is.na(frxvert)) & (!is.na(menyrs)))
data.CL <- transmute(data.CL, age=age, falls=falls,  frx=frxvert, menyrs=menyrs, nfrx=nfrxvert,
                     study=STUDY,  Y=log(dxhhp_24/(1-dxhhp_24)), Z=TRTN)

summary(data.CL)


## 2. Fit naive regression model

# I pool two historical datasets together 
# because, frx == 0 in FIT.clin trial while frx==1 in FIT.vert trial
data.Hist <- filter(data.CL, !study=="ZOL")
data.Cur <- filter(data.CL, study=="ZOL")




# 2.1 divide into 4 subgroups by frx and falls

sublinear.fit <- function(curData){
    curData$age <- as.vector(scale(curData$age))
    #curData$nfrx <- as.vector(scale(curData$nfrx))
    curData$menyrs <- as.vector(scale(curData$menyrs))
    data <- filter(curData, !is.na(Y))
    fit <- lm(Y~Z+age+menyrs, data=data)
    return(fit)
}

fit.Curs <- list()
fit.Curs[[1]] <- sublinear.fit(filter(data.Cur, falls==0 & frx==0))
fit.Curs[[2]] <- sublinear.fit(filter(data.Cur, falls==0 & frx==1))
fit.Curs[[3]] <- sublinear.fit(filter(data.Cur, falls==1 & frx==0))
fit.Curs[[4]] <- sublinear.fit(filter(data.Cur, falls==1 & frx==1))

eparas <- sapply(fit.Curs, function(x)x$coefficients);eparas

fit.Hists <- list()
fit.Hists[[1]] <- sublinear.fit(filter(data.Hist, falls==0 & frx==0))
fit.Hists[[2]] <- sublinear.fit(filter(data.Hist, falls==0 & frx==1))
fit.Hists[[3]] <- sublinear.fit(filter(data.Hist, falls==1 & frx==0))
fit.Hists[[4]] <- sublinear.fit(filter(data.Hist, falls==1 & frx==1))

sapply(fit.Hists, function(x)x$coefficients)

# 2.2 The true parameters of the current trial 
b <- mean(eparas[2, ])
betassMat <- t(rbind(eparas[c(1, 3), ], rep(0, 4), rep(0, 4), eparas[4, ]))
colnames(betassMat) <- c("intercept", "age", "falls", "frx", "menyrs")
betassMat
betass <- list(
   para1=betassMat[1, ],
   para2=betassMat[2, ],
   para3=betassMat[3, ],
   para4=betassMat[4, ]
)
betass
