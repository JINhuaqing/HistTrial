CI.fn <- function(errs){
    m.v <- mean(errs)
    sd.v <- sd(errs)
    se <- sd.v/sqrt(length(errs))
    low <- m.v - 1.96*se
    up <- m.v + 1.96*se
    rvs <- c(low, m.v, up)
    names(rvs) <- c("lower", "mean", "upper")
    rvs
}


post.fn <- function(trts, dlt0=1){
    prob <- mean(trts > dlt0)
    prob
}

post.b.fn <- function(res, dlt0=1){
    sps.trts <- res$sps.trts
    trtsMat <- do.call(rbind, sps.trts)
    prb <- post.fn(trtsMat[, 1], dlt0)
    prb.no <- post.fn(trtsMat[, 2], dlt0)
    c(prb, prb.no)
}

post.b.fn2 <- function(res, dlt0=1){
    sps.trts <- res$sps.trts
    trtsMat <- sps.trts
    #trtsMat <- abs(sps.trts)
    prb <- post.fn(trtsMat[, 1], dlt0)
    prb.no <- post.fn(trtsMat[, 2], dlt0)
    c(prb, prb.no)
}

eps.alig <- function(probs, aprob=0.05){
    can.epss <- seq(0, 1, 0.001)
    tyI.errs <- sapply(can.epss, function(can.eps){mean(probs>can.eps)})
    idx <- which.min(abs(tyI.errs - aprob))
    eps0 <- can.epss[idx]
    tyI.err <- tyI.errs[idx]
    list(eps0=eps0, tyI.err=tyI.err)
}


sub.Z1.fn <- function(data){
    data$subId <- data$X1*2 + data$X2 + 1
    res <- list(Z=data$Z, subId=data$subId)
    res
}

sub.Z1rate.fn <- function(res){
    z1ratesM <- c()
    for (idx in 1:4){
        z1ratesM[idx] <- mean(res$Z[res$subId==idx])
    }
}

vec2code <- function(vec,fct=10){
    vecF <- vec * fct
    pre <- c()
    if (fct == 10){
        for (ivecF in vecF) {
            if (ivecF < fct){
                pre <- c(pre, "0")
            }else{
                pre <- c(pre, "")
            }
        }
    }else if (fct==100){
        for (ivecF in vecF) {
            if (ivecF >= 100){
                pre <- c(pre, "")
            }else if(ivecF >=10){
                pre <- c(pre, "0")
            }else{
                pre <- c(pre, "00")
            }
        }
    }
    idxCode <- paste0(paste0(pre, vecF), collapse="")
    idxCode
}