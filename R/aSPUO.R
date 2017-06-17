#' Adaptive Sum of Powered Score Omnibus tests (aSPUO) test for single trait with multiple weights for each genetic marker
#'
#' It returns p-values.
#'
#' @param U Score vector for the genetic marker set (p by 1 matrix).
#'
#' @param V Corresponding covariance matrix for the score vector (p by p matrix).
#'
#' @param weight Multiple weights for each genetic markers (p by p matrix).
#'
#' @param pow Gamma sets. Different gamma corresponds to different existing test. For example, gamma = 1 equals Sum test; gamma = 2 equals SSU or SKAT with linear kernel; gamma = Inf is similar to minimum p value test. We recommend use pow = c(1:6, Inf) to maintain high power under various scenarios, although other choice may be slightly powerful in some specific situations.
#'
#' @param n.perm number of permutations.
#'
#' @export
#'
#' @return  It returns a list, which contains the following two information:
#'          asy-SPU(1)-O: Theoretical p value for SPU(1)-O. This p value is exactly same as TWAS-Omnibus.
#'          sim-pvalue: a matrix containing p values for SPU/aSPU with different weights, and p values for SPU-O and aSPU-O.
#'
#' @author Chong Wu, Zhiyuan Xu, and Wei Pan

aSPUO <- function(U,V,weight, pow = c(1:8,Inf),n.perm = 1000) {
    
    weight = as.matrix(weight)
    
    # remove SNPs corresponding to zero weight
    weight.tmp = abs(weight)
    index = rowSums(weight.tmp) > 0
    U = as.matrix(U)
    U = U[index,]
    V = V[index,]
    V = V[,index]
    weight = weight[index,]
    
    weight = as.matrix(weight)
    
    SPU1o = combined.score(U,V,weight)
    
    ## calculate p vlaue
    npow = length(pow)
    nweight = dim(weight)[2]
    Ts <- rep(0,npow * nweight)
    
    for(i in 1:nweight) {
        for(j in 1:npow) {
            if (pow[j] < Inf) {
                Ts[(j-1) * nweight + i] <- sum((U * weight[,i])^pow[j])
            } else {
                Ts[(j-1) * nweight + i] <- max(abs(U * weight[,i]))
            }
        }
    }
    
    eV <- eigen(V)
    eV$values[ eV$values < 0 ] = 0
    
    CovSsqrt <- t(eV$vectors %*% (t(eV$vectors) * sqrt(eV$values)))
    pow[pow==Inf] <- 0 # pass 0 as infitiy
    
    T0s.size = npow * nweight
    T0s <- big.matrix(n.perm,T0s.size,type = "double")
    T1s <- big.matrix(n.perm,npow,type = "double")
    
    minp0_sign <- big.matrix(n.perm,nweight,type = "double")
    Ts.abs <- abs(Ts)


    res.tmp = calcT0Wsim4(as.matrix(CovSsqrt), as.matrix(weight),as.matrix(Ts.abs), as.matrix(pow),T0s@address,minp0_sign@address,n.perm)
    
    cov.res = res.tmp$cov
    ########################################################
    ###     Calculate p value of aSPU for every weights  ###
    ########################################################
    
    final.pvalue = matrix(NA,npow + 1,nweight + 1)
    final.pvalue[,1:nweight] = res.tmp$pPerm0
    
    ## Dr.Pan idea
    # calculate test statistics and its corresponding distribution
    Ts.abs <- abs(Ts)
    
    T2s <- big.matrix(n.perm,T0s.size,type = "double")
    
    Res = calc_p_pan(as.matrix(Ts.abs),T0s@address,T2s@address, minp0_sign@address,T0s.size,nweight, n.perm)
    
    Ts2 = Res$pPerm0
    cov.res = Res$cov
    
    for(i in 1:npow) {
        tmp = cov.res[((i-1) * nweight + 1):(i * nweight),((i-1) * nweight + 1):(i * nweight)]
        cov.res[((i-1) * nweight + 1):(i * nweight),((i-1) * nweight + 1):(i * nweight)] = ginv(tmp)
    }
    
    Ts2 = qnorm(1- Ts2/2) * rep(sign(Ts[1:nweight]),npow)
    
    rm(T0s)
    calc_test_pan(cov.res,as.matrix(weight),as.matrix(pow),T1s@address,T2s@address, n.perm)
    
    Ts.chong = calc_test_ts(cov.res, as.matrix(weight), as.matrix(pow), as.matrix(Ts2))
    
    Res = calc_p_ch_pan(as.matrix(pow),as.matrix(Ts.chong),T1s@address, n.perm)
    
    minp0 = Res$minp0
    pPerm0 = Res$pPerm0
    
    Paspu <- (sum(minp0 <= min(pPerm0)) + 1) / (n.perm+1)
    pvs <- c(pPerm0, Paspu)
    pow[pow==0] <- Inf
    
    final.pvalue[,dim(final.pvalue)[2]] = pvs
    rownames(final.pvalue) <- c(paste("SPU(", pow,")", sep=""), "aSPU")
   
    if( is.null(rownames(weight)) ) {
        colnames(final.pvalue) <- c(paste0("weight_",1:nweight),"SPU-O")

    } else {
        colnames(final.pvalue) <- c(paste0("weight_",colnames(weight)),"SPU-O")
    }
    
    output = list("asy_SPU1O" = SPU1o,"sim_pvalue" = final.pvalue )
    return(output)
}
