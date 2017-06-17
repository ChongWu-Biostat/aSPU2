### SSU #################

SumSqU<-function(U, CovS){
    if (is.null(dim(CovS))) {# only one-dim:
        Tscore<- sum(U^2 /CovS)
        if (is.na(Tscore) || is.infinite(Tscore) || is.nan(Tscore)) Tscore<-0
        pTg1<-as.numeric(1-pchisq(Tscore, 1))
    }
    else {
        #it's possible U=0 and Cov(U)=0:
        if (all(abs(U)<1e-20)) pTg1<-1 else{
            Tg1<- t(U) %*% U
            ##distr of Tg1 is sum of cr Chisq_1:
            cr<-eigen(CovS, only.values=TRUE)$values
            ##approximate the distri by alpha Chisq_d + beta:
            alpha1<-sum(cr*cr*cr)/sum(cr*cr)
            beta1<-sum(cr) - (sum(cr*cr)^2)/(sum(cr*cr*cr))
            d1<-(sum(cr*cr)^3)/(sum(cr*cr*cr)^2)
            alpha1<-as.double(alpha1)
            beta1<-as.double(beta1)
            d1<-as.double(d1)
            pTg1<-as.numeric(1-pchisq((Tg1-beta1)/alpha1, d1))
        }
    }
    return(pTg1)
}

##########SumTest########################
Sum<-function(U, CovS){
    #it's possible U=0 and Cov(U)=0:
    if (all(abs(sum(U))<1e-20)) pTsum<-1 else{
        a<-rep(1, length(U))
        Tsum<- sum(U)/(sqrt(as.numeric(t(a) %*% CovS %*% (a))))
        pTsum <- as.numeric( 1-pchisq(Tsum^2, 1) )
    }
    pTsum
}

##########UminP Test########################
UminPd<-function(U, CovS){
    
    if (is.null(dim(CovS))) {# only one-dim:
        Tu<- sum(U^2 /CovS)
        if (is.na(Tu) || is.infinite(Tu) || is.nan(Tu)) Tu<-0
        pTu<-as.numeric(1-pchisq(Tu, 1))
    }
    else{
        ####it's POSSIBLR Ui=0 and CovS[i,i]=0!
        Tu<-as.vector(abs(U)/(sqrt(diag(CovS)) + 1e-20) )
        k<-length(U)
        V <- matrix(0,nrow=k, ncol=k)
        for(i in 1:k){
            for(j in 1:k){
                if (abs(CovS[i,j])>1e-20)
                V[i,j] <- CovS[i,j]/sqrt(CovS[i,i]*CovS[j,j])
                else   V[i,j] <- 1e-20
            }
        }
        pTu <- as.numeric(PowerUniv(Tu,V))
    }
    
    pTu
}


PowerUniv <- function(U,V){
    n <- dim(V)[1]
    
    x <- as.numeric(max(abs(U)))
    TER <- as.numeric(1-pmvnorm(lower=c(rep(-x,n)),upper=c(rep(x,n)),mean=c(rep(0,n)),sigma=V))
    
    return(TER)
}

combined.score <- function(U,V,weight){
    weight = as.matrix(weight)
    cov = t(weight) %*% V %*% weight
    
    df = qr(cov)$rank
    T = t(t(weight) %*% U ) %*% solve(cov) %*% t(weight) %*% U
    
    pvalue = 1 - pchisq(T, df = df)
    return(list(T = T, pvs = pvalue))
}




# sim.aSPUw1, sim.aSPUw2 are for the situations with the number of permutations less than half million
# sim.aSPUw3, sim.aSPUw4 are for the situations with the number of permutations large than half million
sim.aSPUw2 <- function(U,V,weight, pow1 = c(1:8,Inf),pow2 = c(1:8,Inf),n.perm = 1000){
    
    Ts <- rep(0, length(pow1) * length(pow2))
    Ts.tmp <- matrix(0,dim(weight)[2],length(pow1))
    for(i in 1:dim(weight)[2]) {
        for(j in 1:length(pow1)) {
            if (pow1[j] < Inf) {
                a <- sum( (U * weight[,i])^pow1[j] )
                Ts.tmp[i,j] <- sign(a) * ( (abs(a))^(1/pow1[j]) )
            } else {
                Ts.tmp[i,j] <- max(abs(U * weight[,i]))
            }
        }
    }
    
    for(i in 1:length(pow1)) {
        for(j in 1:length(pow2)) {
            if (pow2[j] < Inf) {
                Ts[(j-1) * length(pow1) + i] <- sum(Ts.tmp[,i]^pow2[j])
            } else {
                Ts[(j-1) * length(pow1) + i] <- max(abs(Ts.tmp[,i]))
            }
        }
    }
    
    eV <- eigen(V)
    eV$values[ eV$values < 0 ] = 0
    
    CovSsqrt <- t(eV$vectors %*% (t(eV$vectors) * sqrt(eV$values)))
    pow1[pow1==Inf] <- 0 # pass 0 as infitiy
    pow2[pow2==Inf] <- 0 # pass 0 as infitiy
    
    T0sC <- calcT0Wsim2(as.matrix(CovSsqrt),as.matrix(weight), as.vector(pow1), as.vector(pow2), n.perm)
    
    # We declare the variables before using it.
    T0s <- T0sC$T0s
    T0s[1:5,]
    pPerm0 <- rep(NA,length(pow1) * length(pow2))
    P0s <- rep(NA,n.perm)
    minp0 <- rep(NA,n.perm)
    
    T0s.abs <- abs(T0s)
    Ts.abs <- abs(Ts)
    
    for ( j in 1:(length(pow1)* length(pow2))) {
        pPerm0[j] <- sum( Ts.abs[j] <= T0s.abs[,j] ) / n.perm
        P0s <- ( ( n.perm - rank( T0s.abs[,j] ) ) + 1 ) / (n.perm)
        if (j == 1 ) {
            minp0 <- P0s
        } else {
            minp0[which(minp0>P0s)] <- P0s[which(minp0>P0s)]
        }
    }
    
    Paspu <- (sum(minp0 <= min(pPerm0)) + 1) / (n.perm + 1)
    pvs <- c(pPerm0, Paspu)
    
    Ts <- c(Ts, min(pPerm0))
    
    if(min(pow1) == 0) {
        pow1[which(pow1 == 0 )] = Inf
    }
    
    if(min(pow2) == 0) {
        pow2[which(pow2 == 0)] = Inf
    }
    
    nmvec <- NULL;
    for(ii in pow2) {
    	   for(jj in pow1) {
               nmvec <- c(nmvec, paste("SPU(",jj,",",ii,")",sep=""))
           }
    }
    
    nmvec <- c(nmvec, "daSPU")
    names(pvs) <- nmvec
    names(Ts) <- nmvec
    
    return(pvs)
}




sim.aSPUw4 <- function(U,V,weight, pow1 = c(1:8,Inf),pow2 = c(1:8,Inf),n.perm = 1000){
    
    Ts <- rep(0, length(pow1) * length(pow2))
    Ts.tmp <- matrix(0,dim(weight)[2],length(pow1))
    for(i in 1:dim(weight)[2]) {
        for(j in 1:length(pow1)) {
            if (pow1[j] < Inf) {
                a <- sum( (U * weight[,i])^pow1[j] )
                Ts.tmp[i,j] <- sign(a) * ( (abs(a))^(1/pow1[j]) )
            } else {
                Ts.tmp[i,j] <- max(abs(U * weight[,i]))
            }
        }
    }
    
    for(i in 1:length(pow1)) {
        for(j in 1:length(pow2)) {
            if (pow2[j] < Inf) {
                Ts[(j-1) * length(pow1) + i] <- sum(Ts.tmp[,i]^pow2[j])
            } else {
                Ts[(j-1) * length(pow1) + i] <- max(abs(Ts.tmp[,i]))
            }
        }
    }
    
    eV <- eigen(V)
    eV$values[ eV$values < 0 ] = 0
    
    CovSsqrt <- t(eV$vectors %*% (t(eV$vectors) * sqrt(eV$values)))
    pow1[pow1==Inf] <- 0 # pass 0 as infitiy
    pow2[pow2==Inf] <- 0 # pass 0 as infitiy
    
    len.pow = length(pow1) * length(pow2)
    
    T0s <- big.matrix(n.perm,len.pow,init = 0.0,type = "double")
    
    Ts.abs <- abs(Ts)
    
    
    pvalue = daSPU_calcT0(as.matrix(CovSsqrt),as.matrix(weight), as.vector(pow1), as.vector(pow2), as.vector(Ts.abs), T0s@address, n.perm)
    
    pvalue = as.vector(pvalue)
    # We declare the variables before using it.
    
    if(min(pow1) == 0) {
        pow1[which(pow1 == 0 )] = Inf
    }
    
    if(min(pow2) == 0) {
        pow2[which(pow2 == 0)] = Inf
    }
    
    nmvec <- NULL;
    for(ii in pow2) {
    	   for(jj in pow1) {
               nmvec <- c(nmvec, paste("SPU(",jj,",",ii,")",sep=""))
           }
    }
    
    nmvec <- c(nmvec, "daSPU")
    names(pvalue) <- nmvec
    
    pvalue
}


