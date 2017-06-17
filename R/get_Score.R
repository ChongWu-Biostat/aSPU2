glm_get_sco <- function (Y, X, cov = NULL, model = "binomial") {
  
  model = match.arg(model)
  
  n <- length(Y)
  if (is.null(X) && length(X)>0) X=as.matrix(X, ncol=1)
  k <- ncol(X)
  
  if (is.null(cov)){
    ## NO nuisance parameters:
    Xg <- XUs <- X
    U <- t(Xg) %*% (Y-mean(Y))
    yresids <- Y-mean(Y)
    #  sigma0 = sqrt(sum(yresids^2)/(n-1))
    yfits <- rep(mean(Y), n)
    
    Xgb <- apply(X, 2, function(x)(x-mean(x)) )
    
    if( model == "binomial" ) {
      CovS <- mean(Y) * (1-mean(Y))  * (t(Xgb) %*% Xgb)
    } else {
       CovS <- var(Y)*(t(Xgb) %*% Xgb)
    }
  } else {
    ## with nuisance parameters:
    tdat1 <- data.frame(trait=Y, cov)
    
    if(is.null(colnames(cov))) {
      colnames(tdat1) = c("trait", paste("cov",1:dim(cov)[2],sep=""))
    } else {
      colnames(tdat1) = c("trait", colnames(cov))
    }
    
    fit1 <- glm(trait~., family = model, data=tdat1)
    yfits <- fitted.values(fit1)
    yresids <- Y - yfits
    #       fit1res1<-summary(fit1)
    #       sigma0<-sqrt(fit1res1$dispersion)
    
    Us <- XUs <- matrix(0, nrow=n, ncol=k)
    Xmus = X
    for(i in 1:k){
      tdat2 <- data.frame(X1=X[,i], cov)
      fit2 <- glm(X1~., data=tdat2)
      Xmus[,i] <- fitted.values(fit2)
      XUs[, i] <- (X[,i] - Xmus[,i])
    }
    U <- t(XUs) %*% (Y - yfits)
    
    if( model == "binomial" ) {
      CovS <- mean(yfits*(1-yfits))*(t(XUs) %*% XUs)
    } else {
      CovS <- var(yresids)*(t(XUs) %*% XUs)
    }        
    
    #        CovS<-matrix(0, nrow=k, ncol=k)
    #        for(i in 1:n)
    #            CovS<-CovS + Us[i,] %*% t(Us[i,])
    
  }
  
  return(list(U = U, CovS = CovS))
  
  
  
  
}
 
