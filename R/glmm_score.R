#' GLMM Score
#'
#' It returns the score vector and its corresponding covariance matrix.
#'
#' @param fit0 A generalized linear mixed model fitted under the null hypothesis using the ``glmmkin" function of the GMMAT R package, where p is the number of variants.
#'
#' @param G A n * p genotype matrix. n is the number of samples and the order must correspond with fit0.
#'
#' @export
#'
#' @return score vector (U) and corresponding covariance matrix (V).
#'
#' @author Chong Wu, Jun Young Park, and Wei Pan
#'

glmm_score<-function(fit0, G){
    Y<-fit0$Y
    P<-fit0$P
    
    U = t(G)%*%P%*%Y
    V = t(G)%*%P%*%G
    
    return(list(U=U, V=V))
}

#fit0 :
#G : a
