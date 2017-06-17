#' Doubley Adaptive Sum of Powered Score tests (daSPU) test for single trait with multiple weights for each genetic marker
#'
#' It returns p-values and test statistics
#'
#' @param U Score vector for the marker set we are interested in. (N by 1 matrix)
#'
#' @param V Corresponding covariance matrix for the score vector. (N by N matrix)
#'
#' @param weight Multiple weights for each genetic markers. (N by M matrix)
#'
#' @param pow1 SNP or single CpG sites specific power(gamma values) used in daSPU test.
#'
#' @param pow2 Specific power(gamma values) used for different weight.
#'
#' @param n.perm number of permutations.
#'
#' @export
#' @return P-values for daSPU test.
#'
#' @author Chong Wu and Wei Pan


daSPU <- function(U,V,weight, pow1 = c(1:8,Inf),pow2 = c(1:8,Inf),n.perm = 1000) {
    weight = as.matrix(weight)
    if(n.perm < 500000) {
        sim.aSPUw2(U,V,weight,pow1,pow2,n.perm)
    } else {
        sim.aSPUw4(U,V,weight,pow1,pow2,n.perm)
    }
}

