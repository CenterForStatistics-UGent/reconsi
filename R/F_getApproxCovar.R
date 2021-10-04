#' Obtain a null covariance matrix of binned test statistics
#'
#' @note This is not the covariance matrix of the p test statistic, nor of the
#'    data! It is an approximate covariance matrix of binned test statistics for
#'    visualization and diagnostic purposes.
#'
#' @param statsPerm The pxB matrix of permutation z-values in the columns
#' @param ... passed on to binStats
#'
#' @return The covariance matrix of binned z-values
#' @export
#' @examples
#' p = 200; n = 50; B = 5e1
#' x = rep(c(0,1), each = n/2)
#' mat = cbind(
#' matrix(rnorm(n*p/10, mean = 5+x),n,p/10), #DA
#' matrix(rnorm(n*p*9/10, mean = 5),n,p*9/10) #Non DA
#' )
#' mat = mat = mat + rnorm(n, sd = 0.3) #Introduce some dependence
#' fdrRes = reconsi(mat, x, B = B)
#' corMat = getApproxCovar(fdrRes$statsPerm)
getApproxCovar = function(statsPerm, ...){
    bootYs = binStats(statsPerm, ...)
  tcrossprod(bootYs-rowMeans(bootYs))/(ncol(statsPerm)-1)
}
#' Bin the test statistic into equally sized bins
#' @param z the matrix of permuted test statistics
#' @param nBins an integer, the number of bins
#' @param binEdges A vector of length 2 with the outer bin edges
#' @return Matrix of binned test statistics
binStats = function(z, nBins = 82L, binEdges = c(-4.1,4.1)){
    Breaks = c(-Inf, seq(binEdges[1], binEdges[2], length.out = nBins+1), Inf)
    binnedZ = apply(z, 2, function(xx){table(cut(xx, breaks = Breaks))})
    return(binnedZ)
}
