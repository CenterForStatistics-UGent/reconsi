#' Obtain a null covariance matrix of binned test statistics
#'
#' @note This is not the covariance matrix of the p test statistic, nor of the
#'    data! It is an approximate covariance matrix of binned test statistics for
#'    visualization purposes.
#'
#' @param zValsMat The pxB matrix of permutation z-values in the columns
#' @param nBins an integer, the number of bins
#' @param binEdges A vector of length 2 with the outer bin edges
#'
#' @return The covariance matrix of binned z-values
#'
getCperm = function(zValsMat, nBins = 82L, binEdges = c(-4.1,4.1)){
  Breaks = c(-Inf, seq(binEdges[1], binEdges[2], length.out = nBins+1), Inf)
  bootYs = apply(zValsMat, 2, function(xx){table(cut(xx, breaks = Breaks))})
  tcrossprod(bootYs-rowMeans(bootYs))/(ncol(zValsMat)-1) #The permutation matrix
}
