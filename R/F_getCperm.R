#' A function to obtain a null covariance matrix of binned test statistics
#' @param zValsMat The pxB matrix of permutation z-values in the columns
#' @param nBins an integer, the number of bins
#' @param binEdges A vector of length 2 with the outer bin edges
#'
#' @return The covariance matrix of binned z-values
#'
#' @note This does not calculate the covariance matrix of the p test statistics, but only of the nBins bin counts, for illustrative purposes.
getCperm = function(zValsMat, nBins = 82L, binEdges = c(-4.1,4.1)){
  Breaks = c(-Inf, seq(binEdges[1], binEdges[2], length.out = nBins+1), Inf)
  bootYs = apply(zValsMat, 2, function(xx){table(cut(xx, breaks = Breaks))})
  avCounts = rowMeans(bootYs)
  tcrossprod(bootYs-avCounts)/(ncol(zValsMat)-1) #The permutation matrix
}