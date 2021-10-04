#' Find the dependence pat C1 of the approximate covariance matrix, and extract the ratio of the first eigenvalue to the sum of all positive eigenvalues
#'
#' @param statsPerm Matrix of permuted test statistics
#' @param numEig An integer, number of first eigenvalues
#' @param ... passed onto binStats
#' @return A proportion indicating the ratio of the first eigenvalues to the sum of all eigenvalues
#' @export
#' @importFrom Matrix nearPD
getC1prop = function(statsPerm, numEig = 1, ...){
    bootYs = binStats(statsPerm, ...)
    Ctot = tcrossprod(bootYs-rowMeans(bootYs))/(ncol(statsPerm)-1)
    nu = rowMeans(bootYs)
    C0 = diag(nu)-tcrossprod(nu)/length(statsPerm)
    C1 = Ctot-C0
    eigenValues = nearPD(C1)$eigenvalues #Convert to nearest positive definite matrix, this is sampling error
    sum(eigenValues[seq_len(numEig)])/sum(eigenValues)##The crucial property. this could be used to check the approximation
}
