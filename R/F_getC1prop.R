#' Find the dependence pat C1 of the approximate covariance matrix, and extract the ratio of the first eigenvalue to the sum of all positive eigenvalues
#'
#' @param statsPerm Matrix of permuted test statistics
#' @param numEig An integer, number of first eigenvalues
#' @param ... passed onto binStats
#' @return A proportion indicating the ratio of the first eigenvalues to the sum of all eigenvalues
#' @export
#' @importFrom Matrix nearPD
#' @examples
#' p = 200; n = 50; B = 5e1
#' x = rep(c(0,1), each = n/2)
#' mat = cbind(
#' matrix(rnorm(n*p/10, mean = 5+x),n,p/10), #DA
#' matrix(rnorm(n*p*9/10, mean = 5),n,p*9/10) #Non DA
#' )
#' mat = mat = mat + rnorm(n, sd = 0.3) #Introduce some dependence
#' fdrRes = reconsi(mat, x, B = B)
#' getC1prop(fdrRes$statsPerm)
getC1prop = function(statsPerm, numEig = 1, ...){
    bootYs = binStats(statsPerm, ...)
    Ctot = tcrossprod(bootYs-rowMeans(bootYs))/(ncol(statsPerm)-1)
    nu = rowMeans(bootYs)
    C0 = diag(nu)-tcrossprod(nu)/length(statsPerm)
    C1 = Ctot-C0
    eigenValues = nearPD(C1)$eigenvalues #Convert to nearest positive definite matrix, this is sampling error
    sum(eigenValues[seq_len(numEig)])/sum(eigenValues)##The crucial property. this could be used to check the approximation
}
