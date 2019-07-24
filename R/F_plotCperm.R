#' Plot an approximatio of the correlation structure of the test statistics
#' @examples
#' p = 200; n = 50; B = 5e1
#' x = rep(c(0,1), each = n/2)
#' mat = cbind(
#' matrix(rnorm(n*p/10, mean = 5+x),n,p/10), #DA
#' matrix(rnorm(n*p*9/10, mean = 5),n,p*9/10) #Non DA
#' )
#' mat = mat = mat + rnorm(n, sd = 0.3) #Introduce some dependence
#' fdrRes = rransi(mat, x, B = B)
#' plotCperm(fdrRes$statsPerm)
#' @export
#' @param statsPerm The matrix of resample test statistics
#' @param imageList A list of arguments for the image() function.
#' @param ... passed on to the getCperm function
#' @importFrom graphics image
#' @importFrom grDevices colorRampPalette
#' @details By default, yellow indicates negative correlaton between bin counts,
#' blue positive correlation
#' @note This is not the covariance matrix of the p test statistic, nor of the
#'    data! It is an approximate covariance matrix of binned test statistics for
#'    visualization purposes.
#' @return invisible()
plotCperm = function(statsPerm,
                     imageList = list(
                         col = colorRampPalette(
                         c("yellow","blue"))(12),
                         x = seq(-4.2, 4.2, 0.1),
                         y = seq(-4.2, 4.2, 0.1),
                         xlab = "Z-values", ylab = "Z-values"), ...){
    do.call(what = image, args = c(list(z = getCperm(statsPerm, ...)),
                                   imageList))
}
