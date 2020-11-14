#' Plot an the corvariance matrix of the test statistics estimated through permutations
#' @examples
#' p = 200; n = 50; B = 5e1
#' x = rep(c(0,1), each = n/2)
#' mat = cbind(
#' matrix(rnorm(n*p/10, mean = 5+x),n,p/10), #DA
#' matrix(rnorm(n*p*9/10, mean = 5),n,p*9/10) #Non DA
#' )
#' mat = mat = mat + rnorm(n, sd = 0.3) #Introduce some dependence
#' fdrRes = reconsi(mat, x, B = B)
#' plotCovar(fdrRes)
#' @export
#' @param reconsiFit The reconsi fit
#' @param col,xlab,ylab,... A list of arguments for the image() function.
#' @importFrom graphics image
#' @importFrom grDevices colorRampPalette
#' @details By default, yellow indicates negative correlaton between test statistics,
#' blue positive correlation
#' @note Note the difference with the plotApproxCovar function, where the
#' covariances between binned test statistics are shown to get an idea between
#' covariances between tail and center values of the univariate null distribution.
#' Here the covariance matrix between all test statistics is shown
#' @return invisible()
#' @seealso \link{plotApproxCovar}
plotCovar = function(reconsiFit, col = colorRampPalette(
    c("yellow","blue"))(12),
    xlab = "Test statistic index", ylab = xlab, ...){
    Seq = seq_len(nrow(reconsiFit$statsPerm))
    image(z = var(t(reconsiFit$statsPerm))[,nrow(reconsiFit$statsPerm):1],
          xlab = xlab, ylab = ylab, col = col, x = Seq, y = Seq, ...)
}
