#' Fast estimation of mean and standard deviation of a normal disitrbution, optionally with weights
#' @param y vector of observations
#' @param x design matrix
#' @param w optional weight vector
#' @param p an integer, number of observations
#' @param B an integer, the number of bootstraps
#'
#' @return A vector of length 2 with mean and standard deviation
#' @importFrom stats lm.fit lm.wfit
estNormal = function(y, w = NULL, pbFrac = NULL){
    if(is.null(w)){
        # fit = lm.fit(y = y, x = x)
        # c(mean = fit$coef, sd = sqrt(mean(fit$residuals^2)))
        c(mean = mean(y), sd = sd(y))
    } else {
        # fit = lm.wfit(y = c(y), x = as.matrix(rep.int(1L, p*B)), w = w)
        # c(mean = fit$coef, sd = sqrt(sum(fit$residuals^2*w)/p))
        wmean = weighted.mean(y, w = w)
        c(wmean, sd = sqrt(sum((y-wmean)^2*w)*pbFrac))
    }
}
