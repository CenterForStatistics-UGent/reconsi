#' Fast estimation of mean and standard deviation of a normal distrbution,
#' optionally with weights
#' @param y vector of observations
#' @param w optional weight vector
#' @param p The number of features
#'
#' @return A vector of length 2 with mean and standard deviation
#' @importFrom stats sd weighted.mean
estNormal = function(y, w = NULL, p){
    if(is.null(w)){
        # fit = lm.fit(y = y, x = x)
        # c(mean = fit$coef, sd = sqrt(mean(fit$residuals^2)))
        c(mean = mean(y), sd = sd(y))
        #ML estimate = biased!
    } else {
        # fit = lm.wfit(y = c(y), x = as.matrix(rep.int(1L, p*B)), w = w)
        # c(mean = fit$coef, sd = sqrt(sum(fit$residuals^2*w)/p))
        wmean = weighted.mean(y, w = w)
        c(mean = wmean, sd = sqrt(sum((y-wmean)^2*w)*1/(p-1)))
    }
}