#'Estimate the fraction of true null hypotheses.
#'
#' @param statObs A vector of observed z-values
#' @param fitAll the estimated normal null
#' @param z0quantRange a number of quantiles between 0 and 0.5
#' @param smooth.df degrees of freedom for the spline smoother
#' @param evalVal the value of q at which to evaluate the spline
#' @param assumeNormal A boolean, should normality be assumed for the null distribution?
#'
#' @return The estimated null fraction, the value of the spline evaluated
#'  at the first element of z0quantRange
#'
#' @details A natural spline is used over a range of intervals.
#' Based on the qvalue::qvalue() function and
#'Storey and Tibshirani, 2003
#'
#'@importFrom stats smooth.spline predict
#'@importFrom ks pkde
estP0 = function(statObs, fitAll, z0quantRange, smooth.df, evalVal, assumeNormal){
    pi0 = vapply(z0quantRange, FUN.VALUE = double(1), function(z0Quant) {
    z0Quant = c(z0Quant, 1-z0Quant)
    if(assumeNormal) {
        centralBorders = qnorm(z0Quant, mean = fitAll["mean"], sd = fitAll["sd"])
    } else {
        cumul.prob <- pkde(q = fitAll$eval.points, fhat = fitAll)
        Ord = order(cumul.prob)
        ind <- findInterval(x = z0Quant, vec = cumul.prob[Ord])
        centralBorders = fitAll$eval.points[Ord][ind]
    }
    mean(statObs <= centralBorders[2] &
        statObs >= centralBorders[1])/diff(z0Quant)
    })
    smoothPi0 = smooth.spline(z0quantRange, pi0, df = smooth.df)
    pi0Smooth = predict(smoothPi0, x = evalVal)$y #Conservative
    min(pi0Smooth, 1)
}
