#'Estimate the fraction of true null hypotheses.
#'
#' @param zValObs A vector of observed z-values
#' @param nullDensCum A smoothed cumulative distribution function
#'  of the observed z-values
#' @param zSeq The support of the kernel smoother
#' @param z0quantRange a number of quantiles between 0 and 0.5
#' @param smooth.df degrees of freedom for the spline smoother
#'
#' @return The estimated null fraction, the value of the spline evaluated
#'  at the first element of z0quantRange
#'
#' @details A natural spline is used over a range of intervals.
#' Based on the qvalue::qvalue() function and
#'Storey and Tibshirani, 2003
#'
#'@importFrom stats smooth.spline predict
estP0 = function(zValObs, nullDensCum, zSeq, z0quantRange, smooth.df){
  pi0 = sapply(z0quantRange, function(z0Quant) {
    z0Quant = c(z0Quant, 1-z0Quant)
    centralBorders = zSeq[c(which.max(nullDensCum > z0Quant[1]),
                            which.max(nullDensCum > z0Quant[2]))]
    mean(zValObs <= centralBorders[2] & z
         ValObs >= centralBorders[1])/diff(z0Quant)
  })
  smoothPi0 = smooth.spline(z0quantRange, pi0, df = smooth.df)
  pi0Smooth = predict(smoothPi0, x = z0quantRange)$y
  min(pi0Smooth[1], 1)
}
