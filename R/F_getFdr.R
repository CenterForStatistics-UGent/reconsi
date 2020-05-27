#' Calculate tail-area (Fdr) and local (fdr) false discovery rates,
#'    based on a certain null distribution
#'
#' @param statObs Vector of observed z-values
#' @param fitAll The parameters of the estimated random null
#' @param fdr local false discovery rate, already estimated
#' @param zSeq Support of the density estimation
#' @param p the number of hypotheses
#' @param p0 The estimated fraction of null hypotheses
#' @param zValsDensObs estimated densities of observed test statistics
#' @param smoothObs A boolean, should estimated observed densities of the test
#' statistics be used in estimating the Fdr
#' @param assumeNormal A boolean, should normality be assumed for the null distribution?
#' @param ... more arguments, ignored
#'
#' @importFrom stats ecdf approx
#'
#' @return A list with components
#' \item{Fdr}{Tail are false discovery rate}
#' \item{fdr}{Local false discovery rate}
getFdr = function(statObs, fitAll, fdr, zSeq, p, p0, zValsDensObs, smoothObs,
                  assumeNormal, ...)
{
  statObsNotNA = statObs[!is.na(statObs)]
  #Null
  if(assumeNormal){
      G0 = pnorm(statObsNotNA, mean = fitAll["mean"], sd = fitAll["sd"])
      G0[G0>0.5] = pnorm(statObsNotNA[G0>0.5], mean = fitAll["mean"],
                         sd = fitAll["sd"], lower.tail = FALSE)
  } else {
      G0 = approx(y = cumsum(fitAll/sum(fitAll)), xout = statObsNotNA, x = zSeq)$y
      G0[G0>0.5] = 1-G0[G0>0.5]
  }
  #Observed
  zcum = if(smoothObs) {
    approx(y = cumsum(zValsDensObs/sum(zValsDensObs)),
    xout = statObsNotNA, x = zSeq)$y
  } else {
       ecdf(statObsNotNA)(statObsNotNA)
      }
  zcum[zcum>0.5] = 1-zcum[zcum>0.5]
  #Fdr
  Fdr = G0/zcum*p0
  Fdr[Fdr>1] = 1
  #fdr
  if(is.null(fdr)){
      g0 = if(assumeNormal){
          pnorm(statObsNotNA, mean = fitAll["mean"], sd = fitAll["sd"])
      } else {
          fitAll
      }
    fdr  = g0/approx(y = zValsDensObs, xout = statObsNotNA, x = zSeq)$y*p0
    fdr[fdr>1] = 1
  }
  FdrOut = fdrOut = statObs
  FdrOut[!is.na(statObs)] = Fdr
  fdrOut[!is.na(statObs)] = fdr
return(list(Fdr = FdrOut, fdr = fdrOut))
}
