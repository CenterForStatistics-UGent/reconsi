#' Calculate tail-area (Fdr) and local (fdr) false discovery rates,
#'    based on a certain null distribution
#'
#' @param zValObs Vector of observed z-values
#' @param G0 null distribution function
#' @param g0 Null density function
#' @param fdr local false discovery rate, already estimated
#' @param zSeq Support of the density estimation
#' @param p the number of hypotheses
#' @param p0 The estimated fraction of null hypotheses
#' @param ... more arguments, ignored
#'
#' @importFrom stats ecdf
#'
#' @return A list with components
#' \item{Fdr}{Tail are false discovery rate}
#' \item{fdr}{Local false discovery rate}
#' \item{p0}{The proportion of true null hypotheses}
getFdr = function(zValObs, G0, g0, fdr, zSeq, p, p0, zValsDensObs,
                  smoothObs, ...)
{
  #Null
  G0[G0>0.5] = 1-G0[G0>0.5]
  G0cumDefInterp = approx(x = zSeq, y = G0, xout = zValObs)$y

  #Observed
  zcum = if(smoothObs) {
    approx(y = cumsum(zValsDensObs/sum(zValsDensObs)),
    xout = zValObs, x = zSeq)$y
   } else {ecdf(zValObs)(zValObs)}
  zcum[zcum>0.5] = 1-zcum[zcum>0.5]+1/p


  #Fdr
  Fdr = G0cumDefInterp/zcum*p0

  Fdr[Fdr>1] = 1
return(list(Fdr = Fdr, fdr = fdr, p0 = p0))
}
