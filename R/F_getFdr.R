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
getFdr = function(zValObs, G0, g0, fdr, zSeq, p, p0,...)
{
  g0 = g0/(sum(g0))
  #Null
  G0cumDef = sapply(seq_along(G0), function(yy){
    if(G0[yy]<0.5) {G0[yy]} else {1-G0[yy]+g0[yy]}
    })

  #Observed
  zcum = ecdf(zValObs)(zValObs)
  ZcumDef = sapply(zcum, function(yy){if(yy<0.5) yy else (1-yy+1/p)})
  id = sapply(zValObs, function(rr){which.max(rr<=zSeq)})

  #Fdr
  Fdr = G0cumDef[id]/ZcumDef*p0

  Fdr[Fdr>1] = 1
return(list(Fdr = Fdr, fdr = fdr[id], p0 = min(p0,1)))
}
