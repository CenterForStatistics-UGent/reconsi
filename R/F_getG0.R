#' Obtain the consensus null
#'
#' @param zValObs A vector of lenght p with observed test statistics
#' @param zValsMat  A pxB matrix with permuation z-values
#' @param z0quant a vector of length of quantiles defining the central part R
#'    of the distribution. If a single number is supplied, then
#'    (z0quant, 1-z0quant) will be used
#' @param weightStrat A character vector, indicating the weighting strategy
#' @param gridsize An integer, the gridsize for the density estimation
#' @param maxIter An integer, the maximum number of iterations in determining R
#' @param tol The convergence tolerance.
#' @param estP0args A list of arguments passed on to the estP0args() function
#' @param ... further arguments passed on to the KernSmooth::bkde() function
#'
#' @importFrom KernSmooth bkde
#'
#' @return A list with following entries
#' \item{g0}{The consensus density function}
#' \item{G0}{The consensus distribution function}
#' \item{R}{The central region}
#' \item{zSeq}{The support of the kernel for density estimation}
#' \item{zValsDensObs}{The estimated densities of the observed z-values}
#' \item{zValsDensPerm}{The estimated densities of the permutation z-values}
#' \item{convergence}{A boolean, has the algorithm converged?}
#' \item{weights}{Vector of length B with weights
#'    for the permutation distributions}
#' \item{fdr}{Estimated local false discovery rate along the support
#'    of the kernel}
#' \item{p0}{The estimated fraction of true null hypotheses}
#' \item{iter}{The number of iterations}
getG0 = function(zValObs, zValsMat, z0Quant, weightStrat, gridsize,
                 maxIter, tol, estP0args,...){
if(length(z0Quant)==1) {z0Quant = sort(c(z0Quant, 1-z0Quant))}
  centralBorders = qnorm(z0Quant)
  Range = range(c(zValsMat, zValObs))
  zValsDensPerm = apply(zValsMat, 2, function(zz){
      bkde(zz, range.x = Range, gridsize = gridsize, truncate = FALSE,...)$y})
  zValsDensObs0 = bkde(zValObs, range.x = Range, gridsize = gridsize,
                       truncate = FALSE,...)
  zValsDensObs = zValsDensObs0$y
  zSeq = zValsDensObs0$x #The support
  zValsDensObs[zValsDensObs<=0] =
      zValsDensPerm[zValsDensPerm<=0] =
      .Machine$double.eps #Remove negative densities
  zIndObs = sapply(zValObs, function(rr){which.max(rr<=zSeq)})
  #Indicators for the observed z values in the support of the kernel
  iter = 1L; convergence  = FALSE; p0 = 1
  fdr = as.integer(zSeq >= centralBorders[1] & zSeq <= centralBorders[2]) +
      .Machine$double.eps
  while(iter <= maxIter && !convergence){
    fdrOld = fdr[zIndObs]
    weights = calcWeights(densPerm = zValsDensPerm, zIndObs = zIndObs,
                          weightStrat = weightStrat, fdr = fdr,
                          zIndR = which(zValObs > centralBorders[1] &
                                            zValObs < centralBorders[2]))
    #Null distribution
    g0 = rowSums(rowMultiply(zValsDensPerm, weights))
    G0 = cumsum(g0/sum(g0))
    fdr = g0/zValsDensObs*p0
    fdr[fdr>1] = 1
    centralBorders = zSeq[c(which.max(G0 > z0Quant[1]),
                            which.max(G0 > z0Quant[2]))]
    p0 = do.call(estP0,
                 c(list(zValObs = zValObs, nullDensCum = G0, zSeq = zSeq),
                   estP0args))
    convergence = sqrt(mean((fdrOld-fdr[zIndObs])^2)) < tol
    iter = iter + 1L
  }
  if(!convergence){
      warning("Consensus null estimation did not converge,
              please investigate cause! \n")}
  return(list(g0 = g0, G0 = G0, R = centralBorders, zSeq = zSeq,
              zValsDensObs = zValsDensObs, zValsDensPerm = zValsDensPerm,
              convergence  = convergence, weights = weights, fdr = fdr,
              p0 = p0, iter = iter))
}
