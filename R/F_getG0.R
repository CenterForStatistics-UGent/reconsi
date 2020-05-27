#' Obtain the consensus null
#'
#' @param statObs A vector of lenght p with observed test statistics
#' @param statsPerm  A pxB matrix with permuation z-values
#' @param z0Quant a vector of length of quantiles defining the central part R
#'    of the distribution. If a single number is supplied, then
#'    (z0quant, 1-z0quant) will be used
#' @param gridsize An integer, the gridsize for the density estimation
#' @param maxIter An integer, the maximum number of iterations in determining R
#' @param tol The convergence tolerance.
#' @param estP0args A list of arguments passed on to the estP0args() function
#' @param testPargs A list of arguments passed on to quantileFun
#' @param B an integer, the number of permutations
#' @param p an integer, the number of hypotheses
#' @param pi0 A known fraction of true null hypotheses
#' @param assumeNormal A boolean, should normality be assumed for the null distribution?
#'
#' @importFrom KernSmooth bkde
#' @importFrom stats qnorm dnorm approx quantile
#'
#' @return A list with following entries
#' \item{PermDensFits}{The permutation density fits}
#' \item{zSeq}{The support of the kernel for density estimation}
#' \item{zValsDensObs}{The estimated densities of the observed z-values}
#' \item{convergence}{A boolean, has the algorithm converged?}
#' \item{weights}{Vector of length B with weights
#'    for the permutation distributions}
#' \item{fdr}{Estimated local false discovery rate along the support
#'    of the kernel}
#' \item{p0}{The estimated fraction of true null hypotheses}
#' \item{iter}{The number of iterations}
#' \item{fitAll}{The consensus null fit}
getG0 = function(statObs, statsPerm, z0Quant, gridsize,
                 maxIter, tol, estP0args, testPargs, B, p,
                 pi0, assumeNormal){
  if(length(statObs)!=nrow(statsPerm)){
    stop("Dimensions of observed and permutation test statistics don't match!")
  }
if(length(z0Quant)==1) {
    z0Quant = sort(c(z0Quant, 1-z0Quant))
}
 estPi0 = is.null(pi0) #Should the fraction of nulls be estimated?
  statObs = statObs[!is.na(statObs)] #ignore NA values
  centralBorders = quantile(statObs, probs = z0Quant)
  #Estimate observed densities
  zValsDensObs0 = bkde(statObs, gridsize = gridsize)
  zValsDensObs = zValsDensObs0$y
  zSeq = zValsDensObs0$x #The support
  zValsDensObsInterp = approx(y = zValsDensObs, x = zSeq, xout = statObs)$y
  zValsDensObs[zValsDensObs<=0] =
      .Machine$double.eps #Remove negative densities
  #Estimate permutation densities
  PermDensFits = apply(statsPerm, 2, estNormal)
  LogPermDensEvals = apply(PermDensFits, 2, function(fit){
    dnorm(statObs, mean = fit[1], sd = fit[2], log = TRUE)
  })
  #Indicators for the observed z values in the support of the kernel
  iter = 1L; convergence = FALSE; p0 = 1; fitAll = c("mean" = 0, "sd" = 1)
  fdr = as.integer(statObs >= centralBorders[1] & statObs <= centralBorders[2])
  fdr[fdr==0] = .Machine$double.eps
  while(iter <= maxIter && !convergence){
      fdrOld = fdr; p0old = p0
      weights = calcWeights(logDensPerm = LogPermDensEvals, fdr = fdr)
      #Null distribution
      if(assumeNormal){
        fitAll = estNormal(y = c(statsPerm), w = rep(weights, each = p), p = p)
        g0Obs = dnorm(statObs, mean = fitAll[1], sd = fitAll[2])
      } else {
        g0Obs = wkde(x = c(t(statsPerm)), w = weights/p, u = statObs)
      }
    fdr = g0Obs/zValsDensObsInterp*p0
    fdr[fdr>1] = 1 #Normalize here already!
    p0 = if(estPi0) do.call(estP0, c(list(statObs = statObs, fitAll = fitAll),
                   estP0args)) else pi0
    convergence = all((fdr-fdrOld)^2 < tol) && (p0-p0old)^2 < tol
    iter = iter + 1L
  }
  if(!convergence){
      warning("Consensus null estimation did not converge, please investigate cause! \n")}
  return(list(PermDensFits = PermDensFits, zSeq = zSeq,
              zValsDensObs = zValsDensObs, convergence  = convergence,
              Weights = weights, fdr = fdr,
              p0 = p0, iter = iter, assumeNormal = assumeNormal,
              fitAll = if(assumeNormal) fitAll else wkde(x = c(t(statsPerm)), w = weights/p, u = zSeq)))
}
