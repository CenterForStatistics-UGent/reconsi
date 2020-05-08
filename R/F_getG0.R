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
#' @param pi0 A known fraction of true null hypotheses.
#'
#' @importFrom stats dnorm approx
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
#' \item{G0Z}{The consensus null fit, evaluated at zSeq}
#' \item{PermDensEvals}{The densities of the resampling distribution, evaluated at zSeq}
getG0 = function(statObs, statsPerm, z0Quant, gridsize,
                 maxIter, tol, estP0args, testPargs, B, p, warnConvergence,
                 pi0){
  if(length(statObs)!=nrow(statsPerm)){
    stop("Dimensions of observed and permutation test statistics don't match!")
  }
if(length(z0Quant)==1) {
    z0Quant = sort(c(z0Quant, 1-z0Quant))
}
 estPi0 = is.null(pi0) #Should the fraction of nulls be estimated?
  statObs = statObs[!is.na(statObs)] #ignore NA values
  #The starting values are CRUCIAL!
  Range = range(c(statsPerm, statObs), na.rm = TRUE)
  #Estimate observed densities
  zValsDensObs0 = bkdeStab(statObs, range.x = Range, gridsize = gridsize,
                       truncate = FALSE)
  zValsDensObsInterp = interpKDE(zValsDensObs0, statObs)
  #Estimate permutation densities
  PermDensFits = apply(statsPerm, 2, bkdeStab, range.x = Range, gridsize = gridsize,
                       truncate = FALSE)
  PermDensEvals = vapply(PermDensFits, FUN.VALUE = statObs, function(fit){
    interpKDE(fit, newData = statObs)
  })
  PermDensEvalsZ = vapply(PermDensFits, FUN.VALUE = numeric(gridsize), function(fit){
      fit$y
  })
  LogPermDensEvals = log(PermDensEvals)
  #Starting values
  iter = 1L; convergence = FALSE; p0 = 1; weights = rep(1/B, B)
  while(iter <= maxIter && !convergence){
     p0old = p0; weightsOld = weights
    g0 = rowSums(rowMultiply(PermDensEvals, weights))
    g0Z = rowSums(rowMultiply(PermDensEvalsZ, weights))
    fdr = g0/zValsDensObsInterp*p0
    fdr[fdr>1] = 1 #Normalize here already!
    weights = calcWeights(logDensPerm = LogPermDensEvals, fdr = fdr)
    p0 = if(estPi0) do.call(estP0, c(list(statObs = statObs, zSeq = zValsDensObs0$x,
                                          nullDensCum = cumsum(g0Z)/sum(g0Z)),
                   estP0args)) else pi0
    convergence = (p0-p0old)^2 < tol && sqrt(sum((weightsOld - weights)^2)) < tol
    iter = iter + 1L
  }
  if(!convergence){
      warning("Consensus null estimation did not converge, please investigate cause! \n")}
  return(list(zSeq = zValsDensObs0$x,
              zValsDensObs = zValsDensObs0$y, convergence  = convergence,
              weights = weights, fdr = fdr, PermDensEvals = PermDensEvalsZ,
              p0 = p0, iter = iter, g0Z = g0Z))
}
