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
#' @param quantileFun The quantile function of the test statistic,
#' either as a function or as a character string
#' @param testPargs A list of arguments passed on to quantileFun
#' @param B an integer, the number of permutations
#' @param p an integer, the number of hypotheses
#' @param warnConvergence Should a warning be thrown when the estimation
#' of the random null does not converge?
#'
#' @importFrom KernSmooth bkde
#' @importFrom stats qnorm dnorm approx quantile
#' @importFrom MASS fitdistr
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
                 maxIter, tol, estP0args, quantileFun,
                 testPargs, B, p, warnConvergence){
  if(length(statObs)!=nrow(statsPerm)){
    stop("Dimensions of observed and permutation test statistics don't match!
         ")
  }
if(length(z0Quant)==1) {z0Quant = sort(c(z0Quant, 1-z0Quant))}
  centralBorders = quantile(statObs, probs = c(z0Quant, 1-z0Quant))
  #The starting values are CRUCIAL!
  Range = range(c(statsPerm, statObs))

  #Estimate observed densities
  zValsDensObs0 = bkde(statObs, range.x = Range, gridsize = gridsize,
                       truncate = FALSE)
  zValsDensObs = zValsDensObs0$y
  zSeq = zValsDensObs0$x #The support
  zValsDensObsInterp = approx(y = zValsDensObs, x = zSeq, xout = statObs)$y

  zValsDensPerm = apply(statsPerm, 2, function(zz){
      bkde(zz, range.x = Range, gridsize = gridsize, truncate = FALSE)$y})
  zValsDensObs[zValsDensObs<=0] =
      zValsDensPerm[zValsDensPerm<=0] =
      .Machine$double.eps #Remove negative densities

  #Estimate permutation densities
  PermDensFits = apply(statsPerm, 2, estNormal)
  LogPermDensEvals = apply(PermDensFits, 2, function(fit){
    dnorm(statObs, mean = fit[1], sd = fit[2], log = TRUE)
  })
  #Indicators for the observed z values in the support of the kernel
  iter = 1L; convergence = FALSE; p0 = 1
  fdr = as.integer(statObs >= centralBorders[1] & statObs <= centralBorders[2])
  fdr[fdr==0] = .Machine$double.eps
  while(iter <= maxIter && !convergence){
    fdrOld = fdr; p0old = p0
    weights = calcWeights(logDensPerm = LogPermDensEvals,
                          fdr = fdr)
    #Null distribution
    fitAll = estNormal(y = c(statsPerm), w = rep(weights, each = p), p = p)
    g0 = dnorm(statObs, mean = fitAll[1],
               sd = fitAll[2])
    G0 = pnorm(zSeq, mean = fitAll[1],
               sd = fitAll[2])
    fdr = g0/zValsDensObsInterp*p0
    #fdr = approx(y = fdr, x = zSeq, xout = statObs)$y
    fdr[fdr>1] = 1 #Normalize here already!
    p0 = do.call(estP0,
                 c(list(statObs = statObs, nullDensCum = G0, zSeq = zSeq),
                   estP0args))
    convergence = sqrt(mean((fdrOld-fdr)^2)) < tol &&
      (p0-p0old)^2 < tol
    iter = iter + 1L
  }
  if(!convergence && warnConvergence){
      warning("Consensus null estimation did not converge,
              please investigate cause! \n")}
  return(list(PermDensFits = PermDensFits, zSeq = zSeq,
              zValsDensObs = zValsDensObs, convergence  = convergence,
              weights = weights, fdr = fdr,
              p0 = p0, iter = iter, fitAll = fitAll))
}
