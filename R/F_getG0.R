#' Obtain the consensus null
#'
#' @param statObs A vector of lenght p with observed test statistics
#' @param statsPerm  A pxB matrix with permuation z-values
#' @param z0Quant a vector of length of quantiles defining the central part R
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
#' @importFrom stats qnorm
#' @importFrom MASS fitdistr
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
getG0 = function(statObs, statsPerm, z0Quant, weightStrat, gridsize,
                 maxIter, tol, estP0args, normAsump, normAsumpG0, quantileFun,
                 testPargs, B, p,...){
  if(length(statObs)!=nrow(statsPerm)){
    stop("Dimensions of observed and permutation test statistics don't match!\n")
  }
if(length(z0Quant)==1) {z0Quant = sort(c(z0Quant, 1-z0Quant))}
  centralBorders = do.call(quantileFun, c(list(p = z0Quant), testPargs))
  Range = range(c(statsPerm, statObs))

  #Estimate observed densities
  zValsDensObs0 = bkde(statObs, range.x = Range, gridsize = gridsize,
                       truncate = FALSE,...)
  zValsDensObs = zValsDensObs0$y
  zSeq = zValsDensObs0$x #The support

  #Estimate permuation densities
  x1 = as.matrix(rep.int(1L, p))
  xall = as.matrix(rep.int(1L, p*B))
  zValsDensPerm = if(normAsump){
    apply(statsPerm, 2, function(zz){
      fit = estNormal(y = zz, x =x1, p = p, B=B)
      dnorm(zSeq, mean = fit[1], sd = fit[2])
    })
  } else {apply(statsPerm, 2, function(zz){
    bkde(zz, range.x = Range, gridsize = gridsize, truncate = FALSE,...)$y})}
  zValsDensObs[zValsDensObs<=0] =
      zValsDensPerm[zValsDensPerm<=0] =
      .Machine$double.eps #Remove negative densities

  #Interpolate estimated densities
  #obsDensInterp = approx(xout = statObs, x = zSeq, y = zValsDensObs)$y
  LogPermDensInterp = log(apply(zValsDensPerm, 2, function(dens){
    approx(x = zSeq, y = dens, xout = statObs)$y
  }))

  #zIndObs = sapply(statObs, function(rr){which.min(abs(rr-zSeq))})
  #Indicators for the observed z values in the support of the kernel
  iter = 1L; convergence  = FALSE; p0 = 1
  fdr = as.integer(zSeq >= centralBorders[1] & zSeq <= centralBorders[2]) +
      .Machine$double.eps
  while(iter <= maxIter && !convergence){
    fdrOld = fdr
    weights = calcWeights(logDensPerm = LogPermDensInterp,
                          weightStrat = weightStrat,
                          fdr = if(weightStrat =="LHw") {approx(y = fdr, x = zSeq,
                                       xout = statObs)$y} else {NULL},
                          zIndR = if(weightStrat =="LHw") {NULL
                            } else {which(statObs > centralBorders[1] &
                                            statObs < centralBorders[2])})
    #Null distribution
    if(normAsumpG0){
    fitAll = estNormal(y = c(statsPerm), x = xall,
                    w = rep(weights, each = p), p = p, B=B)
    g0 = dnorm(zSeq, mean = fitAll[1],
               sd = fitAll[2])    #Multiply by total sum of weights
    } else {
      g0 =rowSums(rowMultiply(zValsDensPerm, weights))}
    G0 = cumsum(g0/sum(g0))
    fdr = g0/zValsDensObs*p0
    centralBorders = zSeq[c(which.max(G0 > z0Quant[1]),
                            which.max(G0 > z0Quant[2]))]
    p0 = do.call(estP0,
                 c(list(statObs = statObs, nullDensCum = G0, zSeq = zSeq),
                   estP0args))
    convergence = sqrt(mean((fdrOld-fdr)^2)) < tol
    #    cat(iter, "\n")
    iter = iter + 1L
  }
  fdr = approx(y = fdr, x= zSeq, xout = statObs)$y
  fdr[fdr>1] = 1 #Only normalize here?
  if(!convergence){
      warning("Consensus null estimation did not converge,
              please investigate cause! \n")}
  return(list(g0 = g0, G0 = G0, R = centralBorders, zSeq = zSeq,
              zValsDensObs = zValsDensObs, zValsDensPerm = zValsDensPerm,
              convergence  = convergence, weights = weights, fdr = fdr,
              p0 = p0, iter = iter))
}

estNormal = function(y, x, w  = NULL, p, B){
  if(is.null(w)){
  fit = lm.fit(y = y, x = x)
  c(mean = fit$coef, sd = sqrt(sum(fit$residuals^2)/p))
  } else {
  fit = lm.wfit(y = c(y), x = as.matrix(rep.int(1L, p*B)), w = w)
  c(mean = fit$coef, sd = sqrt(sum(fit$residuals^2*w)/p))
  }

}
