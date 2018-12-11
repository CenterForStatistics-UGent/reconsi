#' Perform simultaneous inference, correcting for correlation between tests.
#' @param Y the matrix of sequencing counts
#' @param x a grouping factor
#' @param B the number of permutations
#' @param test Character string, giving the name of the function
#'  to test for differential absolute abundance.
#'  Must accept the formula interface
#' @param argList A list of arguments, passed on to the testing function
#' @param quantileFun the name of the distribution function of the test statistic
#' @param zValues A boolean, should test statistics be converted to z-values.
#' See details
#' @param testPargs A list of arguments passed on to quantileFun
#' @param z0Quant A vector of length 2 of quantiles of the null distribution,
#'    in between which only null values are expected
#' @param gridsize The number of bins for the kernel density estimates
#' @param weightStrat A character vector, indicating the weighting strategy.
#'    Either "LH" for likelihoods based on the central region,
#'    or "LHw" for weighted likelihoods
#' @param maxIter An integer, the maximum number of iterations in the estimation
#'    of the null distribution
#' @param tol The tolerance for the infinity norm of the central borders
#'    in the iterative procedure
#' @param cPerm A boolean, should the covariance matrix
#'    of the binned counts be returned?
#' @param nBins The number of bins for the binning of the test statistic
#'    in the calculation of the correlation matrix
#' @param binEdges The edges of the furthest bins
#' @param center A boolean, should observations be centered
#'    in each group prior to permuations? See details.
#' @param zVals An optional list of observed (zValObs) and
#' permutation (zValsPerm) z-values. If supplied, the calculation
#'    of the observed and permutation test statistics is skipped
#'    and the algorithm proceeds with calculation
#'    of the consensus null distribution
#' @param estP0args A list of arguments passed on to the estP0 function
#'
#' @details Efron (2007) centers the observations in each group prior
#'  to permutation. As permutations will remove any genuine group differences
#'   anyway, we skip this step by default.\\ If zValues = FALSE,
#'   the density is fitted on the original test statistics rather than converted
#'   to z-values. This unlocks the procedure for test statistics with unknown
#'   distributions, but may be numerically less stable.
#'
#' @return A list with entries
#' \item{zValsMat}{Permutation Z-values}
#' \item{zValObs}{Observed Z-values}
#' \item{Cperm}{(optional) An estimated covariance matrix
#'  of binned test statistics}
#' \item{weightStrat}{The weighting strategy }
#' \item{Fdr, fdr}{The estimated tail-area and local false discovery rates}
#' \item{consensus}{The consensus null distribution}
#' @export
#' @importFrom stats pnorm qnorm
#' @examples
#' p = 100; n = 50
#' x = rep(c(0,1), each = n/2)
#' mat = cbind(
#' matrix(rnorm(n*p/10, mean = 5+x),n,p/10), #DA
#' matrix(rnorm(n*p*9/10, mean = 5),n,p*9/10) #Non DA
#' )
#' fdrRes = fdrCorrect(mat, x)
#' fdrRes$p0
#' #Indeed close to 0.9
#' estFdr = fdrRes$Fdr
#' #The estimated tail-area false discovery rates.
#'
#' #With another type of test. Need to supply quantile function in this case
#' fdrResLm = fdrCorrect(mat, x, B = 5e1,
#' test = function(x, y){
#' c(summary(lm(y~x))$coef["x1","t value"], fit$df.residual)},
#' quantileFun = function(t){pt(a = t[1], df = t[2])})
#'
#' With a test statistic without known null distribution(for small samples)
#' fdrResKruskal = fdrCorrect(mat, x, B = 5e1,
#' test = function(x, y){
#' kruskal.test(y~x)$statistic}, zValue = FALSE)
fdrCorrect = function(Y, x, B = 5e2L, test = "wilcox.test", argList = list(),
                      quantileFun = "pwilcox", zValues = TRUE,
                      testPargs = list(m = table(x)[1],
                                       n = table(x)[2]),
                      z0Quant = pnorm(c(-1,1)), gridsize = 801L,
                      weightStrat = "LH", maxIter = 1000L, tol = 1e-4,
                      cPerm = FALSE, nBins = 82L, binEdges = c(-4.1,4.1),
                      center = FALSE, zVals = NULL,
                      estP0args = list(z0quantRange = seq(0.05,0.45, 0.0125),
                                       smooth.df = 3)){
  if(is.character(test) && test == "t.test"){
  testPargs = list()
  quantileFun = "pt.edit"
  }
  p = ncol(Y)
if(is.null(zVals)){
#Test statistics
testStats = getTestStats(Y = Y, center = center, test = test,
                         x = x, B = B, argList = argList)

#Permuation z-values
zValsMat = apply(testStats$statsPerm,
                 MARGIN = switch(test, "t.test" = c(3,2), 2),
                 function(stats){
  qnormIf(zValues = zValues,
          quantCorrect(do.call(quantileFun,c(list(q = stats), testPargs))))
})
#Observed z-values
zValObs = qnormIf(apply(matrix(testStats$statObs, ncol = p),2, function(stats){
  quantCorrect(do.call(quantileFun, c(list(q = stats), testPargs)))
  }), zValues = zValues)
} else {
    if(!zValues) stop("Z-values supplied, then also let zValues be TRUE!")
  zValObs = zVals$zValObs; zValsMat = zVals$zValsMat
}

#Permuation correlation matrix of binned test statistics
Cperm = if(cPerm){getCperm(zValsMat, nBins = nBins, binEdges = binEdges)
} else {
        NULL
    }

#Consensus distribution estimation
consensus = getG0(zValObs = zValObs, zValsMat = zValsMat, z0Quant = z0Quant,
                  gridsize = gridsize, maxIter = maxIter, tol = tol,
                  weightStrat = weightStrat, estP0args = estP0args)

#False discovery Rates
FdrList = do.call(getFdr, c(list(zValObs = zValObs, p = p), consensus))

names(zValObs) = names(FdrList$Fdr) = names(FdrList$fdr) = colnames(Y)
c(list(zValsMat = zValsMat, zValObs = zValObs, Cperm = Cperm,
       weightStrat = weightStrat), FdrList, consensus)
}
