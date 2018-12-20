#' Perform simultaneous inference, correcting for correlation between tests.
#' @param Y the matrix of sequencing counts
#' @param x a grouping factor
#' @param B the number of permutations
#' @param test Character string, giving the name of the function
#'  to test for differential absolute abundance.
#'  Must accept the formula interface
#' @param argList A list of arguments, passed on to the testing function
#' @param distFun the distribution function of the test statistic,
#' or its name. Must at least accept an argument named q.
#' @param zValues A boolean, should test statistics be converted to z-values.
#' See details
#' @param testPargs A list of arguments passed on to distFun
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
#' \item{zValObsPerm}{Permutation observed z-values}
#' \item{Cperm}{(optional) An estimated covariance matrix
#'  of binned test statistics}
#' \item{weightStrat}{The weighting strategy }
#' \item{Fdr, fdr}{The estimated tail-area and local false discovery rates}
#' \item{consensus}{The consensus null distribution}
#' @export
#' @importFrom stats pnorm qnorm
#' @examples
#' p = 200; n = 50; B = 5e1
#' x = rep(c(0,1), each = n/2)
#' mat = cbind(
#' matrix(rnorm(n*p/10, mean = 5+x),n,p/10), #DA
#' matrix(rnorm(n*p*9/10, mean = 5),n,p*9/10) #Non DA
#' )
#' fdrRes = fdrCorrect(mat, x, B = B)
#' fdrRes$p0
#' #Indeed close to 0.9
#' estFdr = fdrRes$Fdr
#' #The estimated tail-area false discovery rates.
#'
#' #With another type of test. Need to supply quantile function in this case
#' fdrResLm = fdrCorrect(mat, x, B = B,
#' test = function(x, y){
#' fit = lm(y~x)
#' c(summary(fit)$coef["x1","t value"], fit$df.residual)},
#' distFun = function(q){pt(q = q[1], df = q[2])})
#'
#' #With a test statistic without known null distribution(for small samples)
#' fdrResKruskal = fdrCorrect(mat, x, B = B,
#' test = function(x, y){
#' kruskal.test(y~x)$statistic}, zValues = FALSE)
#'
#' #Provide an additional covariate through the 'argList' argument
#' z = rpois(n , lambda = 2)
#' fdrResLmZ = fdrCorrect(mat, x, B = B,
#' test = function(x, y, z){
#' fit = lm(y~x+z)
#' c(summary(fit)$coef["x1","t value"], fit$df.residual)},
#' distFun = function(q){pt(q = q[1], df = q[2])},
#' argList = list(z = z))
fdrCorrect = function(Y, x, B = 1e3L, test = "wilcox.test", argList = list(),
                      distFun = "pwilcox", zValues = TRUE,
                      testPargs = list(),
                      z0Quant = pnorm(c(-1,1)), gridsize = 801L,
                      weightStrat = "LH", maxIter = 1000L, tol = 1e-5,
                      cPerm = FALSE, nBins = 82L, binEdges = c(-4.1,4.1),
                      center = FALSE, zVals = NULL,
                      estP0args = list(z0quantRange = seq(0.05,0.45, 0.0125),
                                       smooth.df = 3), permZvals = FALSE,
                      normAsump = FALSE, smoothObs = FALSE){
  if(is.character(test)){
      if(test == "t.test"){
        testPargs = list()
        distFun = "pt.edit"
        } else if(test=="wilcox.test"){
        testPargs = list(m = table(x)[1], n = table(x)[2])
        }
  }
if(!"q" %in% formalArgs(match.fun(distFun))) {
    stop("Quantile function must accept arguments named 'q'\n")
}
if(!is.matrix(Y)){
    stop(paste("Please provide a data matrix as imput!\n",
               if(is.vector(Y)) "Multiplicity correction only needed when
               testing multiple hypotheses."))
}
if(length(x)!=nrow(Y)){
    stop("Length of grouping factor does not correspond to number of
         rows of data matrix!\n")
}
if(ncol(Y)<30){
 warning("Less than 30 hypotheses tested, multplicity correction may
         be unreliable.\nConsider using p.adjust(..., method ='BH').",
         immediate. = TRUE)
}
  p = ncol(Y)
if(is.null(zVals)){
#Test statistics
testStats = getTestStats(Y = Y, center = center, test = test,
                         x = x, B = B, argList = argList)

#Observed cdf values
cdfValObs = apply(matrix(testStats$statObs, ncol = p), 2, function(stats){
  quantileIf(zValues = TRUE, distFun = distFun, stats = stats,
             testPargs = testPargs)
})

if(!(is.character(test) && test=="wilcox.test" && !zValues)){
#Permutation cdf-values
cdfValsMat =  apply(testStats$statsPerm,
                  MARGIN = if(length(dim(testStats$statsPerm))==3) c(2,3) else 2,
                  function(stats){
                    quantileIf(zValues = zValues, distFun = distFun,
                               stats = stats, testPargs = testPargs)
                  })
}
#Observed z-values
zValObs = if(zValues) {qnorm(
  if(permZvals) quantCorrect(sapply(seq_along(cdfValObs), function(i){
  (sum(cdfValObs[i] > cdfValsMat[,i])+1L)/(B+2L)
})) else cdfValObs
)
} else {testStats$statObs}
#Permuation z-values
zValsMat = if(zValues) {qnorm(
  if(permZvals) sapply(seq_len(B), function(b){
  quantCorrect(sapply(seq_along(cdfValObs), function(i){
    (sum(cdfValsMat[b,i] > cdfValsMat[-b,i])+1L)/(B+2L)
  }))
}) else cdfValsMat
)
} else testStats$statsPerm

} else {
    #if(!zValues) stop("Z-values supplied, then also let zValues be TRUE!")
  zValObs = zVals$zValObs; zValsMat = zVals$zValsMat
  cdfValObs = zVals$cdfValObs
 }

#Permuation correlation matrix of binned test statistics
Cperm = if(cPerm){getCperm(zValsMat, nBins = nBins, binEdges = binEdges)
} else {NULL}

quantileFun = if(zValues) "qnorm" else
  if(test=="t.test") "qt" else if (test=="wilcox.test") "qwilcox"

#Consensus distribution estimation
consensus = getG0(zValObs = zValObs, zValsMat =  zValsMat,
                  z0Quant = z0Quant, gridsize = gridsize, maxIter = maxIter,
                  tol = tol, weightStrat = weightStrat, estP0args = estP0args,
                  normAsump = normAsump, quantileFun = quantileFun,
                  testPargs = testPargs)

#False discovery Rates
FdrList = do.call(getFdr,
                  c(list(zValObs = zValObs,
                         p = p, smoothObs = smoothObs), consensus))

names(zValObs) = names(FdrList$Fdr) = names(FdrList$fdr) = colnames(Y)
c(list(zValsMat = zValsMat, zValObs = zValObs, Cperm = Cperm,
       weightStrat = weightStrat, zValues = zValues, permZvals = permZvals,
       cdfValObs = cdfValObs,
       densFun = if(zValues) "dnorm" else
         if(test=="t.test") "dt" else if (test=="wilcox.test") "dwilcox",
       testPargs = if(zValues) list() else testPargs,
       distFun = if(zValues) "pnorm" else
         if(test=="t.test") "pt" else if (test=="wilcox.test") "pwilcox",
       quantileFun = quantileFun),
  FdrList, consensus)
}
