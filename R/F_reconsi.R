#' Perform simultaneous inference through
#' collapsed resampling null distributions
#' @param Y the matrix of sequencing counts
#' @param x a grouping factor. If provided, this grouping factor is permuted.
#' Otherwise a bootstrap procedure is performed
#' @param B the number of resampling instances
#' @param test Character string, giving the name of the function
#'  to test for differential absolute abundance.
#'  Must accept the formula interface. Features with tests resulting in
#'  observed NA test statistics will be discarded
#' @param argList A list of arguments, passed on to the testing function
#' @param distFun the distribution function of the test statistic, or its name.
#' Must at least accept an argument named 'q', 'p' and 'x' respectively.
#' @param zValues A boolean, should test statistics be converted to z-values.
#' See details
#' @param testPargs A list of arguments passed on to distFun
#' @param z0Quant A vector of length 2 of quantiles of the null distribution,
#'    in between which only null values are expected
#' @param gridsize The number of bins for the kernel density estimates
#' @param maxIter An integer, the maximum number of iterations in the estimation
#'    of the null distribution
#' @param tol The tolerance for the infinity norm of the central borders
#'    in the iterative procedure
#' @param center A boolean, should observations be centered
#'    in each group prior to permuations? See details.
#' @param scale a boolean, should data be scaled prior to resampling
#' @param replace A boolean. Should resampling occur with replacement (boostrap)
#' or without replacement (permutation)
#' @param zVals An optional list of observed (statObs) and
#' resampling (statsPerm) z-values. If supplied, the calculation
#'    of the observed and resampling test statistics is skipped
#'    and the algorithm proceeds with calculation
#'    of the consensus null distribution
#' @param estP0args A list of arguments passed on to the estP0 function
#' @param resamZvals A boolean, should resampling rather than theoretical null
#' distributions be used?
#' @param smoothObs A boolean, should the fitted rather than estimated observed
#' distribution be used in the Fdr calculation?
#' @param tieBreakRan A boolean, should ties of resampling test statistics
#'  be broken randomly? If not, midranks are used
#' @param pi0 A known fraction of true null hypotheses. If provided,
#' the fraction of true null hypotheses will not be estimated.
#' Mainly for oracle purposes.
#'@details Efron (2007) centers the observations in each group prior
#'  to permutation. As permutations will remove any genuine group differences
#'   anyway, we skip this step by default. If zValues = FALSE,
#'   the density is fitted on the original test statistics rather than converted
#'   to z-values. This unlocks the procedure for test statistics with unknown
#'   distributions, but may be numerically less stable.
#' @return A list with entries
#' \item{statsPerm}{Resampling Z-values}
#' \item{statObs}{Observed Z-values}
#' \item{distFun}{Density, distribution and
#' quantile function as given}
#' \item{testPargs}{Same as given}
#' \item{zValues}{A boolean, were z-values used?}
#' \item{resamZvals}{A boolean, were the resampling null distribution used?}
#' \item{cdfValObs}{Cumulative distribution function evaluation
#' of observed test statistics}
#' \item{p0estimated}{A boolean, was the fraction of true null hypotheses
#' estimated from the data?}
#' \item{Fdr,fdr}{Estimates of tail-area and local false discovery rates}
#' \item{p0}{Estimated or supplied fraction of true null hypotheses}
#' \item{iter}{Number of iterations executed}
#' \item{fitAll}{Mean and standard deviation estimated collapsed null}
#' \item{PermDensFits}{Mean and standard deviations of resamples}
#' \item{convergence}{A boolean, did the iterative algorithm converge?}
#' \item{zSeq}{Basis for the evaluation of the densities}
#' \item{weights}{weights of the resampling distributions}
#' \item{zValsDensObs}{Estimated overall densities, evaluated in zSeq}
#' @export
#' @importFrom stats pnorm qnorm
#' @note Ideally, it would be better to only use unique resamples, to avoid
#' unnecesarry replicated calculations of the same test statistics. Yet this
#' issue is almost alwyas ignored in practice; as the sample size grows it also
#' becomes irrelevant. Notice also that this would require to place weights in
#' case of the bootstrap, as some bootstrap samples are more likely than others.
#' @examples
#' #Important notice: low number of resamples B necessary to keep
#' # computation time down, but not recommended. Pray set B at 200 or higher.
#' p = 50; n = 20; B = 5e1
#' x = rep(c(0,1), each = n/2)
#' mat = cbind(
#' matrix(rnorm(n*p/10, mean = 5+x), n, p/10), #DA
#' matrix(rnorm(n*p*9/10, mean = 5), n, p*9/10) #Non DA
#' )
#' fdrRes = reconsi(mat, x, B = B)
#' fdrRes$p0
#' #Indeed close to 0.9
#' estFdr = fdrRes$Fdr
#' #The estimated tail-area false discovery rates.
#'
#' #With another type of test. Need to supply quantile function in this case
#' fdrResLm = reconsi(mat, x, B = B,
#' test = function(x, y){
#' fit = lm(y~x)
#' c(summary(fit)$coef["x","t value"], fit$df.residual)},
#' distFun = function(q){pt(q = q[1], df = q[2])})
#'
#' #With a test statistic without known null distribution(for small samples)
#' fdrResKruskal = reconsi(mat, x, B = B,
#' test = function(x, y){
#' kruskal.test(y~x)$statistic}, zValues = FALSE)
#'
#' #Provide an additional covariate through the 'argList' argument
#' z = rpois(n , lambda = 2)
#' fdrResLmZ = reconsi(mat, x, B = B,
#' test = function(x, y, z){
#' fit = lm(y~x+z)
#' c(summary(fit)$coef["x","t value"], fit$df.residual)},
#' distFun = function(q){pt(q = q[1], df = q[2])},
#' argList = list(z = z))
#'
#' #When nog grouping variable is provided, a bootstrap is performed
#' matBoot = cbind(
#' matrix(rnorm(n*p/10, mean = 1), n, p/10), #DA
#' matrix(rnorm(n*p*9/10, mean = 0), n, p*9/10) #Non DA
#' )
#' fdrResBoot = reconsi(matBoot, B = B,
#' test = function(y, x){testRes = t.test(y, mu = 0, var.equal = TRUE);
#' c(testRes$statistic, testRes$parameter)},
#' distFun = function(q){pt(q = q[1], df = q[2])},
#' center = TRUE, replace = TRUE)
reconsi = function(Y, x = NULL, B = 1e3L, test = "wilcox.test",
                   argList = list(),
                      distFun ="pnorm",
                      zValues = TRUE, testPargs = list(),
                      z0Quant = pnorm(c(-1,1)), gridsize = 801L,
                      maxIter = 1000L, tol = 1e-8,
                      center = FALSE, replace = is.null(x), zVals = NULL,
                      estP0args = list(z0quantRange = seq(0.05,0.45, 0.0125),
                                       smooth.df = 3), resamZvals = FALSE,
                      smoothObs = TRUE, scale = FALSE,
                      tieBreakRan = FALSE, pi0 = NULL){
    #Basic checks
    stopifnot(is.matrix(Y), is.list(argList), is.logical(center),
              is.logical(smoothObs),
              is.logical(resamZvals), is.logical(replace), is.numeric(z0Quant),
              is.numeric(tol), is.numeric(maxIter), is.numeric(gridsize),
              is.numeric(B), is.list(estP0args))
    if(is.function(test)){

    }
    else if(is.character(test)){
        if(test %in% c("wilcox.test","t.test")){
            if(is.null(x) || NCOL(x)!=1){
                stop("Provide single grouping variable for Wilcoxon rank sum test or t-test")
            }
            x = factor(x)
            if(nlevels(x)>2){stop("Wilcoxon rank sum test and t-test only apply to two groups! \n Try 'kruskal.test()' or 'lm()'.")}
        }
      if(test == "t.test"){
        distFun = "ptEdit"
        if(!zValues) {
            stop("With t-tests, the test statistic must be converted to zValues. Please set zValues = TRUE")
        }
        } else if(test=="wilcox.test"){
        testPargs = list(m = table(x)[1], n = table(x)[2])
        distFun = "pwilcox"
        }
    }
    if(!"q" %in% names(formals(distFun))){
        stop("Distribution function must accept arguments named 'q'\n")
    }
    if(!"y" %in% names(formals(test)) && !test %in% c("t.test", "wilcox.test")){
        stop("Test function must accept 'y' argument\n")
    }
    if(!is.matrix(Y)){
        stop("Please provide a data matrix as imput!\n",
               if(is.vector(Y)) "Multiplicity correction only needed when testing multiple hypotheses.")
    }
    if(!is.null(x) && NROW(x)!=nrow(Y)){
        stop("Length of grouping factor does not correspond to number of rows of data matrix!\n")
    }
    if(ncol(Y)<30){
     warning("Less than 30 hypotheses tested, multplicity correction may be unreliable.\nConsider using p.adjust(..., method ='BH').",
             immediate. = TRUE)
    }
    if(nrow(Y)<7){
        warning("Small sample size may not yield sufficient unique resampling instances.\nConsider using p.adjust(..., method ='BH').",
                immediate. = TRUE)
    }
      p = ncol(Y)
    if(is.null(zVals)){
    #Test statistics
    testStats = getTestStats(Y = Y, center = center, test = test,
                             x = x, B = B, argList = argList, scale = scale,
                             tieBreakRan = tieBreakRan, replace = replace)
    #Observed cdf values
    cdfValObs = apply(matrix(testStats$statObs, ncol = p), 2, function(stats){
        quantCorrect(do.call(distFun, c(list(q = stats), testPargs)))
    })
    if(zValues){
    #If procedure works with z-values rather than raw test statistics,
        #convert to z-values
    #Resample cdf-values
    cdfValsMat =  apply(testStats$statsPerm,
                      MARGIN = if(length(dim(testStats$statsPerm))==3) c(2,3)
                      else 2, function(stats){
                          quantCorrect(do.call(distFun,
                                               c(list(q = stats), testPargs)))
                      })
    #Observed z-values
    statObs = qnorm(
        if(resamZvals) quantCorrect(vapply(seq_along(cdfValObs),
                                          FUN.VALUE = double(1),
                                          function(i){
            (sum(cdfValObs[i] > cdfValsMat[i,])+1L)/(B+2L)
        })) else cdfValObs
    )
    #Resample z-values
    statsPerm = qnorm(
        if(resamZvals) vapply(seq_len(B), FUN.VALUE = statObs, function(b){
            quantCorrect(vapply(seq_along(cdfValObs), FUN.VALUE = double(1),
                                function(i){
                (sum(cdfValsMat[i,b] > cdfValsMat[i,-b])+1L)/(B+2L)
            }))
        }) else cdfValsMat
    )
    } else {
        #If procedure works on raw test statistics, not many conversions are needed
    #Observed statistics
    statObs = if(is.matrix(testStats$statObs)) testStats$statObs[1,] else
      testStats$statObs
    #Permuation statistics
    statsPerm = if(length(dim(testStats$statsPerm))==3)
        testStats$statsPerm[1,,] else testStats$statsPerm
    }
    } else {
        #If test statistics give, recover them
      statObs = zVals$statObs; statsPerm = zVals$statsPerm
      cdfValObs = zVals$cdfValObs
     }
      if(zValues) {
          testPargs = list()
          }
    #Consensus distribution estimation
    consensus = getG0(statObs = statObs, statsPerm =  statsPerm,
                      z0Quant = z0Quant, gridsize = gridsize, maxIter = maxIter,
                      tol = tol, estP0args = estP0args, testPargs = testPargs,
                      B = B, p = p, pi0 = pi0)
    #False discovery Rates
    FdrList = do.call(getFdr,
                      c(list(statObs = statObs, p = p, smoothObs = smoothObs),
                        consensus))
    consensus$fdr = NULL
    names(statObs) = names(FdrList$Fdr) = names(FdrList$fdr) = colnames(Y)
    c(list(statsPerm = statsPerm, statObs = statObs, zValues = zValues,
           resamZvals = resamZvals, cdfValObs = cdfValObs,
           testPargs = testPargs, distFun = distFun,
           p0estimated = is.null(pi0),
           resamDesign = if(is.null(zVals)) testStats$resamDesign),
      FdrList, consensus)
}
