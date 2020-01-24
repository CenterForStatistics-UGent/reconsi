#' A function to calculate observed and permuation z-statistics
#' on a n-by-p matrix of observations
#'
#' @param Y The nxp data matrix
#' @param center a boolean, should data be centered prior to permuation
#' @param test A function name, possibly user defined. See details.
#' @param x A vector defining the groups. Will be coerced to factor.
#' @param B an integer, the number of permuations
#' @param argList A list of further arguments passed on to the test function
#' @param tieBreakRan A boolean, should ties of permutation test statistics
#'  be broken randomly? If not, midranks are used
#' @param replace A boolean. If FALSE, samples are permuted
#' (resampled without replacement), if TRUE the samples are bootstrapped
#' (resampled with replacement)
#'
#' @return A list with components
#' \item{statObs}{A vector of length p of observed test statistics}
#' \item{statsPerm}{A p-by-B matrix of permutation test statistics}
#' @importFrom stats runif
#'
#' @details For test "wilcox.test" and "t.test",
#' fast custom implementations are used. Other functions can be supplied
#' but must accept a y outcome variable, a x as grouping variable, and possibly
#' a list of other arguments. It must return all arguments needed to evaluate
#' its quantile function if z-values are to be used.
getTestStats = function(Y, center, test = "wilcox.test", x, B,
                        argList, tieBreakRan = FALSE, replace = FALSE){
  #Sample size
  n = nrow(Y)
  #enumerate B ways to permute/combine
  permDesign = if(replace) matrix(ceiling(n*runif(n*B)),
                                  nrow = n, ncol = B) else
      vapply(integer(B), FUN.VALUE = integer(n),function(x) sample.int(n = n))
  if(center){
    if(replace){
      Ycenter = scale(Y, center = TRUE, scale = FALSE)
    } else {
    for (ii in unique(x)){
      Y[x==ii,] = scale(Y[x==ii,], center = TRUE, scale = FALSE)
      }
    }
  } else {
      Ycenter = Y
      }
  if(is.character(test) && (test %in% c("wilcox.test","t.test"))){
    xLog = x==names(table(x))[1]
    nn = table(x)[2]
    mm = table(x)[1]
  if(test=="wilcox.test"){ #Shortcuts possbile in case of Wilcoxon test
    nFac = mm*(mm + 1)/2
    YRanked = apply(Y, 2, rank)
    #Observed test statistic
    statObs = colSums(YRanked[xLog,]) - nFac
    #Permuation test statistics
    mmSeries = seq_len(mm)
    statsPerm = - nFac + vapply(seq_len(B), function(ii){
      if(tieBreakRan){
          YRanked = apply(Y, 2, rank, ties.method = "random")
      #Random tiebreaking, important to combat test statistic discreteness
      }
        colSums(YRanked[permDesign[mmSeries,ii],])
    }, FUN.VALUE = statObs)
  } else if(test == "t.test"){
    statObs = apply(Y, 2, function(dat){
getTstat(y1 = dat[xLog], y2 = dat[!xLog], mm = mm, nn = nn)
    })
    statsPerm = vapply(seq_len(B), FUN.VALUE = statObs,
                       function(ii){
      xSam = xLog[permDesign[,ii]]
      apply(Ycenter, 2, function(dat){
        getTstat(y1 = dat[xSam], y2 = dat[!xSam], mm = mm, nn = nn)
      })
    })
  }
    } else {
  testFun = match.fun(test)
  statObs = apply(Y,2, function(y){
      do.call(testFun, c(list(y = y, x = x), argList))
  })
  statsPerm = vapply(seq_len(B),
                     FUN.VALUE = statObs,
                     function(ii){
    apply(Ycenter[permDesign[,ii],],2, function(y){
    do.call(testFun, c(list(y = y, x = x), argList))
  })})
    }
  return(list(statObs = statObs, statsPerm = statsPerm))
}
