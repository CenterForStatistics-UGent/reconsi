#' A function to calculate observed and permuation z-statistics
#' on a n-by-p matrix of observations
#'
#' @param Y The nxp data matrix
#' @param center a boolean, should data be centered prior to permuation
#' @param test A function name, possibly user defined. See details.
#' @param x A vector defining the groups. Will be coerced to factor.
#' @param B an integer, the number of permuations
#' @param argList A list of further arguments passed on to the test function
#'
#' @return A list with components
#' \item{statObs}{A vector of length p of observed test statistics}
#' \item{statsPerm}{A p-by-B matrix of permutation test statistics}
#'
#' @importFrom permute shuffleSet
#'
#' @details For test "wilcox.test" and "t.test",
#' fast custom implementations are used. Other functions can be supplied
#' but must accept a y outcome variable, a x as grouping variable, and possibly
#' a list of other arguments. It must return all arguments needed to evaluate
#' its quantile function if z-values are to be used.
getTestStats = function(Y, center, test = "wilcox.test", x, B,
                        argList){
  x = factor(x)
  Ycenter = Y

  #enumerate B unique ways to group

  permDesign = shuffleSet(length(x),B)
  if(center) {
    for (ii in unique(x)){Ycenter[x==ii,] = scale(Ycenter[x==ii,],
                                                  center = TRUE,
                                                  scale = FALSE)}
  }
  if(is.character(test) && (test %in% c("wilcox.test","t.test"))){
    if(nlevels(x)>2){stop("Wilcoxon rank sum test and t-test only apply
                          to two groups! \n Try 'kruskal.test()' or 'lm()'.")}
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
    #YRankedCenter = t(if(center) apply(Ycenter, 2, rank) else YRanked)
    statsPerm = - nFac + vapply(seq_len(B), function(ii){
      Yranked = t(apply(Y, 2, rank, ties.method = "random"))
      #Random tiebreaking, important to combat test statisitc discreteness
      rowSums(Yranked[,permDesign[ii,mmSeries]])
    }, FUN.VALUE = statObs)
  } else if(test == "t.test"){
    statObs = apply(Y, 2, function(dat){
getTstat(y1 = dat[xLog], y2 = dat[!xLog], mm = mm, nn = nn)
    })
    statsPerm = vapply(seq_len(B), FUN.VALUE = statObs,
                       function(ii){
      xSam = xLog[permDesign[ii,]]
      apply(Ycenter, 2, function(dat){
        getTstat(y1 = dat[xSam],y2 = dat[!xSam], mm = mm, nn = nn)
      })
    })
  }
    } else {
  testFun = match.fun(test)
  statObs = apply(Y,2, function(y){
      do.call(testFun, c(list(y = y, x = x), argList))
  })
  n = nrow(Y)
  statsPerm = vapply(seq_len(B),
                     FUN.VALUE = statObs,
                     function(ii){
    apply(Ycenter[permDesign[ii,],],2, function(y){
    do.call(testFun, c(list(y = y, x = x), argList))
  })})
  }

  return(list(statObs = statObs, statsPerm = statsPerm))
}
