#' Convert to cdf-values or not
#' @param zValues a boolean, is conversion to z-values required
#' @param distFun the quantile function
#' @param testPargs additional arguments to the quantile functions
#' @param stats the vector of test statistics
#'
#' @return z-values or original statistics
quantileIf = function(zValues, distFun, testPargs, stats){
  if(zValues) {quantCorrect(do.call(distFun,
                                          c(list(q = stats), testPargs)))
    } else {stats}
}
