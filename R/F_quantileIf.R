#' Convert to cdf-values or not
#' @param zValues a boolean, is conversion to z-values required
#' @param quantileFun the quantile function
#' @param testPargs additional arguments to the quantile functions
#' @param stats the vector of test statistics
#'
#' @return z-values or original statistics
quantileIf = function(zValues, quantileFun, testPargs, stats){
  if(zValues) {quantCorrect(do.call(quantileFun,
                                          c(list(q = stats), testPargs)))
    } else {stats}
}
