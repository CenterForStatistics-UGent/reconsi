#' Convert to z-values or not
#' @param testStats a vector of test statistics
#' @param zValues a boolean, is conversion to z-values required
#'
#' @return z-values or original statistics
qnormIf = function(testStats, zValues){
  if(zValues) qnorm(testStats) else testStats
}
