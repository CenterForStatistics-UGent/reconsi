#' A custom function to calculate the distribution function of the t-test
#'  statistic. For internal use only
#'@importFrom stats pt qt
#'@param q a vector with t-statistic and degrees of freedom
#'@return A value between 0 and 1, the evaluation of the cdf
#'@export
pt.edit = function(q){
  pt(q = q[1], df = q[2])
}
#' A custom function to calculate the quantile function of the t-test
#'  statistic. For internal use only
#'@param p a vector with quantile and degrees of freedom
#'@return the corresponding quantile
qt.edit = function(p){
  qt(p = p[1], df = p[2])
}
