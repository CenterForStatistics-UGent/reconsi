#' A custom function to calculate the distribution function of the t-test
#'  statistic. For internal use only
#'@importFrom stats pt
#'@param q vector with t-statistic and degrees of freedom
#'@param p a vector with quantile and degrees of freedom
#'@return corresponding quantile
pt.edit = function(q){
  pt(q = q[1], df = q[2])
}
qt.edit = function(p){
  qt(p = p[1], df = p[2])
}
