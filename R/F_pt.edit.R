#' A custom function to calculate the distribution function of the t-test
#'  statistic. For internal use only
#'@importFrom stats pt
#'@param q vector with t-statistic and degrees of freedom
#'@return corresponding quantile
pt.edit = function(q){
  pt(q = q[1], df = q[2])
}
