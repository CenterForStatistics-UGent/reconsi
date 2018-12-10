#' A custom function to calculate the distribution function of the t-test
#'  statistic. For internal use only
pt.edit = function(q){
  pt(q = q[1], df = q[2])
}
