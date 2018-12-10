#' A function to obtain a t-test statistic efficiently. For internal use only
#'
#' @param y1,y2 vectors of obsereved values in the two groups
#' @param mm,nn number of observations in the corresponding groups
#'
#' @return A list with items
#' \item{tstat}{The t-test statistic}
#' \item{df}{The degrees of freedom (Welch approximation)}#'
getTstat = function(y1, y2, mm, nn){
  var1 = var(y1)
  var2 = var(y2)
  c(tstat = (mean(y1) - mean(y2))/sqrt(var1/mm + var2/nn),
    df = (var1/mm + var2/nn)^2/((var1/mm)^2/(mm - 1) + (var2/nn)^2/(nn - 1)))
}
