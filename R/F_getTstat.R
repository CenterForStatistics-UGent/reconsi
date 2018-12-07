#' A function to obtain a t-test statistic efficiently. For internal use only
getTstat = function(y1, y2, mm, nn){
  var1 = var(y1)
  var2 = var(y2)
  c(tstat = (mean(y1) - mean(y2))/sqrt(var1/mm + var2/nn),
    df = (var1/mm + var2/nn)^2/((var1/mm)^2/(mm - 1) + (var2/nn)^2/(nn - 1)))
}