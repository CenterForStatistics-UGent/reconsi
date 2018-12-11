#' Correct quantiles by not returning 0 or 1
#' @param quants A vector of quantiles

#'
#' @return The same vector of quantiles but without 0 or 1 values
quantCorrect = function(quants){
    quants[quants==1] = 1-.Machine$double.eps
    quants[quants==0] = .Machine$double.eps
    return(quants)
}
