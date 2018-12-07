#' A function to perform the complete differential absolute abundance analysis
#' @param Y a n-by-p sequence count tabel with taxa in the columns
#' @param FC a vector of length n with total flow cytometry cell counts
#' @param x a grouping factor of length n
#' @param S a vector of library sizes. Will be calculated f
#' @param tieBreak A character string, specifying the tie break doctrine to use.
#' @param ... other arguments passed on to the fdrCorrect function

#'
#' @return see ?fdrCorrect()
#'
#' @details This function breaks ties, and filters empty rows and columns prior to statistical testing
#' @export
#' @examples
testDAA = function(Y, FC, x, S = rowSums(Y), tieBreak = "none", ...){
  idSam = S>0
  tieBrokenY = tiebreak(Y = Y[idSam, ], FC = if(length(FC)==1) FC else FC[idSam], tieBreak = tieBreak, S = S[idSam])
  if(min(table(x[idSam]))<2L) {return(NULL)}
  fdrCorrect(tieBrokenY[,colSums(tieBrokenY)>0], x=x[idSam], ...)
}