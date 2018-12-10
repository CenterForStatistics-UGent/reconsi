#' A function to test for differential absolute abundance on a phyloseq object
#'@param physeq A phyloseq object
#'@param groupName A character string, the name of a variable in physeq
#'indicating the grouping factor
#'@param FCname A character string, the name of a variable in physeq
#'containing the total cell count
#'@param ... passed on to the testDAA() function
#'
#'@return See the fdrCorrect() function
#'@import phyloseq
#'@examples
#'testVanDePutte = testDAA(Vandeputte, "Health.status", "absCountFrozen")
#' @rdname testDAA
#' @export
setGeneric("testDAA", function(Y, ...) standardGeneric("testDAA"))

#'@export
#'@rdname testDAA
setMethod("testDAA", "phyloseq", function(Y, groupName, FCname, ...){
otuTab  = if(taxa_are_rows(Y)) t(as(otu_table(Y), "matrix")) else
    as(otu_table(Y), "matrix")
testDAA(Y = otuTab, FC = get_variable(Y, FCname),
        x = factor(get_variable(Y, groupName)), ...)
})
#' Perform the complete differential absolute abundance analysis
#' @param Y a n-by-p sequence count tabel with taxa in the columns
#' @param FC a vector of length n with total flow cytometry cell counts
#' @param x a grouping factor of length n
#' @param S a vector of library sizes. Will be calculated f
#' @param tieBreak A character string, specifying the tie break doctrine to use.
#' @param ... other arguments passed on to the fdrCorrect function

#' @return see ?fdrCorrect()
#'
#' @details This function breaks ties, and filters empty rows and columns
#' prior to statistical testing
#' @export
#' @rdname testDAA
#' @examples
#' testMat = testDAA(as.matrix(otu_table(Vandeputte)),
#' get_variable(Vandeputte, "Health.status"),
#' get_variable(Vandeputte,"absCountFrozen"))
setMethod("testDAA", "matrix", function(Y, FC, x, S = rowSums(Y), tieBreak = "none", ...){
    idSam = S>0
    tieBrokenY = tiebreak(Y = Y[idSam, ],
                          FC = if(length(FC)==1) FC else FC[idSam],
                          tieBreak = tieBreak, S = S[idSam])
    if(min(table(x[idSam]))<2L) {return(NULL)}
    fdrCorrect(tieBrokenY[,colSums(tieBrokenY)>0], x=x[idSam], ...)
})
