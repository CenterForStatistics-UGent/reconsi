#' A function to test for differential absolute abundance on a phyloseq object
#' @param Y A phyloseq object, or a data matrix with samples in the rows and
#'OTUs in the columns
#' @param groupName A character string, the name of a variable in physeq
#'indicating the grouping factor
#' @param FCname A character string, the name of a variable in physeq
#'containing the total cell count
#' @param FC a vector of length n with total flow cytometry cell counts
#' @param x a grouping factor of length n
#' @param S a vector of library sizes. Will be calculated f
#' @param tieBreak A character string, specifying the tie break doctrine to use.
#' @param ... passed on to the fdrCorrect() function
#'
#'@return See the fdrCorrect() function


#' @rdname testDAA
#' @export
setGeneric("testDAA", function(Y, ...) standardGeneric("testDAA"))

#'@export
#'@import phyloseq
#'@rdname testDAA
#'@examples
#'#Test for phyloseq object
#'testVanDePutte = testDAA(Vandeputte, "Health.status", "absCountFrozen")
setMethod("testDAA", "phyloseq", function(Y, groupName, FCname, ...){
otuTab  = if(taxa_are_rows(Y)) t(as(otu_table(Y), "matrix")) else
    as(otu_table(Y), "matrix")
testDAA(Y = otuTab, FC = get_variable(Y, FCname),
        x = factor(get_variable(Y, groupName)), ...)
})

#' @export
#' @rdname testDAA
#' @examples
#' #Test for matrix
#' library(phyloseq)
#' testMat = testDAA(as(otu_table(Vandeputte), "matrix"),
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
