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
#' @param ... passed on to the reconsi() function
#'
#'@return See the reconsi() function

#' @rdname testDAA
#' @export
setGeneric("testDAA", function(Y, ...) standardGeneric("testDAA"))

#'@export
#'@import phyloseq
#'@import methods
#'@rdname testDAA
#'@examples
#'#Test for phyloseq object
#'library(phyloseq)
#'VandeputtePruned = prune_samples(Vandeputte,
#'samples = sample_names(Vandeputte)[20:40])
#'testVanDePutte = testDAA(VandeputtePruned, "Health.status", "absCountFrozen",
#'B = 15)
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
#' testMat = testDAA(as(otu_table(VandeputtePruned), "matrix"),
#' get_variable(VandeputtePruned, "Health.status"),
#' get_variable(VandeputtePruned,"absCountFrozen"), B = 15)
setMethod("testDAA", "matrix", function(Y, FC, x, S = rowSums(Y), ...){
    idSam = S>0
    if(min(table(x[idSam]))<2L){
        return(NULL)
        }
    reconsi(Y[idSam,colSums(Y)>0], x=x[idSam], ...)
})
