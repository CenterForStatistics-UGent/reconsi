#' A function to test for differential absolute abundance on a phyloseq object
#'@param physeq A phyloseq object
#'@param groupName A character string, the name of a variable in physeq
#'indicating the grouping factor
#'@param FCname A character string, the name of a variable in physeq
#'containing the total cell count
#'@param ... passed on to the testDAA() function
#'
#'@return See the fdrCorrect() function
#'@export
#'@import phyloseq
#'@examples
#'testVanDePutte = testDAA(
testDAA.phylo = function(physeq, groupName, FCname, ...){
otuTab  = if(taxa_are_rows(physeq)) t(as(otu_table(physeq), "matrix")) else
    as(otu_table(physeq), "matrix")
testDAA(Y = otuTab, FC = get_variable(physeq, FCname),
        x = factor(get_variable(physeq, groupName)), ...)
}
