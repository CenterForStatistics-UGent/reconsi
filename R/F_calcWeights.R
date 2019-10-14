#' Obtain weights as posterior probabilities to calculate the consensus null
#'
#' @param logDensPerm A matrix with B rows of logged density estimates of the B
#'    permutation distributions
#' @param fdr A vector of local false discovery rates
#'    along the support of the kernel
#'
#' @return A vector of weights of length B
calcWeights = function(logDensPerm, fdr){
  weights = exp(stabExp(colSums(logDensPerm*fdr)))
  weights/sum(weights)
}
