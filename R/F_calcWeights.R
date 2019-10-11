#' Obtain weights as posterior probabilities to calculate the consensus null
#'
#' @param logDensPerm A matrix with B rows of logged density estimates of the B
#'    permutation distributions
#' @param weightStrat A character string, either "LH" for likelihoods based
#'    on the central region, or "LHw" for weighted likelihoods
#' @param fdr A vector of local false discovery rates
#'    along the support of the kernel
#' @param zIndR An indicator vector indicating the z-values
#'    within the central region
#'
#' @return A vector of weights of length B
calcWeights = function(logDensPerm, fdr){
  weights = exp(stabExp(colSums(logDensPerm*fdr)))
  weights/sum(weights)
}
