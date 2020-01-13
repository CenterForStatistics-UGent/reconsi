#' Obtain weights as posterior probabilities to calculate the consensus null
#'
#' @param logDensPerm A matrix with B rows of logged density estimates of the B
#'    permutation distributions
#' @param fdr A vector of local false discovery rates
#'    along the support of the kernel
#' @param priorProbs A vector of length B with prior probabilities of
#' resampling instances
#'
#' @return A vector of weights of length B
calcWeights = function(logDensPerm, fdr, priorProbs){
  weights = exp(stabExp(colSums(logDensPerm*fdr)))
  if(!is.null(priorProbs)) weights = weights*priorProbs
  weights/sum(weights)
}
