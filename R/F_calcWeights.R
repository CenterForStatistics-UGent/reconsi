#' Obtain weights as posterior probabilities to calculate the consensus null
#'
#' @param logDensPerm A matrix with B rows of logged density estimates of the B
#'    permutation distributions
#' @param fdr A vector of local false discovery rates
#'    along the support of the kernel
#' @param logPriorProbs A vector of length B with logged prior probabilities of
#' resampling instances
#'
#' @return A vector of weights of length B
calcWeights = function(logDensPerm, fdr, logPriorProbs){
  logVec = colSums(logDensPerm*fdr)
  if(!is.null(logPriorProbs)) logVec = logVec + logPriorProbs
  weights = exp(stabExp(logVec))
  weights/sum(weights)
}
