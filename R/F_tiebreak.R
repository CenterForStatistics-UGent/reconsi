#' A function to break the ties because of zeroes in the sequencing counts, and prepare the final matrix
#' @param Y the sequence count matrix
#' @param FC the vector of total cell counts
#' @param S the library sizes. By default calculated from Y
#' @param tieBreak the tiebreaking protocol
#'
tiebreak = function(Y, FC, S = rowSums(Y), tieBreak ="none", tol = 1e-4, maxIter = 1L){
  if(tieBreak=="none"){
    anaMat = Y/S*FC
  } else {
    id = Y==0
    if(tieBreak=="libSizes"){
    anaMat = Y/S*FC
    anaMat[id] = matrix(1/S, nrow(Y), ncol(Y))[id]-1
    } else if(tieBreak=="Binomial"){
    probMat = Y!=0
    captProb = colMeans(probMat)
    iter = 1;convergence  = FALSE
    while(iter<=maxIter && !convergence){
      captProbOld = captProb
      gammas = outer(S, colSums(Y)/sum(Y)*(1-captProb))
      probMat[id] = 1-exp(-gammas)[id]
      captProb = colMeans(probMat)
      iter = iter+1
      convergence = sqrt(mean((captProb-captProbOld)^2)) < tol
    }
    anaMat = Y/S*FC
    anaMat[id] = probMat[id]-1
    } else if(tieBreak=="Bernoulli"){
    px = exp(-outer(S, colSums(Y)/sum(Y)))
    probMat = !id; presProb = colMeans(probMat)
    iter = 1;convergence  = FALSE
    while(iter<=maxIter && !convergence){
      presProbOld = presProb
      pz = matrix(presProb, byrow = TRUE, nrow(Y), ncol(Y))
      probMat[id]=  (px*pz/(pz*(px-1)+1))[id]
      presProb = colMeans(probMat)
      iter = iter+1
      convergence = sqrt(mean((presProb-presProbOld)^2)) < tol
    }
    anaMat = Y/S*FC
    anaMat[id] = probMat[id]-1
    } else if (tieBreak=="FC"){
    anaMat = Y/S*FC
    anaMat[id] = -1/matrix(FC, nrow(Y), ncol(Y))[id]
    } else {stop("Tiebreaking paradigm not known! \n")}
  }
  return(anaMat)
}
