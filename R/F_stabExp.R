#' A function to numerically stabilize an exponentiation. For internal use only
#' @param exps the vector to be exponentiated
#'
#' @return the vector with the maximum subtracted
stabExp = function(exps){
    exps-max(exps)
}
