#' A function to numerically stabilize an exponentiation. For internal use only
stabExp = function(exps){
exps-max(exps)
}