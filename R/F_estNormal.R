#' Fast estimation of mean and standard deviation of a normal distrbution,
#' optionally with weights
#' @param y vector of observations
#' @param w optional weight vector
#' @param p The number of features
#'
#' @return A vector of length 2 with mean and standard deviation
#' @importFrom stats sd weighted.mean
estNormal = function(y, w = NULL, p){
    if(is.null(w)){
        c(mean = mean(y), sd = sd(y))
    } else {
        wmean = weighted.mean(y, w = w)
        c(mean = wmean, sd = sqrt(sum((y-wmean)^2*w)*1/(p-1)))
    }
}
