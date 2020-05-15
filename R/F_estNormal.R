#' Fast estimation of mean and standard deviation of a normal distrbution,
#' optionally with weights
#' @param y vector of observations
#' @param w optional weight vector
#' @param p The number of features
#'
#' @return A vector of length 2 with mean and standard deviation
#' @importFrom stats sd
#' @importFrom matrixStats weightedMean
estNormal = function(y, w = NULL, p = length(y)){
    if(is.null(w)){
        c(mean = mean(y, na.rm = TRUE), sd = sd(y, na.rm = TRUE))
    } else {
        wmean = weightedMean(y, w = w, na.rm = TRUE)
        c(mean = wmean, sd = sqrt(sum((y-wmean)^2*w, na.rm = TRUE)*1/(p-1)))
    }
}
