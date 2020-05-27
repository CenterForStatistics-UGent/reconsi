#' Prepare matrix from eighted kernel density estimation,
#' optimized from mixtools package
#'
#' @param x The observations
#' @param u The point at which to evaluate the density
#' @param bw the bandwith
#'
#' @return densities of the same length as x
wkdePrep = function (x, u = x, bw)
{
    vapply(u, FUN.VALUE = x, function(b) {
        exp(-((x - b)/bw)^2/2)
    })
}
