#' Weighted kernel density estimation, optimized from mixtools package
#'
#' @param x The observations
#' @param u The point at which to evaluate the density
#' @param w the weights, same length as x
#' @param bw the bandwith
#'
#' @return densities of the same length as x
wkde = function (x, u = x, w, bw = bw.nrd0(x))
{
    vapply(u, FUN.VALUE = double(1), function(b) {
        sum(exp(-((x - b)/bw)^2/2)* w)
    })/(bw * sqrt(2 * pi))
}
