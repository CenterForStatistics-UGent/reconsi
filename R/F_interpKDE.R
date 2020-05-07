#' Evaluate a kernel density through interpolation
#' @param kdeObj the kernel density fit
#' @param newData A vector of values at which there is to be interpolated
#'
#' @return the interpolated value vector of the length of newData
interpKDE = function(kdeObj, newData){
    app = approx(y = kdeObj$y, x = kdeObj$x, xout = newData)$y
    app[app<=0] = .Machine$double.eps
    app
}
#'Stabilize the fitting of bkde
#'@param ... passed on to bkde
#'@return see ?bkde
#' @importFrom KernSmooth bkde
bkdeStab = function(...){
    obj = bkde(...)
    obj$y[obj$y<=0] = .Machine$double.eps
    obj
}
