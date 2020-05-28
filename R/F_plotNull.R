#' Plot the obtained null distribution along with a histogram of observed test
#' statistics
#' @param fit an object returned by the reconsi() (or testDAA()) function
#' @param lowColor,highColor The low and high ends of the colour scale
#' @param idNull indices of known null taxa
#' @param nResampleCurves The number of resampling null distributions to plot
#' @param hSize A double, the size of the line of the collapsed null estimate
#'
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom stats density
#' @importFrom ks dkde
#' @export
#' @examples
#'  p = 180; n = 50; B = 1e2
#'  #Low number of resamples keeps computation time down
#'  x = rep(c(0,1), each = n/2)
#'  mat = cbind(
#'  matrix(rnorm(n*p/10, mean = 5+x),n,p/10), #DA
#'  matrix(rnorm(n*p*9/10, mean = 5),n,p*9/10) #Non DA
#'  )
#' #Provide just the matrix and grouping factor, and test using the random null
#' fdrRes = reconsi(mat, x, B = B)
#' plotNull(fdrRes)
plotNull = function(fit, lowColor ="yellow", highColor ="blue",
                    idNull = NULL, nResampleCurves = length(fit$Weights),
                    hSize = 0.5){
    with(fit, {
        zSeq = seq(min(statObs), max(statObs), length.out = 1e3)
        zValsDensPerm = if(resamAssumeNormal){
            apply(PermDensFits, 2, function(Fit){
             dnorm(zSeq, mean = Fit["mean"], sd = Fit["sd"])
            })
            } else {
                vapply(PermDensFits, FUN.VALUE = zSeq, function(Fit){dkde(zSeq, Fit)})
            }
    colnames(zValsDensPerm) = paste0("b", seq_len(ncol(zValsDensPerm)))
    idCurves = Weights >= sort(Weights, decreasing = TRUE)[nResampleCurves]
    #The Curves to plot
    df1 = data.frame(weight = Weights,
                     Curve = colnames(zValsDensPerm))[idCurves,]
    df2 = data.frame(zValsDensPerm[,idCurves], zSeq = zSeq)
    moltdf2 = melt(df2, id.vars = c("zSeq"), variable.name = "Curve",
                   value.name = "Density")
    dfMerged = merge(moltdf2, df1, by = "Curve")
    #Permutation densities
    plot = ggplot(data = dfMerged, aes(x = zSeq, group = Curve, y = Density,
                                col = weight, alpha = weight)) +
        geom_line(linetype = "dashed", size = 0.4) +
        scale_colour_continuous(high = highColor,
                                low = lowColor, name = "Weights") +
        scale_alpha_continuous(guide = FALSE, range = c(0.5,1)) +
        xlab(if(zValues) "z-value" else "Test statistic") +
        ylab("Density/Fdr")  + theme_bw()
    #Histogram of observed z-values
    plot = plot + geom_histogram(data = data.frame(statObs = statObs),
                             aes(x = statObs, y = ..density..),
                             inherit.aes = FALSE, bins = 50, alpha = 0.5,
                             fill = "mediumseagreen")
    # Add density functions
    g0 = if(assumeNormal) dnorm(zSeq,  mean = fitAll["mean"], sd = fitAll["sd"]) else dkde(zSeq, fitAll)
    lfdr = g0/dkde(zSeq, fitObs)*p0
    lfdr[lfdr>1] = 1
    #Only show lfdr for observed z-values
    lfdr[zSeq > (max(statObs)+0.1) | zSeq < (min(statObs)-0.1)] = NA
    dfDens = data.frame("zSeq" = zSeq, "ResampleNull" = g0,
                        "StandardNormal" = dnorm(zSeq)*sum(g0)/sum(dnorm(zSeq)),
                        "fdr" = lfdr)
    if(!is.null(idNull)){
      #If null taxa known, add true normal density
      nullZdens = estNormal(y = statObs[idNull])
      dfDens$OracleNullDensity = dnorm(zSeq, mean = nullZdens["mean"],
                                 sd = nullZdens["sd"])
    }
    dfDensMolt = melt(dfDens, id.vars = "zSeq", value.name = "density",
                      variable.name = "type")
    plot = plot + geom_line(inherit.aes = FALSE, data = dfDensMolt,
                  aes(x = zSeq, y = density, group = type,
                      linetype = type, size = type)) +
        scale_linetype_manual(name = "",
                              values = c("solid", "dashed", "dotdash", if(!is.null(idNull)) "twodash")) +
        scale_size_manual(values = c(0.2, 0.4, 0.4, if(!is.null(idNull)) 0.3), guide = FALSE)
    # Add red dots for Fdr estimates
    dfFdr = data.frame(statObs = statObs, Fdr = Fdr)
    plot = plot + geom_point(inherit.aes = FALSE, data = dfFdr,
                   aes(x = statObs, y = Fdr), col = "red", size = hSize)
    return(plot)
})
}
