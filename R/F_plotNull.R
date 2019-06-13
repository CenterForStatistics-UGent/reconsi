#' Plot the obtained null distribution along with a histogram of observed test
#' statistics
#' @param fit an object returned by the fdrCorrect() (or testDAA()) function
#' @param lowColor,highColor The low and high ends of the colour scale
#' @param dens a boolean, should fdr and Fdr be plotted?
#' @param idDA indices of known null taxa
#'
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
#' @examples
#'  p = 200; n = 50
#'  x = rep(c(0,1), each = n/2)
#'  mat = cbind(
#'  matrix(rnorm(n*p/10, mean = 5+x),n,p/10), #DA
#'  matrix(rnorm(n*p*9/10, mean = 5),n,p*9/10) #Non DA
#'  )
#Provide just the matrix and grouping factor, and test using the random null
#' fdrRes = fdrCorrect(mat, x)
#' plotNull(fdrRes)
plotNull = function(fit, lowColor ="yellow", highColor ="blue", dens = TRUE,
                    idDA = NULL){
    with(fit, {
    colnames(zValsDensPerm) = paste0("b", seq_len(ncol(zValsDensPerm)))
    df1 = data.frame(weight = weights, curve = colnames(zValsDensPerm))
    df2 = data.frame(zValsDensPerm, zSeq = zSeq)
    moltdf2 = melt(df2, id.vars = c("zSeq"), variable.name = "curve",
                   value.name = "Density")
    dfMerged = merge(moltdf2, df1, by = "curve")
    #Permutation densities
    plot = ggplot(data = dfMerged, aes(x =zSeq, group = curve, y = Density,
                                col = weight, alpha = weight)) +
        geom_line(linetype = "dashed", size = 0.5) +
        scale_colour_continuous(high = highColor,
                                low = lowColor, name = "Weights") +
        scale_alpha_continuous(guide = FALSE, range = c(0.5,1)) +
        xlab(if(zValues) "z-value" else "Test statistic") +
        ylab("Density/Fdr")  +
        theme_bw()

    #Histogram of observed z-values
    plot = plot + geom_histogram(data = data.frame(statObs = statObs),
                             aes(x = statObs, y = ..density..),
                             inherit.aes = FALSE, bins = 50, alpha = 0.5,
                             fill = "mediumseagreen")

    # Add density functions
    lfdr = g0/zValsDensObs*sum(zValsDensObs)/sum(g0)*p0
    lfdr[lfdr>1] = 1
    #Only show lfdr for observed z-values
    lfdr[zSeq > (max(statObs)+0.1) | zSeq < (min(statObs)-0.1)] = NA
    dfDens = data.frame(zSeq = zSeq, RandomNull = g0,
                        TheoreticalNull = dnorm(zSeq)*sum(g0)/sum(dnorm(zSeq)),
                        fdr = lfdr)
    if(!is.null(idDA)){
      #If null taxa known, add normal density
      nullZdens = estNormal(y = statObs[idDA], x = matrix(rep.int(1, sum(idDA))),
                            p = sum(idDA))
      dfDens$NullDensity = dnorm(zSeq, mean = nullZdens["mean.x1"], sd = nullZdens["sd"])
    }
    if(!dens){dfDens$fdr =NULL}
    dfDensMolt = melt(dfDens, id.vars ="zSeq", value.name = "density",
                      variable.name = "type")
    plot = plot + geom_line(inherit.aes = FALSE, data = dfDensMolt,
                  aes(x = zSeq, y = density, group = type, linetype = type)) +
        scale_linetype_manual(name = "", values = c("solid", "dashed", "dotdash"))
    if(dens){
    # Add red dots for Fdr estimates
    dfFdr = data.frame(statObs = statObs, Fdr = Fdr)
    plot = plot + geom_point(inherit.aes = FALSE, data = dfFdr,
                   aes(x = statObs, y = Fdr), col = "red", size = 0.75)
    }

    return(plot)
})
}
