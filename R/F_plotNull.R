#' Plot the obtained null distribution along with a histogram of observed test
#' statistics
#' @param fit an object returned by the fdrCorrect() (or testDAA()) function
#'
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
plotNull = function(fit){
    with(fit, {
    colnames(zValsDensPerm) = paste0("b", seq_len(ncol(zValsDensPerm)))
    df1 = data.frame(weight = weights, curve = colnames(zValsDensPerm))
    df2 = data.frame( zValsDensPerm, zSeq = zSeq)
    moltdf2 = melt(df2, id.vars = c("zSeq"), variable.name ="curve",
                   value.name ="Density")
    dfMerged = merge(moltdf2, df1, by = "curve")
    lfdr = g0/zValsDensObs*sum(zValsDensObs)/sum(g0)*p0
    lfdr[lfdr>1] = 1
    dfDens = data.frame(zSeq = zSeq, Observed = zValsDensObs, g0 = g0,
                        Theoretical = dnorm(zSeq)*sum(g0)/sum(dnorm(zSeq)),
                        fdr = lfdr)
    dfDensMolt = melt(dfDens, id.vars ="zSeq", value.name = "density",
                      variable.name = "type")
    dfFdr = data.frame(zValObs = zValObs, Fdr = Fdr)
    plot = ggplot(data = dfMerged, aes(x =zSeq, group = curve, y = Density,
                                col = weight, alpha = weight)) +
        geom_line(linetype = "dashed", size = 0.5) +
        scale_colour_continuous(high = "blue",
                                low = "yellow", name = "Weights") +
        scale_alpha_continuous(guide = FALSE, range = c(0.5,1)) +
        xlab("z-value") +
        ylab("Density/Fdr")  +
        theme_bw()

    # Add permutation densities
    plot = plot + geom_line(inherit.aes = FALSE, data = dfDensMolt,
                  aes(x = zSeq, y = density, group = type, linetype = type)) +
        scale_linetype_discrete(name = "Density type")

    # Add red dots for Fdr estimates
    plot = plot + geom_point(inherit.aes = FALSE, data = dfFdr,
                   aes(x = zValObs, y = Fdr), col = "red", size = 0.75)

    return(plot)
})
}
