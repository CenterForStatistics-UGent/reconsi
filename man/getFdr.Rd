% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/F_getFdr.R
\name{getFdr}
\alias{getFdr}
\title{Calculate tail-area (Fdr) and local (fdr) false discovery rates,
   based on a certain null distribution}
\usage{
getFdr(
  statObs,
  fitAll,
  fdr,
  p,
  p0,
  zValsDensObs,
  smoothObs,
  assumeNormal,
  fitObs,
  ...
)
}
\arguments{
\item{statObs}{Vector of observed z-values}

\item{fitAll}{The parameters of the estimated random null}

\item{fdr}{local false discovery rate, already estimated}

\item{p}{the number of hypotheses}

\item{p0}{The estimated fraction of null hypotheses}

\item{zValsDensObs}{estimated densities of observed test statistics}

\item{smoothObs}{A boolean, should estimated observed densities of the test
statistics be used in estimating the Fdr}

\item{assumeNormal}{A boolean, should normality be assumed for the null distribution?}

\item{fitObs}{The kernel density estimate object of all test statistics}

\item{...}{more arguments, ignored}
}
\value{
A list with components
\item{Fdr}{Tail are false discovery rate}
\item{fdr}{Local false discovery rate}
}
\description{
Calculate tail-area (Fdr) and local (fdr) false discovery rates,
   based on a certain null distribution
}
