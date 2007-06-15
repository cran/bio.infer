\name{coef.west.wt}
\alias{coef.west.wt}
\docType{data}
\title{Weighted regression coefficients for western U.S.}
\description{
  Regression coefficients describing relationships between different
  taxa and percent sands/fines and stream temperature.  
}
\usage{data(coef.west.wt)}
\format{
  A list of with following four elements is provided.
  \describe{
    \item{\code{tnames}}{Character vector with the names of each taxon modelled.}
    \item{\code{csave}}{Matrix of regression coefficients for each taxon.  Number
    of rows is the same as the number of taxa modelled and the number of
    columns is 6.}
  \item{\code{xvar}}{Character vector with the names of variables modelled.}
  \item{\code{xlims}}{List of limits for each of the environmental
    variables.}
  \item{\code{form}}{Regression formula used
    to estimate coefficients.}
  }
}
\details{
  These regression coefficients were computed using logistic regression
  weighted by the number of stream miles represented by each sample.
}
\source{
  Yuan, L.L.
}
\examples{
data(coef.west.wt)
}
\keyword{datasets}
