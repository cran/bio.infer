\name{view.te}
\alias{view.te}
\title{View taxon-environment relationships}
\description{
  Viewer for taxon-environment relationships.
}
\usage{
view.te(coef, plotform = "pdf")
}
\arguments{
  \item{coef}{Taxon-environment coefficient file
    (e.g. \code{\link{coef.west.wt}})}
  \item{plotform}{Specify output plot format.  Options include pdf and
    windows.}
  }
\details{
  Generates contour plots that show taxon-environment relationships
  represented by taxon-environment coefficients.  Models with a single
  explanaotry variable are displayed as line plots, and models with two
  explanatory variables are displayed as contour plots.
  Plots are returned returned the file \code{taxon.env.pdf} is
  \code{plotform = "pdf"}.
}
\value{
  None.  Used for side effects.
}
\author{Lester L. Yuan}
\examples{
\dontrun{
data(coef.west.wt)
view.te(coef.west.wt)}
}
\keyword{models}

