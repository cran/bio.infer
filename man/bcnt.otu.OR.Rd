\name{bcnt.otu.OR}
\alias{bcnt.otu.OR}
\docType{data}
\title{Benthic count data with OTU}
\description{
Benthic count data from western Oregon with operational taxonomic units
assigned by \code{\link{get.otu}}.
}
\format{
  Five fields are included.
  \describe{
  \item{\code{SVN}}{Unique sample identifier}
  \item{\code{Taxon}}{Original taxon name, same as in \code{\link{bcnt.OR}}}
  \item{\code{CountValue}}{Taxon abundance}
  \item{\code{TNAME}}{Taxon name after processing by
    \code{\link{get.taxonomic}}}
  \item{\code{OTU}}{Operational taxonomic unit.  Taxon name selected with
    associated taxon-environment data and selected to eliminate
    double-counting}
}
}
  
\usage{data(bcnt.otu.OR)}
\source{
  Oregon Department of Environmental Quality
}
\examples{
data(bcnt.otu.OR)
}
\keyword{datasets}
