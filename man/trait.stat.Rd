\name{trait.stat}
\alias{trait.stat}
\title{Compute trait metrics}
\description{
Compute trait metrics from benthic count data.
}
\usage{
trait.stat(bcnt.otu, coefs)
}
\arguments{
  \item{bcnt.otu}{Benthic count file with OTU defined consistently with selected trait file}
  \item{coefs}{Trait file}
}
\details{
  OTUs appropriate for the selected trait file must first be assigned
  using \code{\link{get.otu}}.  A list of trait files provided in the
  library can be obtained using \code{\link{sel.coeffile}}.  Script will
  compute richness, relative richness, and relative abundance of each
  trait group.  Note that many taxa in a benthic count file are often
  not assigned to traits, and so relative richness and relative
  abundance values within a given sample may not sum to 1.
}
\value{
  Data frame of trait metrics.
}
\author{Lester L. Yuan}
\keyword{models}
