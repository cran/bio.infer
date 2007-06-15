\name{makess}
\alias{makess}

\title{Make a site-OTU matrix}

\description{
  \code{makess} generates a site-OTU matrix based on a benthic count
  file.  The script assumes that the benthic count file has a field
  entitled 'OTU' on which the matrix is based. 
}

\usage{
makess(bcnt, tname = "OTU", plothist=FALSE, prints = FALSE, nview = 0)
}

\arguments{

  \item{bcnt}{a benthic count file that results from running
    \code{\link{get.otu}} or \code{\link{load.revised.otu}}.
  }
  \item{tname}{Character string specifying the field in bcnt in which
    the name of the taxon is stored.}
  \item{plothist}{Set to \code{TRUE} to plot a histogram of the proportion of
    abundance at each site that is included in the OTU designations}
  \item{prints}{Set to \code{TRUE} to view statistics on the proportion
    of abundance at each site that is included in the OTU designations}
  \item{nview}{Selects the number of sites to view in which the
    proportion of abundance included is low}
}

\value{
  A site-OTU matrix.
}
\examples{
data(bcnt.otu.OR)
ss.OR <- makess(bcnt.otu.OR)
}
\seealso{\code{\link{get.otu}}, \code{\link{load.revised.otu}}}
\author{
Lester L. Yuan
}
\keyword{models}
