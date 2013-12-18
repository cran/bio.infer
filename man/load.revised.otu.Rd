\name{load.revised.otu}
\alias{load.revised.otu}

\title{Loads a user-revised OTU table}

\description{
  If manual correction of the OTU selections is desired, these
  corrections should be made on the "sum.otu.txt" file provided by
  \code{\link{get.otu}}.  The corrections should be saved as
  tab-delimited text, and then \code{load.revised.otu} run.   
}

\usage{
load.revised.otu(bcnt.otu, otufname = "sum.otu.txt")
}
\arguments{
  \item{bcnt.otu}{Original benthic count file with OTU designations.
    Output from \code{\link{get.otu}}
  }
  \item{otufname}{File name of revised OTU assignments.}
  }

\details{
  The benthic count data frame
  is assumed to be the original output from \code{\link{get.otu}}, and
  therefore has a field called \code{TNAME}. 
}

\value{
  The script returns a new benthic count file that incorporates the
  corrections to the OTU table. 
}
\seealso{\code{\link{get.otu}}}
\author{
Lester L. Yuan
}
\keyword{models}
