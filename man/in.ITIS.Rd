\name{in.ITIS}
\alias{in.ITIS}
\title{Check to see if taxon name is in ITIS}
\description{
  Checks to see if taxon names are listed in the ITIS taxon table.
}
\usage{
in.ITIS(df.parse, get.tax.env, col.sel = NULL)
}
\arguments{
  \item{df.parse}{Character vector or matrix character vectors of names
    that need to be checked for presence in ITIS}
  \item{get.tax.env}{Environment variable where the full ITIS taxon
    table is available}
  \item{col.sel}{Integer allowing selection of a particular column
    of the matrix of names to check.}
}
\value{
  Returns a list in which each element of the list is a successive
  character string from the parsed names.
}
\author{Lester L. Yuan}
\keyword{models}
