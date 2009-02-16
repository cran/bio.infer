\name{parse.taxon.name}
\alias{parse.taxon.name}
\title{Parse taxon name}
\description{
  Parses taxon names from benthic count file to facilitate matching with
  ITIS.
}
\usage{
parse.taxon.name(tname.orig)
}
\arguments{
  \item{tname.orig}{Character vector of taxon names}
}
\details{
  Parses compound names typically found in benthic count files into
  a list of strings for each name.
}
\value{
  Returns a matrix in which each row corresponds with a different
  element of the provided character vector, and each column corresponds
  to distinct strings in each element.
}
\author{Lester L. Yuan}
\keyword{models}
