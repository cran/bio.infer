\name{get.valid.names}
\alias{get.valid.names}
\title{Get taxon names from benthic count file}
\description{
Identifies valid names from parsed name matrix.
}
\usage{
get.valid.names(df.parse, get.tax.env)
}
\arguments{
  \item{df.parse}{Parsed name matrix}
  \item{get.tax.env}{Environment for ITIS taxon table}
}
\details{
  Searches for taxon names in ITIS and splits name matrix into list with
  first element being the matrix of valid names and the second element
  being the matrix of invalid names.
}
\value{
  List of two taxon name matrices.
}
\author{Lester L. Yuan}
\keyword{models}
