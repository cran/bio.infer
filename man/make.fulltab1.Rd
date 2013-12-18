\name{make.fulltab1}
\alias{make.fulltab1}
\title{Make full taxonomic table}
\description{
  Merge full ITIS taxonomy with names in benthic count file
}
\usage{
make.fulltab1(df.parse, get.tax.env)
}
\arguments{
  \item{df.parse}{Matrix with valid taxon names in the second column}
  \item{get.tax.env}{Environment variable where the full ITIS taxon
    table is available}
}
\details{
  Initial merge of taxon names with full ITIS table.
}
\value{
  Returns the current list of valid taxon names merged with the full
  taxonomic hierarchy.
}
\author{Lester L. Yuan}
\keyword{models}
